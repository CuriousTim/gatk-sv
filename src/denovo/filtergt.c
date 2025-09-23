#define _POSIX_C_SOURCE 202405L // getline, strdup
#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <string.h>

#include <htslib/vcf.h>
#include <htslib/hts.h>
#include <htslib/khash.h>
#include <htslib/kstring.h>

#define err(fmt, ...)                                              \
	do {                                                       \
		fprintf(stderr, "error: " fmt "\n", __VA_ARGS__);  \
		exit(EXIT_FAILURE);                                \
	} while (0)

KHASH_SET_INIT_STR(set)
KHASH_MAP_INIT_STR(map, khash_t(set)*)

// pool of truth variant IDs
khash_t(set) *tvids_pool;
// pool of sample IDs of samples with alt allele at variant site
khash_t(set) *samples_pool;
// map between truth variants and samples with alt allele at variant site
khash_t(map) *genotypes;
// map between eval variants and truth variants
khash_t(map) *concordance;

const char *pool_get_string(khash_t(set) *h, char *s);
khiter_t set_insert(khash_t(set) *h, const char *s);

void usage(FILE *fp)
{
	fprintf(fp, "usage: filtergt <inbcf> <concordance_vcf> <genotypes> <outbcf> <filter_flag>\n");
}

/**
 * Split the string in `s` by `sep`, save a copy of each field in `pool` if
 * it doesn't exist there and insert the pointer from `pool` into `set`.
 */
void split_and_insert(char *s, khash_t(set) *pool, khash_t(set) *set,
		char sep)
{
	char *p = s;
	const char *f;
	while ((p = strchr(s, sep))) {
		*p = '\0';
		f = pool_get_string(pool, s);
		set_insert(set, f);
		s = p + 1;
	}
	f = pool_get_string(pool, s);
	set_insert(set, f);
}

/**
 * Insert a string into a set, crashing if the operations fails. Return an
 * iterator to the inserted element.
 */
khiter_t set_insert(khash_t(set) *h, const char *s)
{
	int ret;
	khiter_t tmp = kh_put(set, h, s, &ret);
	if (ret < 0)
		err("%s", "failed to insert into set");

	return tmp;
}

/**
 * Insert a string as a key into a map, crashing if the operations fails.
 * Return an iterator to the inserted element.
 */
khiter_t map_insert_key(khash_t(map) *h, const char *s)
{
	int ret;
	khiter_t tmp = kh_put(map, h, s, &ret);
	if (ret < 0)
		err("%s", "failed to insert key into map");

	return tmp;
}

/**
 * Get a string from a string pool, copying and inserting the string first if
 * it doesn't exist.
 */ 
const char *pool_get_string(khash_t(set) *h, char *s)
{
	khiter_t i = kh_get(set, h, s);
	if (i == kh_end(h)) {
		const char *tmp = strdup(s);
		if (!tmp)
			err("%s", "OOM");
		i = set_insert(h, tmp);
	}

	return kh_key(h, i);
}

/**
 * Add `key` to the hash map `h` and add the values in `vals` to the hash
 * set associated with the key. `vals` is a string of fields separated by
 * `sep` and `pool` is the string pool to use for the values. If `cp` is
 * `true`, then `key` will be copied before inserting it into the map,
 * otherwise the string pointer will be passed directly to the map.
 */
void map_insert(khash_t(map) *h, const char *key, char *vals,
		khash_t(set) *pool, char sep, bool cp)
{
	khiter_t p = kh_get(map, h, key);
	khash_t(set) *s;
	if (p == kh_end(h)) {
		const char *tmp = key;
		if (cp) {
			tmp = strdup(key);
			if (!tmp)
				err("%s", "OOM");
		}
		p = map_insert_key(h, tmp);
		s = kh_init(set);
		if (!s)
			err("%s", "OOM");
		kh_val(h, p) = s;
	}

	s = kh_val(h, p);
	split_and_insert(vals, pool, s, sep);
}

/**
 * Add a record from the concordance VCF to the global hash map.
 */
void add_concordance_record(bcf_hdr_t *hdr, bcf1_t *rec)
{
	// assume record was unpacked
	char *vid = rec->d.id;
	// ID is missing so ignore
	if (!vid || (vid[0] == '.' && vid[1] == '\0'))
		return;

	char *dst = 0;
	int ndst = 0;
	int ret = bcf_get_info_string(hdr, rec, "TRUTH_VID", &dst, &ndst);
	switch(ret) {
	case 0:
		return;
	case -1:
		err("%s", "TRUTH_VID tag is missing from header");
		break;
	case -2:
		err("%s", "mismatch TRUTH_VID types between header and record");
		break;
	case -3: // tag is not present in the VCF record
		return;
	case -4:
		err("%s", "OOM");
		break;
	}

	map_insert(concordance, vid, dst, tvids_pool, ',', true);
	hts_free(dst);
}

/**
 * Load all the TRUTH_VID records from the concordance VCF.
 */
void load_concordance(const char *path)
{
	htsFile *fp = hts_open(path, "r");
	if (!fp)
		err("%s", "could not open concordance VCF");
	
	bcf_hdr_t *hdr = bcf_hdr_read(fp);
	if (!hdr)
		err("%s", "could not read concordance VCF header");

	bcf1_t *rec = bcf_init();
	while (bcf_read(fp, hdr, rec) == 0) {
		bcf_unpack(rec, BCF_UN_INFO);
		add_concordance_record(hdr, rec);
	}
	bcf_destroy(rec);
	bcf_hdr_destroy(hdr);
	hts_close(fp);
}

/**
 * Add a record from the genotypes file to the global hash map.
 */
void add_genotype_record(char *line)
{
	char *q = strchr(line, '\t');
	// first field is variant ID, rest are sample IDs
	if (!q)
		return;

	*q = '\0';
	khiter_t k = kh_get(set, tvids_pool, line);
	// only add genotypes for truth variants from concordance
	if (k == kh_end(tvids_pool))
		return;

	map_insert(genotypes, kh_key(tvids_pool, k), q + 1, samples_pool, '\t', false);
}

/**
 * Load all the records from the genotypes file.
 */
void load_genotypes(const char *path)
{
	htsFile *fp = hts_open(path, "r");
	if (!fp)
		err("%s", "failed to open genotypes file");

	kstring_t s = KS_INITIALIZE;
	while (hts_getline(fp, '\n', &s) >= 0) {
		add_genotype_record(s.s);
	}
	ks_free(&s);
	hts_close(fp);
}

/**
 * Does the sample `sid` have an alt allele at the truth variant `vid`?
 */
bool sample_alt_in_truth(const char *vid, const char *sid)
{
	khiter_t p = kh_get(map, genotypes, vid);
	if (p == kh_end(genotypes))
		return false;

	khash_t(set) *samples = kh_val(genotypes, p);

	return kh_get(set, samples, sid) != kh_end(samples);
}

/**
 * Does the sample `sid` have support for an alt genotype at the eval variant
 * `vid`?
 */
bool gt_supported_in_truth(const char *vid, const char *sid)
{
	khiter_t p = kh_get(map, concordance, vid);
	if (p == kh_end(concordance))
		return false;

	khash_t(set) *tvids = kh_val(concordance, p);
	for (khiter_t i = kh_begin(tvids); i != kh_end(tvids); ++i) {
		if (kh_exist(tvids, i) && sample_alt_in_truth(kh_key(tvids, i), sid))
			return true;
	}

	return false;
}

/**
 * Update genotypes in a BCF record according to concordance.
 *
 * @param hdr                  BCF header.
 * @param rec                  BCF record.
 * @param null_unsupported_gt  If `true`, unsupported genotypes will be set to
 *   null. If `false`, supported genotypes will be set to null. The point of
 *   the `false` case to is nullify genotypes in offspring that have support in
 *   the parents when calling de novo variants.
 */
void update_genotypes(bcf_hdr_t *hdr, bcf1_t *rec, bool null_unsupported_gt)
{
	const char *vid = rec->d.id;
	int nsample = bcf_hdr_nsamples(hdr);
	int32_t *gt_arr = 0;
	int ngt_arr = 0;
	int ngt = bcf_get_genotypes(hdr, rec, &gt_arr, &ngt_arr);
	int ploidy = ngt / nsample;
	if (ploidy != 2)
		err("%s", "sample ploidy is not 2");

	for (int i = 0; i < nsample; ++i) {
		int32_t *p = gt_arr + i * ploidy;
		// we want to ignore genotypes that don't have alt allele
		// assume no genotypes are 1/. or ./1
		if (bcf_gt_is_missing(p[0])
				|| bcf_gt_is_missing(p[1])
				|| (bcf_gt_allele(p[0]) == 0 && bcf_gt_allele(p[1]) == 0))
			continue;
		bool gt_supported = gt_supported_in_truth(vid, hdr->id[BCF_DT_SAMPLE][i].key);
		if (gt_supported != null_unsupported_gt) {
			p[0] = bcf_gt_missing;
			p[1] = bcf_gt_missing;
		}
	}
	bcf_update_genotypes(hdr, rec, gt_arr, ngt);
	hts_free(gt_arr);
}

bool is_cnv(const bcf_hdr_t *hdr, bcf1_t *rec)
{
	bcf_unpack(rec, BCF_UN_INFO);
	char *svtype = 0;
	int n = 0;
	int ret = bcf_get_info_values(hdr, rec, "SVTYPE", (void **)&svtype, &n, BCF_HT_STR);
	if (ret < 0)
		err("%s", "record is missing an SVTYPE");

	bool tmp = strcmp(svtype, "CNV") == 0;
	hts_free(svtype);

	return tmp;
}

void free_hashes(void)
{
	for (khiter_t i = kh_begin(concordance); i != kh_end(concordance); ++i) {
		if (kh_exist(concordance, i)) {
			free((void *)kh_key(concordance, i));
			kh_destroy(set, kh_val(concordance, i));
		}
	}
	kh_destroy(map, concordance);

	for (khiter_t i = kh_begin(genotypes); i != kh_end(genotypes); ++i) {
		if (kh_exist(genotypes, i))
			kh_destroy(set, kh_val(genotypes, i));
	}
	kh_destroy(map, genotypes);

	for (khiter_t i = kh_begin(tvids_pool); i != kh_end(tvids_pool); ++i) {
		if (kh_exist(tvids_pool, i))
			free((void *)kh_key(tvids_pool, i));
	}
	kh_destroy(set, tvids_pool);


	for (khiter_t i = kh_begin(samples_pool); i != kh_end(samples_pool); ++i) {
		if (kh_exist(samples_pool, i))
			free((void *)kh_key(samples_pool, i));
	}
	kh_destroy(set, samples_pool);
}

int main(int argc, char *argv[])
{
	if (argc != 6) {
		usage(stderr);
		return EXIT_FAILURE;
	}

	bool null_unsupported_gt = argv[5][0] == '1';

	tvids_pool = kh_init(set);
	samples_pool = kh_init(set);
	genotypes = kh_init(map);
	concordance = kh_init(map);

	load_concordance(argv[2]);
	load_genotypes(argv[3]);

	htsFile *infp = hts_open(argv[1], "r");
	if (!infp)
		err("%s", "failed to open input BCF");

	htsFile *outfp = hts_open(argv[4], "wb");
	if (!outfp)
		err("%s", "failed to open output BCF");

	bcf_hdr_t *hdr = bcf_hdr_read(infp);
	if (bcf_hdr_write(outfp, hdr) != 0)
		err("%s", "failed to write header to output BCF");

	bcf1_t *rec = bcf_init();
	while (bcf_read(infp, hdr, rec) == 0) {
		if (is_cnv(hdr, rec))
			continue;
		bcf_unpack(rec, BCF_UN_ALL);
		update_genotypes(hdr, rec, null_unsupported_gt);
		if (bcf_write(outfp, hdr, rec) != 0)
			err("%s", "failed to write record to output BCF");
	}

	bcf_destroy(rec);
	bcf_hdr_destroy(hdr);
	hts_close(infp);
	hts_close(outfp);
	free_hashes();
	
	return EXIT_SUCCESS;
}
