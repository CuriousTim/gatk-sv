version 1.0

import "Structs.wdl"
import "DeNovoSvsProcessOffspringVcf.wdl"
import "DeNovoSvsSvConcordancePerContig.wdl"

###########################
# MAIN WORKFLOW DEFINITION
###########################

workflow DeNovoSvs {
  input {
    File pedigree
    # One family ID per line to call de novo in subset of families
    File? family_ids

    Float? max_cohort_af
    Float? max_gnomad_af
    Int? large_cnv_size
    Int? depth_only_size
    Array[File]? exclude_regions
    Float? exclude_regions_ovp
    File gd_regions
    Float? gd_regions_ovp

    # Either a single VCF or an array of VCFs with each one containing a single
    # contig. In the case of a single VCF, it is expected that all the contigs
    # in the input contigs are present. In the case of multiple VCFs, all VCFs
    # must contain the exact same set of samples.
    Array[File]+ vcfs
    Array[File]+ vcf_indices

    Array[String]+ contigs = ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6",
      "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15",
      "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX"]

    File reference_dict

    # Raw data
    Array[String] batch_name_list            # batch IDs
    Array[File] batch_sample_lists           # samples in each batch (filtered set)
    Array[String] batch_bincov_matrix
    Array[String] batch_bincov_matrix_index
    Array[String]? clustered_manta_vcf
    Array[String]? clustered_melt_vcf
    Array[String]? clustered_wham_vcf
    Array[String]? clustered_scramble_vcf
    Array[String] clustered_depth_vcf

    # Dockers
    String linux_docker
    String sv_base_mini_docker
    String svconcordance_keep_all_docker
    String denovo_docker

    RuntimeAttr? runtime_override_make_manifests
    RuntimeAttr? runtime_override_subset_vcf_by_contig
    RuntimeAttr? runtime_override_subset_samples
    RuntimeAttr? runtime_override_match_vcf_to_contig
    RuntimeAttr? runtime_override_subset_bcf_by_samples
    RuntimeAttr? runtime_override_remove_uncalled_svtypes
    RuntimeAttr? runtime_override_filter_offspring_sites
    RuntimeAttr? runtime_override_group_offspring_by_batch
    RuntimeAttr? runtime_override_concat_raw_evidence
    RuntimeAttr? runtime_override_svconcordance
  }

  output {
    Array[Array[File]] offspring_concordance = select_all(offspring_v_offspring.concordance_vcfs)
    Array[Array[File]] father_concordance = select_all(offspring_v_father.concordance_vcfs)
    Array[Array[File]] mother_concordance = select_all(offspring_v_mother.concordance_vcfs)
  }

  call MakeManifests {
    input:
      batch_name_list = batch_name_list,
      batch_sample_lists = batch_sample_lists,
      batch_bincov_matrix = batch_bincov_matrix,
      batch_bincov_matrix_index = batch_bincov_matrix_index,
      clustered_manta_vcf = clustered_manta_vcf,
      clustered_melt_vcf = clustered_melt_vcf,
      clustered_wham_vcf = clustered_wham_vcf,
      clustered_scramble_vcf = clustered_scramble_vcf,
      clustered_depth_vcf = clustered_depth_vcf,
      denovo_docker = denovo_docker,
      runtime_attr_override = runtime_override_make_manifests
  }

  if (length(vcfs) == 1) {
    scatter (i in range(length(contigs))) {
      call SubsetVcfByContig {
        input:
          vcf = vcfs[0],
          vcf_index = vcf_indices[0],
          contig = contigs[i],
          sv_base_mini_docker = sv_base_mini_docker,
          runtime_attr_override = runtime_override_subset_vcf_by_contig
      }
    }
    Array[File] subset_vcfs = select_all(SubsetVcfByContig.subset_vcf)
    Array[File] subset_vcf_indices = select_all(SubsetVcfByContig.subset_vcf_index)
  }

  call SubsetSamples {
    input:
      ped = pedigree,
      fams = family_ids,
      vcf = vcfs[(length(vcfs) - 1)],
      sample_manifest = MakeManifests.sample_manifest,
      sv_base_mini_docker = sv_base_mini_docker,
      runtime_attr_override = runtime_override_subset_samples
  }

  Array[File] split_vcfs = select_first([subset_vcfs, vcfs])
  Array[File] split_vcf_indices = select_first([subset_vcf_indices, vcf_indices])
  scatter (i in range(length(split_vcfs))) {
    call MatchVcfToContig {
      input:
        vcf = split_vcfs[i],
        vcf_index = split_vcf_indices[i],
        contigs = contigs,
        sv_base_mini_docker = sv_base_mini_docker,
        runtime_attr_override = runtime_override_match_vcf_to_contig
    }
  }

  Array[File] matched_vcfs = select_all(MatchVcfToContig.matched_vcf)
  Array[String] kept_contigs = select_all(MatchVcfToContig.matched_contig)
  scatter (i in range(length(matched_vcfs))) {
    call DeNovoSvsProcessOffspringVcf.DeNovoSvsProcessOffspringVcf as process_offspring_vcf {
      input:
        vcf = matched_vcfs[i],
        offspring = SubsetSamples.offspring,
        max_cohort_af = max_cohort_af,
        max_gnomad_af = max_gnomad_af,
        large_cnv_size = large_cnv_size,
        depth_only_size = depth_only_size,
        exclude_regions = exclude_regions,
        exclude_regions_ovp = exclude_regions_ovp,
        gd_regions = gd_regions,
        gd_regions_ovp = gd_regions_ovp,
        batches = SubsetSamples.batch_subset,
        pedigree = SubsetSamples.ped_subset,
        sample_manifest = MakeManifests.sample_manifest,
        linux_docker = linux_docker,
        sv_base_mini_docker = sv_base_mini_docker,
        runtime_override_subset_bcf_by_samples = runtime_override_subset_bcf_by_samples,
        runtime_override_remove_uncalled_svtypes = runtime_override_remove_uncalled_svtypes,
        runtime_override_filter_offspring_sites = runtime_override_filter_offspring_sites,
        runtime_override_group_offspring_by_batch = runtime_override_group_offspring_by_batch
    }
  }

  # [
  #   [ contig0-batch0, contig0-batch1, contig0-batch2, ...],
  #   [ contig1-batch0, contig1-batch1, contig1-batch2, ...],   <--+
  #   ...                                                           \
  # ]                                                               |
  #                                                           transpose
  # [                                                               |
  #   [ contig0-batch0, contig1-batch0, contig2-batch0, ...],       /
  #   [ contig0-batch1, contig1-batch1, contig2-batch1, ...],   <--+
  #   ...
  # ]
  Array[Array[File]] offspring_batch_grouped_bcfs = transpose(process_offspring_vcf.offspring_batch_grouped_bcfs)
  Array[Array[File]] father_batch_grouped_bcfs = transpose(process_offspring_vcf.father_batch_grouped_bcfs)
  Array[Array[File]] mother_batch_grouped_bcfs = transpose(process_offspring_vcf.mother_batch_grouped_bcfs)
  Array[Array[File]] offspring_batch_offspring_ids = transpose(process_offspring_vcf.offspring_batch_offspring_ids)
  Array[Array[File]] father_batch_father_ids = transpose(process_offspring_vcf.father_batch_father_ids)
  Array[Array[File]] mother_batch_mother_ids = transpose(process_offspring_vcf.mother_batch_mother_ids)

  Array[String] kept_batches = read_lines(SubsetSamples.batch_subset)
  scatter (i in range(length(kept_batches))) {
    String current_batch = kept_batches[i]
    if (size(offspring_batch_grouped_bcfs[i]) > 0) {
      call DeNovoSvsSvConcordancePerContig.DeNovoSvsSvConcordancePerContig as offspring_v_offspring {
        input:
          truth_vcfs = MakeManifests.raw_vcf_map[current_batch],
          truth_samples = offspring_batch_offspring_ids[i][0],
          eval_bcfs = offspring_batch_grouped_bcfs[i],
          contigs = kept_contigs,
          reference_dict = reference_dict,
          batch = current_batch,
          sv_base_mini_docker = sv_base_mini_docker,
          svconcordance_keep_all_docker = svconcordance_keep_all_docker,
          linux_docker = linux_docker,
          runtime_override_concat_raw_evidence = runtime_override_concat_raw_evidence,
          runtime_override_svconcordance = runtime_override_svconcordance
      }
    }

    if (size(father_batch_grouped_bcfs[i]) > 0) {
      call DeNovoSvsSvConcordancePerContig.DeNovoSvsSvConcordancePerContig as offspring_v_father {
        input:
          truth_vcfs = MakeManifests.raw_vcf_map[current_batch],
          truth_samples = father_batch_father_ids[i][0],
          eval_bcfs = father_batch_grouped_bcfs[i],
          contigs = kept_contigs,
          reference_dict = reference_dict,
          batch = current_batch,
          sv_base_mini_docker = sv_base_mini_docker,
          svconcordance_keep_all_docker = svconcordance_keep_all_docker,
          linux_docker = linux_docker,
          runtime_override_concat_raw_evidence = runtime_override_concat_raw_evidence,
          runtime_override_svconcordance = runtime_override_svconcordance
      }
    }

    if (size(mother_batch_grouped_bcfs[i]) > 0) {
      call DeNovoSvsSvConcordancePerContig.DeNovoSvsSvConcordancePerContig as offspring_v_mother {
        input:
          truth_vcfs = MakeManifests.raw_vcf_map[current_batch],
          truth_samples = mother_batch_mother_ids[i][0],
          eval_bcfs = mother_batch_grouped_bcfs[i],
          contigs = kept_contigs,
          reference_dict = reference_dict,
          batch = current_batch,
          sv_base_mini_docker = sv_base_mini_docker,
          svconcordance_keep_all_docker = svconcordance_keep_all_docker,
          linux_docker = linux_docker,
          runtime_override_concat_raw_evidence = runtime_override_concat_raw_evidence,
          runtime_override_svconcordance = runtime_override_svconcordance
      }
    }
  }
}

# Create manifests of the paths to the raw evidence files.
task MakeManifests {
  input {
    Array[String] batch_name_list
    Array[File] batch_sample_lists
    # These are not really optional. They should correspond to the raw algorithms
    # used to generate the callset used to run the workflow, but making them optional
    # is the only way to allow for different algorithms in different versions of
    # GATK-SV.
    Array[String]? clustered_manta_vcf
    Array[String]? clustered_melt_vcf
    Array[String]? clustered_wham_vcf
    Array[String]? clustered_scramble_vcf
    Array[String] clustered_depth_vcf
    Array[String] batch_bincov_matrix
    Array[String] batch_bincov_matrix_index
    String denovo_docker
    RuntimeAttr? runtime_attr_override
  }

  parameter_meta {
    batch_name_list: "GATK-SV batch IDs."
    batch_sample_lists: "One file per batch listing the samples, one per line, in that batch."
    clustered_manta_vcf: "The paths to the clustered Manta VCFs for each batch."
    clustered_melt_vcf: "The paths to the clustered MELT VCFs for each batch."
    clustered_wham_vcf: "The paths to the clustered Wham VCFs for each batch."
    clustered_scramble_vcf: "The paths to the clustered Scramble VCFs for each batch."
    clustered_depth_vcf: "The paths to the clustered depth VCFs for each batch."
    batch_bincov_matrix: "The paths to the merged bincov matrix for each batch."
    batch_bincov_matrix_index: "The paths to the merged bincov matrix index for each batch."
    denovo_docker: "The corresponding Docker image from GATK-SV."
    runtime_attr_override: "Runtime attribute overrides."
  }

  output {
    File sample_manifest = "sample_manifest.tsv"
    Map[String, Array[String]] raw_vcf_map = read_json("raw_manifest.json")["raw_vcf"]
    Map[String, String] bincov_map = read_json("bincov_manifest.json")["bincov"]
    Map[String, String] bincov_index_map = read_json("bincov_index_manifest.json")["bincov_index"]
  }

  Float input_size = size(batch_sample_lists, "GB")
  RuntimeAttr default_attr = object {
    mem_gb: 1,
    cpu_cores: 1,
    disk_gb: ceil(input_size * 3) + 16,
    boot_disk_gb: 8,
    preemptible_tries: 3,
    max_retries: 1,
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  Array[String] manta_vcfs = select_first([clustered_manta_vcf, []])
  Array[String] melt_vcfs = select_first([clustered_melt_vcf, []])
  Array[String] wham_vcfs = select_first([clustered_wham_vcf, []])
  Array[String] scramble_vcfs = select_first([clustered_scramble_vcf, []])

  runtime {
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GB"
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    docker: denovo_docker
  }

  command <<<
    set -euxo pipefail

    batch_names='~{write_lines(batch_name_list)}'

    paste "${batch_names}" '~{write_lines(clustered_depth_vcf)}' > 'depth_manifest.tsv'

    manta='~{if length(manta_vcfs) > 0 then write_lines(manta_vcfs) else ""}'
    melt='~{if length(melt_vcfs) > 0 then write_lines(melt_vcfs) else ""}'
    wham='~{if length(wham_vcfs) > 0 then write_lines(wham_vcfs) else ""}'
    scramble='~{if length(scramble_vcfs) > 0 then write_lines(scramble_vcfs) else ""}'

    : > pesr_manifest.tsv
    if [[ "${manta}" ]]; then
        paste "${batch_names}" "${manta}" >> pesr_manifest.tsv
    fi
    if [[ "${melt}" ]]; then
        paste "${batch_names}" "${melt}" >> pesr_manifest.tsv
    fi
    if [[ "${wham}" ]]; then
        paste "${batch_names}" "${wham}" >> pesr_manifest.tsv
    fi
    if [[ "${scramble}" ]]; then
        paste "${batch_names}" "${scramble}" >> pesr_manifest.tsv
    fi

    if [[ ! -s pesr_manifest.tsv ]]; then
        printf 'at least one non-empty list of PESR evidence VCFs should be provided\n' >&2
        exit 1
    fi

    cat depth_manifest.tsv pesr_manifest.tsv > raw_manifest.tsv

    paste "${batch_names}" '~{write_lines(batch_sample_lists)}' \
      | awk -F'\t' '{while((getline line < $2) > 0) {print $1 "\t" line}}' \
      | sort -u -k2,2 > 'sample_manifest.tsv'

    paste "${batch_names}" '~{write_lines(batch_bincov_matrix)}' \
      '~{write_lines(batch_bincov_matrix_index)}' > 'bincov_manifest.tsv'

    cut -f 1,2 'bincov_manifest.tsv' > 'bincov.tsv'
    cut -f 1,3 'bincov_manifest.tsv' > 'bincov_index.tsv'

duckdb <<'EOF'
COPY (
  SELECT json_group_object(batch, vcfs) AS raw_vcf
  FROM (
    SELECT batch, list(vcf) AS vcfs
    FROM read_csv('raw_manifest.tsv',
                  delim = '\t',
                  header = false,
                  names = ['batch', 'vcf'])
    GROUP BY batch
  )
) TO 'raw_manifest.json' (FORMAT JSON);
COPY (
  SELECT json_group_object(batch, bincov) AS bincov
  FROM read_csv('bincov.tsv',
                delim = '\t',
                header = false,
                names = ['batch', 'bincov'])
) TO 'bincov_manifest.json' (FORMAT JSON);
COPY (
  SELECT json_group_object(batch, bincov_index) AS bincov_index
  FROM read_csv('bincov_index.tsv',
                delim = '\t',
                header = false,
                names = ['batch', 'bincov_index'])
) TO 'bincov_index_manifest.json' (FORMAT JSON);
EOF
  >>>
}

# Retrieve a single contig from a VCF.
task SubsetVcfByContig {
  input {
    File vcf
    File vcf_index
    String contig
    String sv_base_mini_docker
    RuntimeAttr? runtime_attr_override
  }

  parameter_meta {
    vcf: "VCF to subset."
    vcf_index: "Index file of VCF to subset."
    contig: "Contig to extract."
    sv_base_mini_docker: "The corresponding Docker image from GATK-SV."
    runtime_attr_override: "Runtime attribute overrides."
  }

  output {
    File? subset_vcf = contig_vcf
    File? subset_vcf_index = "${contig_vcf}.tbi"
  }

  Float input_size = size([vcf, vcf_index], "GB")
  RuntimeAttr default_attr = object {
    mem_gb: 1,
    cpu_cores: 1,
    disk_gb: ceil(input_size * 1.1) + 16,
    boot_disk_gb: 8,
    preemptible_tries: 3,
    max_retries: 1,
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  runtime {
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GB"
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    docker: sv_base_mini_docker
  }

  String contig_vcf = "${contig}.vcf.gz"

  command <<<
    set -euxo pipefail

    bcftools view --regions '~{contig}' --output-type z --output '~{contig_vcf}' \
      '~{vcf}'

    read -r nrec < <(bcftools head --header 0 --records 1 '~{contig_vcf}' | wc -l)
    if (( nrec == 0 )); then
      rm '~{contig_vcf}'
      exit 0
    fi
    bcftools index --tbi '~{contig_vcf}'
  >>>
}

# Subset the pedigree, samples and batches so they are all synchronized.
# 1. Subset the pedigree to trios.
# 2. Subset the pedigree by families, if present.
# 3. Subset the pedigree to trios for which all members are in the VCF.
# 4. Subset the samples in the VCF to samples that are also in the pedigree.
# 5. Subset the batch manifest to batches with samples in the list from 4.
task SubsetSamples {
  input {
    File ped
    File? fams
    File vcf
    File sample_manifest
    String sv_base_mini_docker
    RuntimeAttr? runtime_attr_override
  }

  parameter_meta {
    ped: "Pedigree."
    fams: "Family IDs, one per line, to use to subset the pedigree."
    vcf: "VCF to use to synchronize."
    sample_manifest: "TSV with batches in the first column and samples in the second."
    sv_base_mini_docker: "The corresponding Docker image from GATK-SV."
    runtime_attr_override: "Runtime attribute overrides."
  }

  output {
    File ped_subset = "subset.ped"
    File offspring = "offspring.list"
    File batch_subset = "batch_subset.list"
  }

  Float input_size = size(select_all([ped, fams, vcf, sample_manifest]), "GB")
  RuntimeAttr default_attr = object {
    mem_gb: 2,
    cpu_cores: 1,
    disk_gb: ceil(input_size) + 16,
    boot_disk_gb: 8,
    preemptible_tries: 3,
    max_retries: 1,
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  runtime {
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GB"
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    docker: sv_base_mini_docker
  }

  command <<<
    set -euxo pipefail

    fam_ids='~{if defined(fams) then fams else ""}'
    awk -F'\t' '$2 && $3 && $4' '~{ped}' > trios.ped
    ped=trios.ped

    if [[ -n "${fam_ids:-}" && -s "${fam_ids}" ]]; then
      awk -F'\t' 'NR==FNR{a[$1]} NR>FNR && ($1 in a)' "${fam_ids}" "${ped}" > fam_subset.ped
      ped=fam_subset.ped
    fi

    bcftools query --list-samples '~{vcf}' | LC_ALL=C sort > vcf_samples.list
    awk -F'\t' 'NR==FNR{a[$1]} NR>FNR && ($2 in a) && ($3 in a) && ($4 in a)' \
      vcf_samples.list "${ped}" > subset.ped

    awk -F'\t' '{print $2; print $3; print $4}' subset.ped \
      | LC_ALL=C sort > ped_samples.list

    LC_ALL=C comm -12 ped_samples.list vcf_samples.list > sample_subset.list

    awk -F'\t' 'NR==FNR{a[$1]} NR>FNR && ($2 in a){print $1}' \
      sample_subset.list '~{sample_manifest}' \
      | LC_ALL=C sort -u > batch_subset.list

    read -r ped_n _ < <(wc -l subset.ped)
    if (( ped_n == 0 )); then
      printf 'pedigree is empty\n' >&2
      exit 1
    fi

    read -r sample_n _ < <(wc -l sample_subset.list)
    if (( sample_n == 0 )); then
      printf 'sample set is empty\n' >&2
      exit 1
    fi

    read -r batch_n _ < <(wc -l batch_subset.list)
    if (( batch_n == 0 )); then
      printf 'batch set is empty\n' >&2
      exit 1
    fi

    awk -F'\t' '{print $2}' subset.ped | sort -u > offspring.list
  >>>
}

# Find out which contig, among a set, a VCF contains. If the VCF contains more
# than one contig, the task will error. If the VCF does not contain any of the
# given contigs, it will not produce an output.
task MatchVcfToContig {
  input {
    File vcf
    File vcf_index
    Array[String] contigs
    String sv_base_mini_docker
    RuntimeAttr? runtime_attr_override

    # NOT AN INPUT! Only exists to create an optional type for use in outputs.
    String? null
  }

  parameter_meta {
    vcf: "VCF to match."
    vcf_index: "Index file of VCF to match."
    contigs: "Contigs to match."
    sv_base_mini_docker: "The corresponding Docker image from GATK-SV."
    runtime_attr_override: "Runtime attribute overrides."
  }

  # There doesn't seem to be way in the WDL language to specify an optional
  # `String` output without this hack of using a fake input for its optional
  # type. The `read_string` function will always return `String`, not
  # `String?` because it will error if its argument file does not exist. This
  # is problematic because this task relies on outputting optional values to
  # indicate that a VCF did not match a contig.
  output {
    File? matched_vcf = vcf_bn
    String? matched_contig = if read_string("matched_contig.list") == "" then null else read_string("matched_contig.list")
  }

  Float input_size = size([vcf, vcf_index], "GB")
  RuntimeAttr default_attr = object {
    mem_gb: 1,
    cpu_cores: 1,
    disk_gb: ceil(input_size) + 16,
    boot_disk_gb: 8,
    preemptible_tries: 3,
    max_retries: 1,
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  runtime {
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb])+ " GB"
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    docker: sv_base_mini_docker
  }

  String vcf_bn = basename(vcf)

  command <<<
    set -euxo pipefail

    bcftools index --stats '~{vcf}' | cut -f 1 > contigs_in_vcf.list
    read -r vcf_contigs_count _ < <(wc -l contigs_in_vcf.list)
    if (( vcf_contigs_count != 1 )); then
      printf 'VCF must contain exactly 1 contig. Found %d\n' "${vcf_contigs_count}" >&2
      exit 1
    fi

    awk 'NR==FNR{a=$1} NR>FNR && (a==$1){print a; exit 0}' \
      contigs_in_vcf.list \
      '~{write_lines(contigs)}' > matched_contig.list

    # If the contig in the VCF matched one of the input contigs, then
    # matched_contig.list should contain the matched contig and have a non-zero
    # filesize.
    mv '~{vcf}' '~{vcf_bn}'
    if [[ ! -s matched_contig.list ]]; then
      printf '\n' > matched_contig.list
      rm '~{vcf_bn}'
    fi
  >>>
}
