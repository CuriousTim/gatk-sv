version 1.0

###########################
# IMPORT TOOLS
###########################

import "Structs.wdl"

###########################
# MAIN WORKFLOW DEFINITION
###########################

workflow DeNovoSvs {
  input {
    File pedigree
    # One family ID per line to call de novo in subset of families
    File? family_ids

    #' VCF filter parameters
    Float max_cohort_af = 0.02
    Float max_gnomad_af = 0.01
    # Minimum fraction of SV overlapped by GD region
    Float gd_overlap = 0.5
    File gd_regions
    # BED3 files with hg38 exclude regions
    Array[File]? exclude_regions
    Float exclude_regions_ovp = 0.5

    Int large_cnv_size = 1000
    Int depth_only_size = 5000

    # Either a single VCF or an array of VCFs with each one containing a single
    # contig. In the case of a single VCF, it is expected that all the contigs
    # in the input contigs are present. In the case of multiple VCFs, all VCFs
    # must contain the exact same set of samples.
    Array[File]+ vcfs
    Array[File]+ vcf_indices

    Array[String]+ contigs = ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6",
      "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15",
      "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX"]

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
    String sv_pipeline_docker
    String denovo_docker

    RuntimeAttr? runtime_override_make_manifests
    RuntimeAttr? runtime_override_subset_vcf_by_contig
    RuntimeAttr? runtime_override_subset_samples
    RuntimeAttr? runtime_override_match_vcf_to_contig
    RuntimeAttr? runtime_override_filter_bcf_sites
    RuntimeAttr? runtime_override_split_bcf_by_samples
    RuntimeAttr? runtime_override_filter_proband_sites
    RuntimeAttr? runtime_override_filter_proband_genotypes
    RuntimeAttr? runtime_override_split_proband_bcf_by_batch
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
          denovo_docker = denovo_docker,
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
      batch_manifest = MakeManifests.sample_manifest,
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
  scatter (vcf in matched_vcfs) {
    call FilterVcfSites {
      input:
        vcf = vcf,
        sv_base_mini_docker = sv_base_mini_docker,
        runtime_attr_override = runtime_override_filter_bcf_sites
    }

    call SplitBcfBySamples {
      input:
        bcf = FilterVcfSites.filtered_bcf,
        proband_ids = SubsetSamples.probands,
        parent_ids = SubsetSamples.parents,
        sv_base_mini_docker = sv_base_mini_docker,
        runtime_attr_override = runtime_override_split_bcf_by_samples
    }

    call FilterProbandSites {
      input:
        bcf = SplitBcfBySamples.proband_bcf,
        max_cohort_af = max_cohort_af,
        max_gnomad_af = max_gnomad_af,
        large_cnv_size = large_cnv_size,
        depth_only_size = depth_only_size,
        exclude_regions = exclude_regions,
        exclude_regions_ovp = exclude_regions_ovp,
        gd_regions = gd_regions,
        gd_overlap = gd_overlap,
        sv_base_mini_docker = sv_base_mini_docker,
        runtime_attr_override = runtime_override_filter_proband_sites
    }

    call FilterProbandGenotypes {
      input:
        bcf = FilterProbandSites.filtered_bcf,
        sv_base_mini_docker = sv_base_mini_docker,
        runtime_attr_override = runtime_override_filter_proband_genotypes
    }

    call SplitProbandBcfByBatch {
      input:
        bcf = FilterProbandGenotypes.filtered_bcf,
        batch_n = length(batch_name_list),
        sample_manifest = MakeManifests.sample_manifest,
        sv_base_mini_docker = sv_base_mini_docker,
        runtime_attr_override = runtime_override_split_proband_bcf_by_batch
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
  Array[Array[File]] batched_proband_bcfs = transpose(SplitProbandBcfByBatch.split_bcfs)
  # scatter (batch_vcfs in batched_proband_vcfs) {
  #   if (size(batch_vcfs) > 0) {
  #     String batch_id = basename(batch_vcfs[0], ".vcf.gz")
  #     call CheckProbandRawEvidence {
  #       input:
  #         proband_vcfs = batch_vcfs,
  #         clustered_pesr_vcfs = MakeManifests.pesr_map[batch_id],
  #         clustered_depth_vcf = MakeManifests.depth_map[batch_id],
  #         sv_pipeline_docker = sv_pipeline_docker,
  #         runtime_attr_override = runtime_override_check_proband_raw_evidence
  #     }
  #   }
  # }

  output {
    Array[Array[File]] probands = batched_proband_bcfs
  }
}

###########################
# TASK DEFINITIONS
###########################

# Convert the input arrays into files that are needed by other tasks in the workflow.
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
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb])+ " GB"
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

    paste "${batch_names}" '~{write_lines(batch_sample_lists)}' \
      | awk -F'\t' '{while((getline line < $2) > 0) {print $1 "\t" line}}' > 'sample_manifest.tsv'

    paste "${batch_names}" '~{write_lines(batch_bincov_matrix)}' \
      '~{write_lines(batch_bincov_matrix_index)}' > 'bincov_manifest.tsv'

duckdb <<'EOF'
COPY (
  SELECT json_group_object(batch, vcfs) AS pesr
  FROM (
    SELECT batch, list(vcf) AS vcfs
    FROM read_csv('pesr_manifest.tsv',
                  delim = '\t',
                  header = false,
                  names = ['batch', 'vcf'])
    GROUP BY batch
  )
) TO 'pesr_manifest.json' (FORMAT JSON);
COPY (
  SELECT json_group_object(batch, vcf) AS depth
  FROM read_csv('depth_manifest.tsv',
                delim = '\t',
                header = false,
                names = ['batch', 'vcf'])
) TO 'depth_manifest.json' (FORMAT JSON);
COPY (
  SELECT json_group_object(batch, bincov) AS bincov
  FROM read_csv('bincov_manifest.tsv',
                delim = '\t',
                header = false,
                names = ['batch', 'bincov'])
) TO 'bincov_manifest.json' (FORMAT JSON);
EOF
  >>>

  # depth_manifest: BATCH_ID CLUSTERED_DEPTH_VCF_PATH
  # pesr_manifest: BATCH_ID CLUSTERED_*_VCF_PATH
  # bincov_manifest: BATCH_ID BINCOV_MATRIX_PATH BINCOV_MATRIX_INDEX_PATH
  # sample_manifest: BATCH_ID SAMPLE_ID
  output {
    File depth_manifest = "depth_manifest.tsv"
    File pesr_manifest = "pesr_manifest.tsv"
    File bincov_manifest = "bincov_manifest.tsv"
    File sample_manifest = "sample_manifest.tsv"
    Map[String, String] depth_map = read_json("depth_manifest.json")["depth"]
    Map[String, Array[String]] pesr_map = read_json("pesr_manifest.json")["pesr"]
    Map[String, String] bincov_map = read_json("bincov_manifest.json")["bincov"]
  }
}

# Retrieve a single contig from a VCF.
task SubsetVcfByContig {
  input {
    File vcf
    File vcf_index
    String contig
    String denovo_docker
    RuntimeAttr? runtime_attr_override
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
    docker: denovo_docker
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

  output {
    File? subset_vcf = contig_vcf
    File? subset_vcf_index = contig_vcf + ".tbi"
  }
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
    File batch_manifest
    String sv_base_mini_docker
    RuntimeAttr? runtime_attr_override
  }

  Float input_size = size(select_all([ped, fams, vcf, batch_manifest]), "GB")
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
      sample_subset.list '~{batch_manifest}' \
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

    awk -F'\t' '{print $2 > "probands.list"; print $3 > "parents.list"; print $4 > "parents.list"}' \
      subset.ped
    sort -u probands.list > probands.list.tmp
    mv probands.list.tmp probands.list
    sort -u parents.list > parents.list.tmp
    mv parents.list.tmp parents.list
  >>>

  output {
    File ped_subset = "subset.ped"
    File probands = "probands.list"
    File parents = "parents.list"
    File sample_subset = "sample_subset.list"
    File batch_subset = "batch_subset.list"
  }
}

# Remove BND and mCNV sites and create a BCF file without any CPX or CTX sites
# and a VCF with only the CPX and CTX sites
task FilterVcfSites {
  input {
    File vcf
    String sv_base_mini_docker
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
    mem_gb: 4,
    cpu_cores: 4,
    disk_gb: ceil(size(vcf, "GB") * 3)  + 16,
    boot_disk_gb: 8,
    preemptible_tries: 3,
    max_retries: 1,
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
  Int cpus = select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])

  runtime {
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GB"
    cpu: cpus
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    docker: sv_base_mini_docker
  }

  String output_name = basename(vcf, ".vcf.gz") + ".bcf"
  String cpx_output_name = "cpx_ctx-" + basename(vcf)

  command <<<
    set -euxo pipefail
    bcftools view --exclude 'SVTYPE = "BND" || SVTYPE = "CNV"' \
      --threads ~{cpus} --output-type b --output tmp.bcf '~{vcf}'

    bcftools view --exclude 'SVTYPE = "CPX" || SVTYPE = "CTX"' \
      --threads ~{cpus} --output-type b --output '~{output_name}' tmp.bcf
    bcftools view --include 'SVTYPE = "CPX" || SVTYPE = "CTX"' \
      --threads ~{cpus} --output-type z --output '~{cpx_output_name}' tmp.bcf
  >>>

  output {
    File filtered_bcf = output_name
    File cpx_vcf = cpx_output_name
  }
}

task SplitBcfBySamples {
  input {
    File bcf
    File proband_ids
    File parent_ids
    String sv_base_mini_docker
    RuntimeAttr? runtime_attr_override
  }

  Float bcf_size = size(bcf, "GB")
  Float other_size = size([proband_ids, parent_ids], "GB")
  RuntimeAttr default_attr = object {
    mem_gb: 4,
    cpu_cores: 4,
    disk_gb: ceil(bcf_size * 4 + other_size) + 16,
    boot_disk_gb: 8,
    preemptible_tries: 3,
    max_retries: 1,
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
  Int cpus = select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])

  runtime {
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GB"
    cpu: cpus
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    docker: sv_base_mini_docker
  }

  String proband_output = "proband-" + basename(bcf)
  String parent_output = "parent-" + basename(bcf)

  command <<<
    set -euxo pipefail

    bcftools view --output-type b --output '~{proband_output}' \
      --threads ~{cpus} --no-update --samples-file '~{proband_ids}' '~{bcf}'
    bcftools view --output-type b --output '~{parent_output}' \
      --threads ~{cpus} --no-update --samples-file '~{parent_ids}' '~{bcf}'
  >>>

  output {
    File proband_bcf = proband_output
    File parental_bcf = parent_output
  }
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

  # There doesn't seem to be way in the WDL language to specify an optional
  # `String` output without this hack of using a fake input for its optional
  # type. The `read_string` function will always return `String`, not
  # `String?` because it will error if its argument file does not exist. This
  # problematic because this task relies on outputting optional values to
  # indicate that a VCF did not match a contig.
  output {
    File? matched_vcf = vcf_bn
    String? matched_contig = if read_string("matched_contig.list") == "" then null else read_string("matched_contig.list")
  }
}

# Filter sites in a proband BCF for potential de novos
# 1. Remove all sites that:
#    a. have an cohort or gnomAD allele frequency greater than the input
#       thresholds
#    b. are overlapped by exclude regions by a minimum of
#       `exclude_regions_ovp` fraction of the SV
#    c. small CNVs that are SR-only and don't have BOTHSIDES_SUPPORT
#    d. are depth-only DUPs and are smaller than the depth-only size threshold
#    e. are not covered by genomic disorder regions by a minimum of
#       `gd_overlap` fraction of the SV (any site meeting this criteria will be
#       kept, even if it would otherwise excluded by the previous criteria)
task FilterProbandSites {
  input {
    File bcf
    Float max_cohort_af
    Float max_gnomad_af
    Int large_cnv_size
    Int depth_only_size
    Array[File]? exclude_regions
    Float exclude_regions_ovp
    File gd_regions
    Float gd_overlap
    String sv_base_mini_docker
    RuntimeAttr? runtime_attr_override
  }

  Float bcf_size = size(bcf, "GB")
  Float other_size = size(gd_regions, "GB") + (if defined(exclude_regions) then size(select_first([exclude_regions]), "GB") else 0)
  RuntimeAttr default_attr = object {
    mem_gb: 4,
    cpu_cores: 2,
    disk_gb: ceil(bcf_size * 4 + other_size) + 16,
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

  String output_bcf = basename(bcf)
  Array[File] bl = if defined(exclude_regions) then select_first([exclude_regions]) else []

  command <<<
    set -euxo pipefail

    cat2() {
      if [[ "$1" = *.gz ]]; then
        zcat "$1"
      else
        cat "$1"
      fi
    }

    bcftools view --drop-genotypes --output-type b --output sites_only.bcf '~{bcf}'
    # TODO should missing gnomAD AF mean filter out?
    bcftools query \
      --include 'AF > ~{max_cohort_af} || gnomad_v4.1_sv_AF = "." || gnomad_v4.1_sv_AF > ~{max_gnomad_af}' \
      --format '%ID\n' \
      sites_only.bcf > af_fail

    bcftools head sites_only.bcf | grep '^##' > headers.txt

    # Older GATK-SV VCFs have BOTHSIDES_SUPPORT in the FILTER field while
    # newer ones have it in the INFO field
    if grep -qF '##INFO=<ID=BOTHSIDES_SUPPORT,' headers.txt; then
      bothsides_filter='INFO/BOTHSIDES_SUPPORT = 1'
    elif grep -qF '##FILTER=<ID=BOTHSIDES_SUPPORT,' headers.txt; then
      bothsides_filter='FILTER ~ "BOTHSIDES_SUPPORT"'
    else
      printf 'BOTHSIDES_SUPPORT not found in BCF\n' >&2
      exit 1
    fi
    bcftools view \
      --include '(SVTYPE = "DEL" || SVTYPE = "DUP") && (EVIDENCE ~ "^RD,SR$" || EVIDENCE = "SR") && SVLEN < ~{large_cnv_size}' \
      --output-type u \
      sites_only.bcf \
      | bcftools query --exclude "${bothsides_filter}" --format '%ID\n' > bothsides_fail
    bcftools query \
      --include 'SVTYPE = "DUP" && ALGORITHMS = "depth" && SVLEN < ~{depth_only_size}' \
      --format '%ID\n' \
      sites_only.bcf > depth_only_fail

    bcftools query --format '%CHROM\t%POS0\t%END\t%ID\n' sites_only.bcf > sites.bed
    : > exclude_regions_fail
    er_paths='~{if defined(exclude_regions) then write_lines(select_first([exclude_regions])) else ""}'
    if [[ -n "${er_paths:-}" ]]; then
      while read -r f; do cat2 "${f}"; done < "${er_paths}" \
        | LC_ALL=C sort -k1,1 -k2,2n > er_merged.bed

        bedtools coverage -a sites.bed -b er_merged.bed -sorted \
          | awk -F'\t' '$8 >= ovp {print $4}' ovp=~{exclude_regions_ovp} >> exclude_regions_fail
    fi

    bedtools coverage -a sites.bed -b '~{gd_regions}' \
      | awk -F'\t' '$8 >= ~{gd_overlap} {print $4}' > gd_pass

    sort -u gd_pass > whitelist
    cat af_fail bothsides_fail depth_only_fail exclude_regions_fail | sort -u > blacklist
    comm -13 whitelist blacklist > blacklist_clean

    bcftools view --exclude 'ID = "@blacklist_clean"' --output-type u \
      --output '~{output_bcf}' '~{bcf}'
  >>>

  output {
    File filtered_bcf = output_bcf
  }
}

# Filter genotypes in a proband BCF for potential de novos
# 1. Set genotypes to missing where:
#    a. SV type is INS and algorithm is Manta or MELT and evidence is SR-only
#       and GQ = 0
#    b. evidence is Wham-only and GQ = 1
#    c. SV type is DEL and RD_CN is 2 or 3 and evidence is PE
task FilterProbandGenotypes {
  input {
    File bcf
    String sv_base_mini_docker 
    RuntimeAttr? runtime_attr_override
  }

  Float bcf_size = size(bcf, "GB")
  RuntimeAttr default_attr = object {
    mem_gb: 4,
    cpu_cores: 2,
    disk_gb: ceil(bcf_size * 3) + 16,
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

  String output_bcf = basename(bcf)

  command <<<
    set -euxo pipefail

    bcftools head '~{bcf}' | grep '^##' > headers.txt

    # Older GATK-SV VCFs have HIGH_SR_BACKGROUND in the FILTER field while
    # newer ones have it in the INFO field
    if grep -qF '##INFO=<ID=HIGH_SR_BACKGROUND,' headers.txt; then
      high_sr_filter='INFO/HIGH_SR_BACKGROUND = 1'
    elif grep -qF '##FILTER=<ID=HIGH_SR_BACKGROUND,' headers.txt; then
      high_sr_filter='FILTER ~ "HIGH_SR_BACKGROUND"'
    else
      printf 'HIGH_SR_BACKGROUND not found in BCF\n' >&2
      exit 1
    fi
    ins_filter='SVTYPE = "INS" & (ALGORITHMS = "manta" | ALGORITHMS = "melt") & (EVIDENCE ~ "^RS,SR$" | EVIDENCE = "SR") & GQ = 0'
    bcftools plugin setGT --output-type u '~{bcf}' -- \
      --target-gt q --new-gt . \
      --include "${ins_filter} & ${high_sr_filter}" \
      | bcftools plugin setGT --output-type u - -- \
          --target-gt q --new-gt . \
          --include 'EVIDENCE = "wham" & GQ = 1' \
      | bcftools plugin setGT --output-type b --output '~{output_bcf}' - -- \
          --target-gt q --new-gt . \
          --include 'SVTYPE = "DEL" & (RD_CN = 2 | RD_CN = 3) & EVIDENCE = "PE"'
  >>>

  output {
    File filtered_bcf = output_bcf
  }
}

# Split a BCF of proband samples into one BCF per batch
task SplitProbandBcfByBatch {
  input {
    File bcf
    Int batch_n
    File sample_manifest
    String sv_base_mini_docker
    RuntimeAttr? runtime_attr_override
  }

  Float input_size = size([bcf, sample_manifest], "GB")
  Int default_cpus = if batch_n < 8 then batch_n else 8
  RuntimeAttr default_attr = object {
    mem_gb: 1,
    cpu_cores: default_cpus,
    disk_gb: ceil(input_size * 3) + 16,
    boot_disk_gb: 8,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
  Int cpus = select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
  Float mem = select_first([runtime_attr.mem_gb, cpus * 2])

  runtime {
    memory: mem + " GB"
    cpu: cpus
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    docker: sv_base_mini_docker
  }

  command <<<
    set -euxo pipefail

    bcftools query --list-samples '~{bcf}' > bcf_samples
    mkdir batches bcfs
    # We need to unconditionally create a BCF for each batch for transpose to work
    cut -f 1 '~{sample_manifest}' | sort -u | xargs touch "bcfs/{}.bcf"
    awk -F'\t' 'NR==FNR{a[$2]=$1} NR>FNR && ($1 in a){print $1 > ("batches/" a[$1])}' \
      '~{sample_manifest}' bcf_samples

    find batches -type f -exec basename '{}' \; \
      | xargs -L 1 -P '~{cpus}' bcftools --output-type z --output "bcfs/{}.bcf" \
          --no-update --samples-files "batches/{}"
  >>>

  output {
    Array[File] split_bcfs = glob("bcfs/*.bcf")
  }
}
