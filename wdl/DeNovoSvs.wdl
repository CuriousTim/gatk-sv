version 1.0

import "Structs.wdl"
import "DeNovoSvsGroupOffspringBcf.wdl"
import "DeNovoSvsSvConcordancePerContig.wdl"

###########################
# MAIN WORKFLOW DEFINITION
###########################

workflow DeNovoSvs {
  input {
    File pedigree
    # One family ID per line to call de novo in subset of families
    File? family_ids

    Float max_cohort_af = 0.02
    Float max_gnomad_af = 0.01
    Int large_cnv_size = 1000
    Int depth_only_size = 5000
    Array[File]? exclude_regions
    Float exclude_regions_ovp = 0.5
    File gd_regions
    Float gd_regions_ovp = 0.5

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
    RuntimeAttr? runtime_override_make_offspring_bcf
    RuntimeAttr? runtime_override_remove_uncalled_svtypes
    RuntimeAttr? runtime_override_filter_offspring_sites
    RuntimeAttr? runtime_override_match_bcf_to_contig
    RuntimeAttr? runtime_override_group_offspring_by_batch
    RuntimeAttr? runtime_override_group_bcf_by_family_batch
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
  if (size(SubsetSamples.offspring) > 0) {
    scatter (i in range(length(split_vcfs))) {
      call MakeOffspringBcf {
        input:
          vcf = split_vcfs[i],
          offspring = SubsetSamples.offspring,
          bcf_prefix = "offspring_${i}",
          sv_base_mini_docker = sv_base_mini_docker,
          runtime_attr_override = runtime_override_make_offspring_bcf
      }

      call RemoveUncalledSvtypes {
        input:
          bcf = MakeOffspringBcf.offspring_bcf,
          sv_base_mini_docker = sv_base_mini_docker,
          runtime_attr_override = runtime_override_remove_uncalled_svtypes
      }

      call FilterOffspringSites {
        input:
          bcf = RemoveUncalledSvtypes.filtered_bcf,
          max_cohort_af = max_cohort_af,
          max_gnomad_af = max_gnomad_af,
          large_cnv_size = large_cnv_size,
          depth_only_size = depth_only_size,
          exclude_regions = exclude_regions,
          exclude_regions_ovp = exclude_regions_ovp,
          gd_regions = gd_regions,
          gd_regions_ovp = gd_regions_ovp,
          sv_base_mini_docker = sv_base_mini_docker,
          runtime_attr_override = runtime_override_filter_offspring_sites
      }

      call MatchBcfToContig {
        input:
          bcf = FilterOffspringSites.filtered_bcf,
          contigs = contigs,
          sv_base_mini_docker = sv_base_mini_docker,
          runtime_attr_override = runtime_override_match_bcf_to_contig
      }
    }
  }

  Array[File] matched_bcfs = select_all(select_first([MatchBcfToContig.matched_bcf, []]))
  Array[String] kept_contigs = select_all(select_first([MatchBcfToContig.matched_contig, []]))

  scatter (bcf in matched_bcfs) {
    call DeNovoSvsGroupOffspringBcf.DeNovoSvsGroupOffspringBcf as group_offspring_bcf {
      input:
        bcf = bcf,
        offspring = SubsetSamples.offspring,
        batches = SubsetSamples.batch_subset,
        pedigree = SubsetSamples.ped_subset,
        sample_manifest = MakeManifests.sample_manifest,
        linux_docker = linux_docker,
        sv_base_mini_docker = sv_base_mini_docker,
        runtime_override_group_offspring_by_batch = runtime_override_group_offspring_by_batch,
        runtime_override_group_bcf_by_family_batch = runtime_override_group_bcf_by_family_batch
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
  Array[Array[File]] offspring_batch_grouped_bcfs = transpose(group_offspring_bcf.offspring_batch_grouped_bcfs)
  Array[Array[File]] father_batch_grouped_bcfs = transpose(group_offspring_bcf.father_batch_grouped_bcfs)
  Array[Array[File]] mother_batch_grouped_bcfs = transpose(group_offspring_bcf.mother_batch_grouped_bcfs)
  Array[Array[File]] offspring_batch_offspring_ids = transpose(group_offspring_bcf.offspring_batch_offspring_ids)
  Array[Array[File]] father_batch_father_ids = transpose(group_offspring_bcf.father_batch_father_ids)
  Array[Array[File]] mother_batch_mother_ids = transpose(group_offspring_bcf.mother_batch_mother_ids)

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

# Subset a VCF to offspring samples and convert to BCF.
task MakeOffspringBcf {
  input {
    File vcf
    File offspring
    String bcf_prefix
    String sv_base_mini_docker
    RuntimeAttr? runtime_attr_override
  }

  parameter_meta {
    vcf: "VCF file to subset."
    offspring: "Offspring samples, one per line."
    bcf_prefix: "Prefix to use for the output."
    sv_base_mini_docker: "The corresponding Docker image from GATK-SV."
    runtime_attr_override: "Runtime attribute overrides."
  }

  output {
    File offspring_bcf = offspring_name
  }

  RuntimeAttr default_attr = object {
    mem_gb: 4,
    cpu_cores: 1,
    disk_gb: ceil(size(vcf, "GB") * 2 + size(offspring, "GB")) + 16,
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

  String offspring_name = "${bcf_prefix}.bcf"

  command <<<
    set -euxo pipefail

    bcftools view --no-update --samples-file '~{offspring}' --output-type b \
      --output '~{offspring_name}' '~{vcf}'
  >>>
}

# Remove SV types from a BCF of offspring sites that are not handled in the de novo pipeline. All
# BND and CNV sites are removed. Then the remaining sites are split into two: those that are CPX or
# CTX and those that are not. The file with CPX and CTX events will be a VCF and the file without
# those SV types will be BCF.
task RemoveUncalledSvtypes {
  input {
    File bcf
    String sv_base_mini_docker
    RuntimeAttr? runtime_attr_override
  }

  parameter_meta {
    bcf: "BCF file to filter."
    sv_base_mini_docker: "The corresponding Docker image from GATK-SV."
    runtime_attr_override: "Runtime attribute overrides."
  }

  output {
    File filtered_bcf = filtered_bcf_name
    File cpx_vcf = cpx_vcf_name
  }

  RuntimeAttr default_attr = object {
    mem_gb: 4,
    cpu_cores: 1,
    disk_gb: ceil(size(bcf, "GB") * 3)  + 16,
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

  String filtered_bcf_name = "svtypes_filtered-" + basename(bcf)
  String cpx_vcf_name = "cpx_ctx-" + basename(bcf, ".bcf") + ".vcf.gz"

  command <<<
    set -euxo pipefail

    bcftools view --exclude 'INFO/SVTYPE = "BND" || INFO/SVTYPE = "CNV"' \
      --output-type b --output tmp.bcf '~{bcf}'

    bcftools view --exclude 'INFO/SVTYPE = "CPX" || INFO/SVTYPE = "CTX"' \
      --output-type b --output '~{filtered_bcf_name}' tmp.bcf
    bcftools view --include 'INFO/SVTYPE = "CPX" || INFO/SVTYPE = "CTX"' \
      --output-type z --output '~{cpx_vcf_name}' tmp.bcf
  >>>
}

# Filter sites in an offspring BCF for potential de novos
# 1. Remove all sites that:
#    a. have an cohort or gnomAD allele frequency greater than the input
#       thresholds
#    b. are overlapped by exclude regions by a minimum of
#       `exclude_regions_ovp` fraction of the SV
#    c. small CNVs that are SR-only and don't have BOTHSIDES_SUPPORT
#    d. are depth-only DUPs and are smaller than the depth-only size threshold
#    e. are not covered by genomic disorder regions by a minimum of
#       `gd_regions_ovp` fraction of the SV (any site meeting this criteria will be
#       kept, even if it would otherwise excluded by the previous criteria)
task FilterOffspringSites {
  input {
    File bcf
    Float max_cohort_af
    Float max_gnomad_af
    Int large_cnv_size
    Int depth_only_size
    Array[File]? exclude_regions
    Float exclude_regions_ovp
    File gd_regions
    Float gd_regions_ovp
    String sv_base_mini_docker
    RuntimeAttr? runtime_attr_override
  }

  parameter_meta {
    bcf: "BCF with offspring samples."
    max_cohort_af: "Maximum cohort allele frequency allowed."
    max_gnomad_af: "Maximum gnomAD allele frequency allowed. The value should be in the INFO field with the key 'gnomad_v4.1_sv_AF'."
    large_cnv_size: "Minimum size, in bases, of a large CNV."
    depth_only_size: "Minimum size, in bases, of a DUP that is allowed to a depth-only call."
    exclude_regions: "BED3 files of genomic regions to exclude. The files are concatenated before testing for coverage."
    exclude_regions_ovp: "Fraction of SV that must be covered by exclude regions to be dropped."
    gd_regions: "BED3 files of genomic disorder regions."
    gd_regions_ovp: "Fraction of SV that must be covered by genomic disorder regions to be dropped."
    sv_base_mini_docker: "The corresponding Docker image from GATK-SV."
    runtime_attr_override: "Runtime attribute overrides."
  }

  output {
    File filtered_bcf = filtered_bcf_name
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

  String filtered_bcf_name = "sites_filtered-${basename(bcf)}"

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
    bcftools query \
      --include 'AF > ~{max_cohort_af} || (gnomad_v4.1_sv_AF != "." && gnomad_v4.1_sv_AF > ~{max_gnomad_af})' \
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
      | awk -F'\t' '$8 >= ~{gd_regions_ovp} {print $4}' > gd_pass

    sort -u gd_pass > whitelist
    cat af_fail bothsides_fail depth_only_fail exclude_regions_fail | sort -u > blacklist
    comm -13 whitelist blacklist > blacklist_clean

    bcftools view --exclude 'ID = "@blacklist_clean"' --output-type u \
      --output '~{filtered_bcf_name}' '~{bcf}'
  >>>
}


# Find out which contig, among a set, a BCF contains. If the BCF contains more
# than one contig, the task will error. If the BCF does not contain any of the
# given contigs, it will not produce an output.
task MatchBcfToContig {
  input {
    File bcf
    Array[String] contigs
    String sv_base_mini_docker
    RuntimeAttr? runtime_attr_override

    # NOT AN INPUT! Only exists to create an optional type for use in outputs.
    String? null
  }

  parameter_meta {
    bcf: "BCF to match."
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
    File? matched_bcf = bcf_bn
    String? matched_contig = if read_string("matched_contig.list") == "" then null else read_string("matched_contig.list")
  }

  Float input_size = size(bcf, "GB")
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
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb])+ " GB"
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    docker: sv_base_mini_docker
  }

  String bcf_bn = basename(bcf)

  command <<<
    set -euxo pipefail

    bcftools index --stats '~{bcf}' | cut -f 1 > contigs_in_bcf.list
    read -r bcf_contigs_count _ < <(wc -l contigs_in_bcf.list)
    if (( bcf_contigs_count != 1 )); then
      printf 'BCF must contain exactly 1 contig. Found %d\n' "${vcf_contigs_count}" >&2
      exit 1
    fi

    awk 'NR==FNR{a=$1} NR>FNR && (a==$1){print a; exit 0}' \
      contigs_in_bcf.list \
      '~{write_lines(contigs)}' > matched_contig.list

    # If the contig in the BCF matched one of the input contigs, then
    # matched_contig.list should contain the matched contig and have a non-zero
    # filesize.
    mv '~{bcf}' '~{bcf_bn}'
    if [[ ! -s matched_contig.list ]]; then
      printf '\n' > matched_contig.list
      rm '~{bcf_bn}'
    fi
  >>>
}
