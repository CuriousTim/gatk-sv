version 1.0

import "Structs.wdl"
import "DeNovoSvsGroupOffspringBcf.wdl"

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
    Float max_gd_af = 0.1

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
    RuntimeAttr? runtime_override_group_offspring_by_batch
    RuntimeAttr? runtime_override_make_offspring_bcf
    RuntimeAttr? runtime_override_remove_uncalled_svtypes
    RuntimeAttr? runtime_override_filter_offspring_sites
    RuntimeAttr? runtime_override_match_bcf_to_contig
    RuntimeAttr? runtime_override_merge_offspring_sites
    RuntimeAttr? runtime_override_group_bcf_by_family_batch
    RuntimeAttr? runtime_override_merge_clustered_batch_vcfs
    RuntimeAttr? runtime_override_sv_concordance
    RuntimeAttr? runtime_override_filter_genotypes
    RuntimeAttr? runtime_override_make_denovo_calls
    RuntimeAttr? runtime_override_merge_denovo_calls
  }

  output {
    File merged_denovos = MergeDeNovoCalls.merged_denovos
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

  Array[File] contig_vcfs = select_first([subset_vcfs, vcfs])
  call GroupOffspringByBatch {
    input:
      offspring = SubsetSamples.offspring,
      batches = SubsetSamples.batch_subset,
      pedigree = SubsetSamples.ped_subset,
      sample_manifest = MakeManifests.sample_manifest,
      linux_docker = linux_docker,
      runtime_attr_override = runtime_override_group_offspring_by_batch
  }

  scatter (i in range(length(contig_vcfs))) {
    call MakeOffspringBcf {
      input:
        vcf = contig_vcfs[i],
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
        max_gd_af = max_gd_af,
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

  Array[File] by_offspring_batch = GroupOffspringByBatch.by_offspring
  Array[File] by_father_batch = GroupOffspringByBatch.by_father
  Array[File] by_mother_batch = GroupOffspringByBatch.by_mother
  Array[File] father_ids = GroupOffspringByBatch.fathers
  Array[File] mother_ids = GroupOffspringByBatch.mothers
  Array[File] matched_bcfs = select_all(MatchBcfToContig.matched_bcf)
  Array[File] sites_only_matched_bcf = select_all(MatchBcfToContig.sites_only_matched_bcf)
  Array[String] kept_contigs = select_all(MatchBcfToContig.matched_contig)
  Array[String] kept_batches = read_lines(SubsetSamples.batch_subset)

  call MergeOffspringSites {
    input:
      bcfs = sites_only_matched_bcf,
      sv_base_mini_docker = sv_base_mini_docker,
      runtime_attr_override = runtime_override_merge_offspring_sites
  }

  scatter (bcf in matched_bcfs) {
    call DeNovoSvsGroupOffspringBcf.DeNovoSvsGroupOffspringBcf as group_offspring_bcf {
      input:
        bcf = bcf,
        by_offspring_batch = by_offspring_batch,
        by_father_batch = by_father_batch,
        by_mother_batch = by_mother_batch,
        batch_ids = kept_batches,
        sv_base_mini_docker = sv_base_mini_docker,
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

  scatter (i in range(length(by_offspring_batch))) {
    String current_batch = basename(by_offspring_batch[i])
    call MergeClusteredBatchVcfs {
      input:
        vcfs = MakeManifests.raw_vcf_map[current_batch],
        batch_id = current_batch,
        offspring_ids = by_offspring_batch[i],
        father_ids = father_ids[i],
        mother_ids = mother_ids[i],
        sv_base_mini_docker = sv_base_mini_docker,
        runtime_attr_override = runtime_override_merge_clustered_batch_vcfs
    }

    call SVConcordance {
      input:
        truth_vcf = MergeClusteredBatchVcfs.merged_vcf,
        truth_vcf_index = MergeClusteredBatchVcfs.merged_vcf_index,
        eval_vcf = MergeOffspringSites.merged_vcf,
        eval_vcf_index = MergeOffspringSites.merged_vcf_index,
        concordance_prefix = "${current_batch}-concordance",
        reference_dict = reference_dict,
        svconcordance_keep_all_docker = svconcordance_keep_all_docker,
        runtime_attr_override = runtime_override_sv_concordance
    }

    call FilterGenotypes {
      input:
        batch_id = current_batch,
        contigs = kept_contigs,
        by_offspring_batch_bcfs = offspring_batch_grouped_bcfs[i],
        by_father_batch_bcfs = father_batch_grouped_bcfs[i],
        by_mother_batch_bcfs = mother_batch_grouped_bcfs[i],
        offspring_genotypes = MergeClusteredBatchVcfs.offspring_genotypes,
        father_genotypes = MergeClusteredBatchVcfs.father_genotypes,
        mother_genotypes = MergeClusteredBatchVcfs.mother_genotypes,
        strict_concordance_vcf = SVConcordance.strict_concordance_vcf,
        lenient_concordance_vcf = SVConcordance.lenient_concordance_vcf,
        pedigree = SubsetSamples.ped_subset,
        denovo_docker = denovo_docker,
        runtime_attr_override = runtime_override_filter_genotypes
    }
  }

  Array[Array[File]] offspring_gt_filtered_tsvs = transpose(FilterGenotypes.self_filtered_tsvs)
  Array[Array[File]] father_gt_filtered_tsvs = transpose(FilterGenotypes.father_filtered_tsvs)
  Array[Array[File]] mother_gt_filtered_tsvs = transpose(FilterGenotypes.mother_filtered_tsvs)

  scatter (i in range(length(kept_contigs))) {
    call MakeDeNovoCalls {
      input:
        offspring_filtered_tsvs = offspring_gt_filtered_tsvs[i],
        father_filtered_tsvs = father_gt_filtered_tsvs[i],
        mother_filtered_tsvs = mother_gt_filtered_tsvs[i],
        contig = kept_contigs[i],
        linux_docker = linux_docker,
        runtime_attr_override = runtime_override_make_denovo_calls
    }
  }

  call MergeDeNovoCalls {
    input:
      denovos = MakeDeNovoCalls.denovos_tsv,
      linux_docker = linux_docker,
      runtime_attr_override = runtime_override_merge_denovo_calls
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
    memory: "${select_first([runtime_attr.mem_gb, default_attr.mem_gb])} GB"
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    disks: "local-disk ${select_first([runtime_attr.disk_gb, default_attr.disk_gb])} HDD"
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
    memory: "${select_first([runtime_attr.mem_gb, default_attr.mem_gb])} GB"
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    disks: "local-disk ${select_first([runtime_attr.disk_gb, default_attr.disk_gb])} HDD"
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
    memory: "${select_first([runtime_attr.mem_gb, default_attr.mem_gb])} GB"
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    disks: "local-disk ${select_first([runtime_attr.disk_gb, default_attr.disk_gb])} HDD"
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

# Group offspring IDs by different batches.
# The offspring IDs are grouped once by offspring batch, once by father batch and once by mother
# batch. The parental IDs are also grouped by their batches.
task GroupOffspringByBatch {
  input {
    File offspring
    File batches
    File pedigree
    File sample_manifest
    String linux_docker
    RuntimeAttr? runtime_attr_override
  }

  parameter_meta {
    offspring: "Offspring sample IDs, one per line."
    batches: "IDs of batches being processed, one per line."
    pedigree: "Pedigree."
    sample_manifest: "TSV of all samples (including parents) and batches being processed. First column is batch ID, second is sample ID."
    linux_docker: "A Linux Docker image."
    runtime_attr_override: "Runtime attribute overrides."
  }

  output {
    Array[File] by_offspring = glob("by_offspring/*")
    Array[File] by_father = glob("by_father/*")
    Array[File] fathers = glob("fathers/*")
    Array[File] by_mother = glob("by_mother/*")
    Array[File] mothers = glob("mothers/*")
  }

  Float input_size = size([offspring, batches, pedigree, sample_manifest], "GB")

  RuntimeAttr default_attr = object {
    mem_gb: 2,
    cpu_cores: 1,
    disk_gb: ceil(input_size * 2) + 16,
    boot_disk_gb: 8,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  runtime {
    memory: "${select_first([runtime_attr.mem_gb, default_attr.mem_gb])} GB"
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    disks: "local-disk ${select_first([runtime_attr.disk_gb, default_attr.disk_gb])} HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    docker: linux_docker
  }

  command <<<
    set -euxo pipefail

    mkdir by_offspring by_father fathers by_mother mothers
    # We need to unconditionally create a file for each batch for transpose to work
    while read -r b; do
      touch "by_offspring/${b}" "by_father/${b}" "fathers/${b}" "by_mother/${b}" "mothers/${b}"
    done < <(sort -u '~{batches}')

    awk -F'\t' 'FILENAME == ARGV[1] {a[$2]=$1}
                FILENAME == ARGV[2] {b[$2]=$3;c[$2]=$4}
                FILENAME == ARGV[3] {
                  print $1 > ("by_offspring/" a[$1])
                  print $1 > ("by_father/" a[b[$1]])
                  print b[$1] | ("sort -u > fathers/" a[b[$1]])
                  print $1 > ("by_mother/" a[c[$1]])
                  print c[$1] | ("sort -u > mothers/" a[c[$1]])
                }' '~{sample_manifest}' '~{pedigree}' '~{offspring}'
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
    memory: "${select_first([runtime_attr.mem_gb, default_attr.mem_gb])} GB"
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    disks: "local-disk ${select_first([runtime_attr.disk_gb, default_attr.disk_gb])} HDD"
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
    memory: "${select_first([runtime_attr.mem_gb, default_attr.mem_gb])} GB"
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    disks: "local-disk ${select_first([runtime_attr.disk_gb, default_attr.disk_gb])} HDD"
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
    Float max_gd_af
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
    gd_regions_ovp: "Fraction of SV that must be covered by genomic disorder regions to bypass site filters."
    max_gd_af: "Maximum allele frequency of a site to be considered for genomic disorder regions overlap."
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
    memory: "${select_first([runtime_attr.mem_gb, default_attr.mem_gb])} GB"
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    disks: "local-disk ${select_first([runtime_attr.disk_gb, default_attr.disk_gb])} HDD"
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

    bcftools query --include '(INFO/AF = "." || INFO/AF <= ~{max_gd_af}) && (INFO/SVTYPE = "DEL" || INFO/SVTYPE = "DUP")' \
      --format '%CHROM\t%POS0\t%END\t%ID\n' sites_only.bcf > gd_candidates.bed
    bedtools coverage -a gd_candidates.bed -b '~{gd_regions}' \
      | awk -F'\t' '$8 >= ~{gd_regions_ovp} {print $4}' > gd_pass

    sort -u gd_pass > whitelist
    cat af_fail bothsides_fail depth_only_fail exclude_regions_fail | sort -u > blacklist
    comm -13 whitelist blacklist > blacklist_clean

    bcftools view --exclude 'ID = @blacklist_clean' --output-type u \
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

    String? null
  }

  parameter_meta {
    bcf: "BCF to match."
    contigs: "Contigs to match."
    sv_base_mini_docker: "The corresponding Docker image from GATK-SV."
    runtime_attr_override: "Runtime attribute overrides."
    null: "Not an input. Only exists to create an optional type for use in outputs."
  }

  # There doesn't seem to be way in the WDL language to specify an optional
  # `String` output without this hack of using a fake input for its optional
  # type. The `read_string` function will always return `String`, not
  # `String?` because it will error if its argument file does not exist. This
  # is problematic because this task relies on outputting optional values to
  # indicate that a VCF did not match a contig.
  output {
    File? matched_bcf = bcf_name
    File? sites_only_matched_bcf = sites_only_bcf_name
    String? matched_contig = if read_string("matched_contig.list") == "" then null else read_string("matched_contig.list")
  }

  RuntimeAttr default_attr = object {
    mem_gb: 8,
    cpu_cores: 1,
    disk_gb: ceil(size(bcf, "GB") * 4) + 16,
    boot_disk_gb: 8,
    preemptible_tries: 3,
    max_retries: 1,
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  Float mem = select_first([runtime_attr.mem_gb, default_attr.mem_gb])
  Float max_sort_mem = mem * 0.8
  runtime {
    memory: "${mem} GB"
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    disks: "local-disk ${select_first([runtime_attr.disk_gb, default_attr.disk_gb])} HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    docker: sv_base_mini_docker
  }

  String bcf_name = basename(bcf)
  String sites_only_bcf_name = "sites_only-${bcf_name}"

  command <<<
    set -euxo pipefail

    bcftools view --drop-genotypes --output-type u '~{bcf}' \
      | bcftools sort --max-mem '~{max_sort_mem}G' --output '~{sites_only_bcf_name}' \
          --output-type b
    bcftools index '~{sites_only_bcf_name}'
    bcftools index --stats '~{sites_only_bcf_name}' | cut -f 1 > contigs_in_bcf.list
    read -r bcf_contigs_count _ < <(wc -l contigs_in_bcf.list)
    if (( bcf_contigs_count != 1 )); then
      printf 'BCF must contain exactly 1 contig. Found %d\n' "${bcf_contigs_count}" >&2
      exit 1
    fi

    awk 'NR==FNR{a=$1} NR>FNR && (a==$1){print a; exit 0}' \
      contigs_in_bcf.list \
      '~{write_lines(contigs)}' > matched_contig.list

    # If the contig in the BCF matched one of the input contigs, then
    # matched_contig.list should contain the matched contig and have a non-zero
    # filesize.
    mv '~{bcf}' '~{bcf_name}'
    if [[ ! -s matched_contig.list ]]; then
      printf '\n' > matched_contig.list
      rm '~{bcf_name}'
      rm '~{sites_only_bcf_name}'
    fi
  >>>
}

# Merge sites-only, per-contig BCFs into a single VCF.
task MergeOffspringSites {
  input {
    Array[File] bcfs
    String sv_base_mini_docker
    RuntimeAttr? runtime_attr_override
  }

  parameter_meta {
    bcfs: "BCFs to merge. Each file should be a sites-only BCF with a single contig."
    sv_base_mini_docker: "The corresponding Docker image from GATK-SV."
    runtime_attr_override: "Runtime attribute overrides."
  }

  output {
    File merged_vcf = vcf_name
    File merged_vcf_index = "${vcf_name}.tbi"
  }

  RuntimeAttr default_attr = object {
    mem_gb: 4,
    cpu_cores: 1,
    disk_gb: ceil(size(bcfs, "GB") * 4) + 16,
    boot_disk_gb: 8,
    preemptible_tries: 3,
    max_retries: 1,
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  runtime {
    memory: "${select_first([runtime_attr.mem_gb, default_attr.mem_gb])} GB"
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    disks: "local-disk ${select_first([runtime_attr.disk_gb, default_attr.disk_gb])} HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    docker: sv_base_mini_docker
  }

  String vcf_name = "merged.vcf.gz"

  command <<<
    set -euxo pipefail

    bcftools concat --file-list '~{write_lines(bcfs)}' --output '~{vcf_name}' \
      --output-type z
    bcftools index --tbi '~{vcf_name}'
  >>>
}

# Merge raw evidence VCFs.
task MergeClusteredBatchVcfs {
  input {
    Array[File]+ vcfs
    String batch_id
    File offspring_ids
    File father_ids
    File mother_ids
    String sv_base_mini_docker
    RuntimeAttr? runtime_attr_override
  }

  parameter_meta {
    vcfs: "VCFs from a single batch to merge."
    batch_id: "Batch ID of this batch."
    offspring_ids: "IDs of offspring in this batch."
    father_ids: "IDs of fathers in this batch."
    mother_ids: "IDs of mothers in this batch."
    sv_base_mini_docker: "The corresponding Docker image from GATK-SV."
    runtime_attr_override: "Runtime attribute overrides."
  }

  output {
    File merged_vcf = "${batch_id}-sites_only-merged.vcf.gz"
    File merged_vcf_index = "${batch_id}-sites_only-merged.vcf.gz.tbi"
    File offspring_genotypes = "${batch_id}-offspring_genotypes.tsv.gz"
    File father_genotypes = "${batch_id}-father_genotypes.tsv.gz"
    File mother_genotypes = "${batch_id}-mother_genotypes.tsv.gz"
  }

  RuntimeAttr default_attr = object {
    mem_gb: 4,
    cpu_cores: 1,
    disk_gb: ceil(size(vcfs, "GB") * 4) + 16,
    boot_disk_gb: 8,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  Float mem = select_first([runtime_attr.mem_gb, default_attr.mem_gb])
  Float max_sort_mem = mem * 0.8

  runtime {
    memory: "${mem} GB"
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    disks: "local-disk ${select_first([runtime_attr.disk_gb, default_attr.disk_gb])} HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    docker: sv_base_mini_docker
  }

  command <<<
    set -euxo pipefail

    : > groups.tsv
    if [[ -s '~{offspring_ids}' ]]; then
      awk '{print $1"\t-\toffspring.bcf"}' '~{offspring_ids}' >> groups.tsv
    fi
    if [[ -s '~{father_ids}' ]]; then
      awk '{print $1"\t-\tfather.bcf"}' '~{father_ids}' >> groups.tsv
    fi
    if [[ -s '~{mother_ids}' ]]; then
      awk '{print $1"\t-\tmother.bcf"}' '~{mother_ids}' >> groups.tsv
    fi

    vcfs='~{write_lines(vcfs)}'
    declare -i i=0
    while read -r vcf; do
      dest_dir="algo_${i}"
      mkdir "${dest_dir}"
      bcftools +split --groups-file groups.tsv --output "${dest_dir}" --output-type b \
        --exclude 'INFO/SVTYPE = "BND" || INFO/SVTYPE = "CPX" || INFO/SVTYPE = "CTX" || INFO/SVTYPE = "CNV"' \
        "${vcf}"
      bcftools view --drop-genotypes --output-type b \
        --exclude 'INFO/SVTYPE = "BND" || INFO/SVTYPE = "CPX" || INFO/SVTYPE = "CTX" || INFO/SVTYPE = "CNV"' \
        --output "${dest_dir}/sites_only.bcf" "${vcf}"
      bcftools index "${dest_dir}/sites_only.bcf"
      i=$((i + 1))
    done < "${vcfs}"

    find . -type f -name 'offspring.bcf' \
      | xargs -L 1 bcftools query --include 'GT="alt"' --format '%ID[\t%SAMPLE]\n' \
      | gzip -c > '~{batch_id}-offspring_genotypes.tsv.gz'

    find . -type f -name 'father.bcf' \
      | xargs -L 1 bcftools query --include 'GT="alt"' --format '%ID[\t%SAMPLE]\n' \
      | gzip -c > '~{batch_id}-father_genotypes.tsv.gz'

    find . -type f -name 'mother.bcf' \
      | xargs -L 1 bcftools query --include 'GT="alt"' --format '%ID[\t%SAMPLE]\n' \
      | gzip -c > '~{batch_id}-mother_genotypes.tsv.gz'

    bcftools concat --file-list <(find . -type f -name 'sites_only.bcf') --allow-overlaps \
      --output-type u \
      | bcftools sort --max-mem '~{max_sort_mem}G' --output '~{batch_id}-sites_only-merged.vcf.gz' \
        --output-type z
    bcftools index --tbi '~{batch_id}-sites_only-merged.vcf.gz'
  >>>
}

task SVConcordance {
  input {
    File truth_vcf
    File truth_vcf_index
    File eval_vcf
    File eval_vcf_index
    String concordance_prefix
    File reference_dict
    String svconcordance_keep_all_docker
    RuntimeAttr? runtime_attr_override
  }

  parameter_meta {
    truth_vcf: "VCF against which to match variants."
    truth_vcf_index: "Index file of the truth VCF."
    eval_vcf: "VCF to annotate with variants matched in the truth VCF."
    eval_vcf_index: "Index file of the evaluation VCF."
    concordance_prefix: "Prefix of the output VCF."
    reference_dict: "Sequence dictionary in the form of a '.dict' file."
    svconcordance_keep_all_docker: "Docker with a build of GATK that supports the `--keep-all` option of SVConcordance."
    runtime_attr_override: "Runtime attribute overrides."
  }

  output {
    File strict_concordance_vcf = "${strict_concordance_name}"
    File lenient_concordance_vcf = "${lenient_concordance_name}"
  }

  RuntimeAttr default_attr = object {
    mem_gb: 4,
    cpu_cores: 1,
    disk_gb: ceil(size([truth_vcf, eval_vcf], "GB") * 2) + 16,
    boot_disk_gb: 8,
    preemptible_tries: 3,
    max_retries: 1,
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  Float mem = select_first([runtime_attr.mem_gb, default_attr.mem_gb])
  Int jvm_mem = floor(mem * 800)
  runtime {
    memory: "${mem} GB"
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    disks: "local-disk ${select_first([runtime_attr.disk_gb, default_attr.disk_gb])} HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    docker: svconcordance_keep_all_docker
  }

  String strict_concordance_name = "${concordance_prefix}-strict.vcf.gz"
  String lenient_concordance_name = "${concordance_prefix}-lenient.vcf.gz"

  command <<<
    set -euxo pipefail

    printf 'NAME\tSVTYPE\tMIN_SIZE\tMAX_SIZE\tTRACKS\n' > stratify.tsv
    printf 'DEL_small\tDEL\t-1\t5000\tNULL\n' >> stratify.tsv
    printf 'DUP_small\tDUP\t-1\t5000\tNULL\n' >> stratify.tsv
    printf 'DEL_large\tDEL\t5000\t-1\tNULL\n' >> stratify.tsv
    printf 'DUP_large\tDUP\t5000\t-1\tNULL\n' >> stratify.tsv
    printf 'INV_small\tINV\t-1\t5000\tNULL\n' >> stratify.tsv
    printf 'INV_large\tINV\t5000\t-1\tNULL\n' >> stratify.tsv
    printf 'INS\tINS\t-1\t-1\tNULL\n' >> stratify.tsv

    printf 'NAME\tRECIPROCAL_OVERLAP\tSIZE_SIMILARITY\tBREAKEND_WINDOW\tSAMPLE_OVERLAP\n' > cluster_strict.tsv
    printf 'DEL_small\t0.1\t0.5\t300\t0\n' >> cluster_strict.tsv
    printf 'DUP_small\t0.1\t0.5\t300\t0\n' >> cluster_strict.tsv
    printf 'DEL_large\t0.8\t0\t100000000\t0\n' >> cluster_strict.tsv
    printf 'DUP_large\t0.8\t0\t100000000\t0\n' >> cluster_strict.tsv
    printf 'INV_small\t0.1\t0.5\t300\t0\n' >> cluster_strict.tsv
    printf 'INV_large\t0.8\t0\t100000000\t0\n' >> cluster_strict.tsv
    printf 'INS\t0\t0.5\t300\t0\n' >> cluster_strict.tsv

    printf 'NAME\tRECIPROCAL_OVERLAP\tSIZE_SIMILARITY\tBREAKEND_WINDOW\tSAMPLE_OVERLAP\n' > cluster_lenient.tsv
    printf 'DEL_small\t0.1\t0\t500\t0\n' >> cluster_lenient.tsv
    printf 'DUP_small\t0.1\t0\t500\t0\n' >> cluster_lenient.tsv
    printf 'DEL_large\t0.5\t0\t100000000\t0\n' >> cluster_lenient.tsv
    printf 'DUP_large\t0.5\t0\t100000000\t0\n' >> cluster_lenient.tsv
    printf 'INV_small\t0.1\t0\t500\t0\n' >> cluster_lenient.tsv
    printf 'INV_large\t0.5\t0\t100000000\t0\n' >> cluster_lenient.tsv
    printf 'INS\t0\t0\t500\t0\n' >> cluster_lenient.tsv

    gatk --java-options '-Xmx~{jvm_mem}M' SVConcordance \
      --keep-all \
      --sequence-dictionary '~{reference_dict}' \
      --eval '~{eval_vcf}' \
      --truth '~{truth_vcf}'\
      --output '~{strict_concordance_name}' \
      --stratify-config stratify.tsv \
      --clustering-config cluster_strict.tsv

    gatk --java-options '-Xmx~{jvm_mem}M' SVConcordance \
      --keep-all \
      --sequence-dictionary '~{reference_dict}' \
      --eval '~{eval_vcf}' \
      --truth '~{truth_vcf}'\
      --output '~{lenient_concordance_name}' \
      --stratify-config stratify.tsv \
      --clustering-config cluster_lenient.tsv
  >>>
}

# Filter genotypes per batch
task FilterGenotypes {
  input {
    String batch_id
    Array[String]+ contigs
    Array[File] by_offspring_batch_bcfs
    Array[File] by_father_batch_bcfs
    Array[File] by_mother_batch_bcfs
    File offspring_genotypes
    File father_genotypes
    File mother_genotypes
    File strict_concordance_vcf
    File lenient_concordance_vcf
    File pedigree
    String denovo_docker
    RuntimeAttr? runtime_attr_override
  }

  parameter_meta {
    batch_id: "Batch ID of this batch."
    contigs: "Contigs to keep."
    by_offspring_batch_bcfs: "Offsprings BCFs grouped by offspring batch, split by contig."
    by_father_batch_bcfs: "Offspring BCFs grouped by father batch, split by contig."
    by_mother_batch_bcfs: "Offspring BCFs grouped by mother batch, split by contig."
    offspring_genotypes: "Offspring genotypes from the MergeClusteredBatchVcfs."
    father_genotypes: "Father genotypes from the MergeClusteredBatchVcfs."
    mother_genotypes: "Mother genotypes from the MergeClusteredBatchVcfs."
    strict_concordance_vcf: "SVConcordance VCF between offspring sites and ClusterBatch sites with strict parameters."
    lenient_concordance_vcf: "SVConcordance VCF between offspring sites and ClusterBatch sites with lenient parameters."
    pedigree: "Cohort pedigree."
    denovo_docker: "The corresponding Docker image from GATK-SV."
    runtime_attr_override: "Runtime attribute overrides."
  }

  output {
    Array[File] self_filtered_tsvs = glob("offspring/*.tsv.gz")
    Array[File] father_filtered_tsvs = glob("father/*.tsv.gz")
    Array[File] mother_filtered_tsvs = glob("mother/*.tsv.gz")
  }

  Float disk_size = size(by_offspring_batch_bcfs, "GB")
    + size(by_father_batch_bcfs, "GB")
    + size(by_mother_batch_bcfs, "GB")
    + size(offspring_genotypes, "GB")
    + size(father_genotypes, "GB")
    + size(mother_genotypes, "GB")
    + size(strict_concordance_vcf, "GB")
    + size(lenient_concordance_vcf, "GB")

  RuntimeAttr default_attr = object {
    mem_gb: 4,
    cpu_cores: 1,
    disk_gb: ceil(disk_size * 3) + 16,
    boot_disk_gb: 8,
    preemptible_tries: 3,
    max_retries: 1,
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  runtime {
    memory: "${select_first([runtime_attr.mem_gb, default_attr.mem_gb])} GB"
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    disks: "local-disk ${select_first([runtime_attr.disk_gb, default_attr.disk_gb])} HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    docker: denovo_docker
  }

  command <<<
    set -euxo pipefail

    write_gts() {
      # a file needs to be written for each contig unconditionally for transpose
      # to work
      # the numeric prefix is for glob to keep the order of the contigs
      awk 'NR == FNR {
        a[$1] = sprintf("%s/%03d-%s.tsv.gz", dir, i++, bid)
        system("touch " a[$1])
      }
      NR > FNR && ($1 in a) {
        cmd = "gzip -c > " a[$1]
        print $0 | cmd
      }' dir="$1" bid='~{batch_id}' '~{write_lines(contigs)}' -
    }

    bcftools concat --file-list '~{write_lines(by_offspring_batch_bcfs)}' \
     --output by_offspring.bcf --output-type b
    bcftools concat --file-list '~{write_lines(by_father_batch_bcfs)}' \
      --output by_father.bcf --output-type b
    bcftools concat --file-list '~{write_lines(by_mother_batch_bcfs)}' \
      --output by_mother.bcf --output-type b

    filtergt by_offspring.bcf '~{strict_concordance_vcf}' \
      '~{offspring_genotypes}' \
      'self_filtered.bcf'
    cut -f2,3 '~{pedigree}' > fathers.tsv
    filtergt by_father.bcf '~{lenient_concordance_vcf}' \
      '~{father_genotypes}' \
      'father_filtered.bcf' \
      fathers.tsv
    cut -f2,4 '~{pedigree}' > mothers.tsv
    filtergt by_mother.bcf '~{lenient_concordance_vcf}' \
      '~{mother_genotypes}' \
      'mother_filtered.bcf' \
      mothers.tsv

    mkdir offspring father mother

    bcftools query --include 'GT == "alt"' \
      --format '[%CHROM\t%POS\t%INFO/END\t%INFO/SVLEN\t%ID\t%INFO/SVTYPE\t%SAMPLE\n]' \
      self_filtered.bcf \
      | write_gts offspring
    bcftools query --include 'GT == "alt"' \
      --format '[%CHROM\t%POS\t%INFO/END\t%INFO/SVLEN\t%ID\t%INFO/SVTYPE\t%SAMPLE\n]' \
      father_filtered.bcf\
      | write_gts father
    bcftools query --include 'GT == "alt"' \
      --format '[%CHROM\t%POS\t%INFO/END\t%INFO/SVLEN\t%ID\t%INFO/SVTYPE\t%SAMPLE\n]' \
      mother_filtered.bcf \
      | write_gts mother
  >>>
}

task MakeDeNovoCalls {
  input {
    Array[File] offspring_filtered_tsvs
    Array[File] father_filtered_tsvs
    Array[File] mother_filtered_tsvs
    String contig
    String linux_docker
    RuntimeAttr? runtime_attr_override
  }

  parameter_meta {
    offspring_filtered_tsvs: "Offspring genotypes filtered against offspring."
    father_filtered_tsvs: "Offspring genotypes filtered against father."
    mother_filtered_tsvs: "Offspring genotypes filtered against mother."
    contig: "Contig being processed."
    linux_docker: "The corresponding Docker image from GATK-SV."
    runtime_attr_override: "Runtime attribute overrides."
  }

  output {
    File denovos_tsv = "${contig}-denovos.tsv.gz"
  }

  Float disk_size = size(offspring_filtered_tsvs, "GB")
    + size(father_filtered_tsvs, "GB")
    + size(mother_filtered_tsvs, "GB")

  RuntimeAttr default_attr = object {
    mem_gb: 4,
    cpu_cores: 1,
    disk_gb: ceil(disk_size * 3) + 16,
    boot_disk_gb: 8,
    preemptible_tries: 3,
    max_retries: 1,
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  runtime {
    memory: "${select_first([runtime_attr.mem_gb, default_attr.mem_gb])} GB"
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    disks: "local-disk ${select_first([runtime_attr.disk_gb, default_attr.disk_gb])} HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    docker: linux_docker
  }

  command <<<
    set -euxo pipefail

    mkdir tmp
    mkdir tmp/offspring tmp/father tmp/mother
    cat '~{write_lines(offspring_filtered_tsvs)}' \
      | xargs mv -t tmp/offspring
    cat '~{write_lines(father_filtered_tsvs)}' \
      | xargs mv -t tmp/father
    cat '~{write_lines(mother_filtered_tsvs)}' \
      | xargs mv -t tmp/mother
    printf 'chr\tstart\tend\tsvlen\tname\tsvtype\tsample\n' \
      | gzip -c > '~{contig}-denovos.tsv.gz'
    find tmp -type f '!' -empty -exec gzip -cd '{}' \; \
      | LC_ALL=C sort \
      | uniq -c \
      | awk '$1==3{sub(/^ *[0-9]+ /, ""); print}' \
      | gzip -c >> '~{contig}-denovos.tsv.gz'
  >>>
}

task MergeDeNovoCalls {
  input {
    Array[File]+ denovos
    String linux_docker
    RuntimeAttr? runtime_attr_override
  }

  parameter_meta {
    denovos: "TSVs with de novo calls."
    linux_docker: "The corresponding Docker image from GATK-SV."
    runtime_attr_override: "Runtime attribute overrides."
  }

  output {
    File merged_denovos  = "denovo_svs.tsv.gz"
  }

  Float disk_size = size(denovos, "GB")

  RuntimeAttr default_attr = object {
    mem_gb: 4,
    cpu_cores: 1,
    disk_gb: ceil(disk_size * 2) + 16,
    boot_disk_gb: 8,
    preemptible_tries: 3,
    max_retries: 1,
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  runtime {
    memory: "${select_first([runtime_attr.mem_gb, default_attr.mem_gb])} GB"
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    disks: "local-disk ${select_first([runtime_attr.disk_gb, default_attr.disk_gb])} HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    docker: linux_docker
  }

  command <<<
    set -euxo pipefail

    manifest='~{write_lines(denovos)}'
    cp "$(head -n 1 "${manifest}")" 'denovo_svs.tsv.gz'
    awk 'NR>1' "${manifest}" \
      | while read -r f; do gzip -cd "${f}" | awk 'NR>1'; done \
      | gzip -c >> 'denovo_svs.tsv.gz'
  >>>
}
