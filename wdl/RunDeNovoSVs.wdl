version 1.0

###########################
# IMPORT TOOLS
###########################

import "Structs.wdl"
import "TasksMakeCohortVcf.wdl" as miniTasks
import "DeNovoSVsScatter.wdl" as runDeNovo

###########################
# MAIN WORKFLOW DEFINITION
###########################

workflow DeNovoSV {
  input {
    # Core input files
    File pedigree
    # One family ID per line to call de novo in subset of families
    File? family_ids

    # VCF filter parameters
    Float max_cohort_af = 0.02
    Float max_gnomad_af = 0.01
    # Minimum fraction of SV overlapped by GD region
    Float gd_overlap = 0.05
    File gd_regions

    # Either a single VCF or an array of VCFs with each one containing a single
    # contig. In the case of a single VCF, it is expected that all the contigs
    # in the input contigs are present. In the case of multiple VCFs, all VCFs
    # must contain the exact same set of samples.
    Array[File] vcfs
    Array[File] vcf_indicies
    # BED3 files with hg38 blacklist regions
    Array[File]? blacklists
    Float blacklist_overlap = 0.5

    # Running parameters
    String output_prefix
    Int records_per_shard
    Array[String] contigs = ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6",
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

    # Parameters for denovo_svs.py
    Int? small_cnv_size
    Int? intermediate_cnv_size
    Int? depth_only_size
    Int? exclude_parent_cnv_size
    Float? gnomad_af
    Float? parents_af
    Float? cohort_af
    Float? large_raw_overlap
    Float? small_raw_overlap
    Float? parents_overlap
    Int? nearby_insertion
    Int? coverage_cutoff
    Int? gq_min
    String? gnomad_col
    String? alt_gnomad_col

    # Parameters for denovo_outliers.py with default values
    Int denovo_outlier_factor = 3

    # Dockers
    String linux_docker
    String sv_base_mini_docker
    String sv_pipeline_docker
    String variant_interpretation_docker


    RuntimeAttr? runtime_override_make_manifests
    RuntimeAttr? runtime_override_subset_vcf_by_contig
    RuntimeAttr? runtime_override_filter_vcf
    RuntimeAttr? runtime_override_match_vcf_to_contig
    RuntimeAttr? runtime_override_subset_samples
    RuntimeAttr? runtime_override_subset_vcf
    RuntimeAttr? runtime_override_subset_manifests
    RuntimeAttr? runtime_override_clustered_vcf_to_bed
    RuntimeAttr? runtime_override_merge_clustered_bed
    RuntimeAttr? runtime_override_reformat_bed
    RuntimeAttr? runtime_override_shard_vcf
    RuntimeAttr? runtime_override_denovo
    RuntimeAttr? runtime_override_denovo_merge_bed
    RuntimeAttr? runtime_override_vcf_to_bed
    RuntimeAttr? runtime_override_merge_denovo_bed
    RuntimeAttr? runtime_override_create_plots
    RuntimeAttr? runtime_override_call_outliers
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
      linux_docker = linux_docker,
      runtime_attr_override = runtime_override_make_manifests
  }

  if (length(vcfs) == 1) {
    scatter (i in range(length(contigs))) {
      call SubsetVcfByContig {
        input:
          vcf = vcfs[0],
          vcf_index = vcf_indicies[0],
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
      batch_manifest = MakeManifests.sample_manifest,
      sv_base_mini_docker = sv_base_mini_docker,
      runtime_attr_override = runtime_override_subset_samples
  }

  Array[File] to_filter_vcfs = select_first([subset_vcfs, vcfs])
  scatter (vcf in to_filter_vcfs) {
    call PreFilterVcf {
      input:
        vcf = vcf,
        max_cohort_af = max_cohort_af,
        max_gnomad_af = max_gnomad_af,
        blacklists = blacklists,
        blacklist_overlap = blacklist_overlap,
        gd_regions = gd_regions,
        samples = SubsetSamples.sample_subset,
        sv_base_mini_docker = sv_base_mini_docker,
        runtime_attr_override = runtime_override_filter_vcf
    }
  }

  Array[File] filtered_vcfs = select_all(PreFilterVcf.filtered_vcf)
  Array[File] filtered_vcf_indicies = select_all(PreFilterVcf.filtered_vcf_index)

  scatter (i in range(length(filtered_vcfs))) {
    call MatchVcfToContig {
      input:
        vcf = filtered_vcfs[i],
        vcf_index = filtered_vcf_indicies[i],
        contigs = contigs,
        sv_base_mini_docker = sv_base_mini_docker,
        runtime_attr_override = runtime_override_match_vcf_to_contig
    }
  }

  Array[File] split_vcfs = select_all(MatchVcfToContig.matched_vcf)
  Array[File] split_vcf_indicies = select_all(MatchVcfToContig.matched_vcf_index)
  Array[String] kept_contigs = select_all(MatchVcfToContig.matched_contig)

  call SubsetManifestsByBatches {
    input:
      pesr_manifest = MakeManifests.pesr_manifest,
      depth_manifest = MakeManifests.depth_manifest,
      sample_manifest = MakeManifests.sample_manifest,
      bincov_manifest = MakeManifests.bincov_manifest,
      batches = SubsetSamples.batch_subset,
      linux_docker = linux_docker,
      runtime_attr_override = runtime_override_subset_manifests
  }

  Array[String] all_pesr_vcfs = transpose(read_tsv(SubsetManifestsByBatches.subset_pesr_manifest))[1]
  Array[String] all_depth_vcfs = transpose(read_tsv(SubsetManifestsByBatches.subset_depth_manifest))[1]
  scatter (vcf in all_pesr_vcfs) {
    call BatchVcfToBed as PesrVcfToBed {
      input:
        vcf = vcf,
        samples = SubsetSamples.sample_subset,
        sv_pipeline_docker = sv_pipeline_docker,
        runtime_attr_override = runtime_override_clustered_vcf_to_bed
    }
  }

  scatter (vcf in all_depth_vcfs) {
    call BatchVcfToBed as DepthVcfToBed {
      input:
        vcf = vcf,
        samples = SubsetSamples.sample_subset,
        sv_pipeline_docker = sv_pipeline_docker,
        runtime_attr_override = runtime_override_clustered_vcf_to_bed
    }
  }

  # The BatchVcfToBed task outputs a bgzip compressed file, but the format
  # allows compressed files to be concatenated.
  call MergeBatchBedsToContigs as MergePesrBed {
    input:
      beds = PesrVcfToBed.bed,
      contigs = kept_contigs,
      sv_base_mini_docker = sv_base_mini_docker,
      runtime_attr_override = runtime_override_merge_clustered_bed
  }

  call MergeBatchBedsToContigs as MergeDepthBed {
    input:
      beds = DepthVcfToBed.bed,
      contigs = kept_contigs,
      sv_base_mini_docker = sv_base_mini_docker,
      runtime_attr_override = runtime_override_merge_clustered_bed
  }

  scatter (contig in kept_contigs) {
    call ReformatContigBed as ReformatPesrBed {
      input:
        bed = MergePesrBed.contig_beds[contig],
        contig = contig,
        type = "",
        pedigree = SubsetSamples.ped_subset,
        sv_base_mini_docker = sv_base_mini_docker,
        runtime_attr_override = runtime_override_reformat_bed
    }

    call ReformatContigBed as ReformatDepthBed {
      input:
        bed = MergeDepthBed.contig_beds[contig],
        contig = contig,
        type = "depth",
        pedigree = SubsetSamples.ped_subset,
        sv_base_mini_docker = sv_base_mini_docker,
        runtime_attr_override = runtime_override_reformat_bed
    }
  }

  # Scatter the following tasks across chromosomes: miniTasks.ScatterVCf,
  # and runDeNovo.DeNovoSVsScatter.
  scatter (i in range(length(split_vcfs))) {
    # Shards vcf
    call miniTasks.ScatterVcf as ScatterVcf {
      input:
        vcf = split_vcfs[i],
        vcf_index = split_vcf_indicies[i],
        prefix = output_prefix,
        records_per_shard = records_per_shard,
        sv_pipeline_docker = sv_pipeline_docker,
        runtime_attr_override = runtime_override_shard_vcf
    }

    # Runs the de novo calling python script on each shard and outputs a per chromosome list of de novo SVs
    call runDeNovo.DeNovoSVsScatter as GetDeNovo {
      input:
        pedigree = SubsetSamples.ped_subset,
        vcfs = ScatterVcf.shards,
        chromosome = kept_contigs[i],
        raw_proband = ReformatPesrBed.reformatted_proband_bed[i],
        raw_parents = ReformatPesrBed.reformatted_parents_bed[i],
        raw_depth_proband = ReformatDepthBed.reformatted_proband_bed[i],
        raw_depth_parents = ReformatDepthBed.reformatted_parents_bed[i],
        sample_batches = SubsetManifestsByBatches.subset_sample_manifest,
        batch_bincov_index = SubsetManifestsByBatches.subset_bincov_manifest,
        small_cnv_size = small_cnv_size,
        intermediate_cnv_size = intermediate_cnv_size,
        depth_only_size = depth_only_size,
        exclude_parent_cnv_size = exclude_parent_cnv_size,
        parents_af = parents_af,
        large_raw_overlap = large_raw_overlap,
        small_raw_overlap = small_raw_overlap,
        parents_overlap = parents_overlap,
        nearby_insertion = nearby_insertion,
        coverage_cutoff = coverage_cutoff,
        gq_min = gq_min,
        variant_interpretation_docker = variant_interpretation_docker,
        runtime_override_denovo = runtime_override_denovo,
        runtime_override_vcf_to_bed = runtime_override_vcf_to_bed,
        runtime_override_merge_bed = runtime_override_denovo_merge_bed
    }
  }

  # Merges the per-chromosome final de novo SV outputs
  call MergeDenovoBedFiles {
    input:
      bed_files = GetDeNovo.per_chromosome_final_output_file,
      variant_interpretation_docker=variant_interpretation_docker,
      runtime_attr_override = runtime_override_merge_denovo_bed
  }

  # Outputs a final callset of de novo SVs as well as outlier de novo SV calls
  call CallOutliers {
    input:
      bed_file = MergeDenovoBedFiles.merged_denovo_output,
      denovo_outlier_factor = denovo_outlier_factor,
      variant_interpretation_docker=variant_interpretation_docker,
      runtime_attr_override = runtime_override_call_outliers
  }

  # Generates plots for QC
  call CreatePlots {
    input:
      bed_file = CallOutliers.final_denovo_nonOutliers_output,
      ped_file = SubsetSamples.ped_subset,
      variant_interpretation_docker=variant_interpretation_docker,
      runtime_attr_override = runtime_override_create_plots
  }

  output {
    File cleaned_ped = SubsetSamples.ped_subset
    File final_denovo_nonOutliers = CallOutliers.final_denovo_nonOutliers_output
    File final_denovo_outliers = CallOutliers.final_denovo_outliers_output
    File final_denovo_nonOutliers_plots = CreatePlots.output_plots
    Array [File] denovo_output_annotated = GetDeNovo.per_chromosome_annotation_output_file
  }
}

###########################
# TASK DEFINITIONS
###########################

# Filter the input VCF(s).
# 1. Remove all sites that:
#    a. Are BND or CNV
#    b. Have an allele frequency or gnomAD allele frequency greater than the
#       input thresholds
#    c. Are covered by at least 50% by blacklist regions
#    d. Not covered by at least 50% by genomic disorder regions (any site
#       meeting this criteria will be kept, even if it would be
#       otherwise excluded by the previous criteria)
# 2. Remove all samples from the VCF that are not in the input samples file
task PreFilterVcf {
  input {
    File vcf
    Float max_cohort_af
    Float max_gnomad_af
    Array[File]? blacklists
    Float blacklist_overlap
    File gd_regions
    Float gd_overlap
    File samples
    String sv_base_mini_docker
    RuntimeAttr? runtime_attr_override
  }

  Float vcf_size = size(vcf, "GB")
  Float other_size = size([gd_regions, samples], "GB") + (if defined(blacklists) then size(select_first([blacklists]), "GB") else 0)
  RuntimeAttr default_attr = object {
    mem_gb: 4,
    cpu_cores: 1,
    disk_gb: ceil(vcf_size * 2 + other_size) + 16,
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

  String output_vcf = "filtered-${basename(vcf)}"
  String output_vcf_index = "${output_vcf}.tbi"
  Array[File] bl = if defined(blacklists) then select_first([blacklists]) else []

  command <<<
    set -euo pipefail

    cat2() {
      if [[ "$1" = *.gz ]]; then
        zcat "$1"
      else
        cat "$1"
      fi
    }

    bcftools query --exclude 'SVTYPE = "BND" || SVTYPE = "CNV"' \
      --format '%CHROM\t%POS0\t%END\t%ID\t%INFO/AF\t%INFO/gnomad_v4.1_sv_AF\n' \
      '~{vcf}' \
      | LC_ALL=C sort -k1,1 -k2,2n > sites.tsv
    awk -F'\t' '$5 > af || ($6 != "." && $6 > gaf){print $4}' \
      af=~{max_cohort_af} gaf=~{max_gnomad_af} sites.tsv > af_fail.list

    cut -f 1-4 sites.tsv > sites.bed
    
    # Use coverage instead of intersect
    : > blacklist_fail.list
    bl_paths='~{if defined(blacklists) then write_lines(select_first([blacklists])) else ""}'
    if [[ -n "${bl_paths:-}" ]]; then
      while read -r f; do cat2 "$f"; done < "${bl_paths}" \
        | LC_ALL=C sort -k1,1 -k2,2n > bl_merged.bed

        bedtools coverage -a sites.bed -b bl_merged.bed -sorted \
        | awk -F'\t' '$8 >= ~{blacklist_overlap} {print $4}' > blacklist_fail.list
    fi

    bedtools coverage -a sites.bed -b '~{gd_regions}' \
      | awk -F'\t' '$8 >= ~{gd_overlap} {print $4}' > gd_pass.list

    awk 'FILENAME==ARGV[1]{a[$1]} FILENAME!=ARGV[1] && !($1 in a){print}' \
      gd_pass.list af_fail.list blacklist_fail.list > exclude.list

    # we do not want to force samples because in the context of this workflow, the
    # list of samples should only contain samples that are present in the last VCF
    # in the array of input VCFs so if any other VCF does not have the exact same
    # set of samples, there is an error
    bcftools view --exclude 'SVTYPE ="BND" || SVTYPE = "CNV" || ID = @exclude.list' \
      --output-type u '~{vcf}' \
      | bcftools view --output-type z --output '~{output_vcf}' \
        --no-update --samples-file '~{samples}'

    read -r nrec _ < <(bcftools head --header 0 --records 1 "~{output_vcf}" | wc -l)
    if (( nrec == 0 )); then
      echo 'filtering VCF removed all sites'
      rm "~{output_vcf}"
    else
      bcftools index --tbi '~{output_vcf}'
    fi
  >>>

  output {
    File? filtered_vcf = output_vcf
    File? filtered_vcf_index = output_vcf_index
  }
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
    set -euo pipefail

    bcftools view --regions '~{contig}' --output-type z --output '~{contig_vcf}' \
      '~{vcf}'

    read -r nrec < <(bcftools head --header 0 --records 1 "${contig_vcf}" | wc -l)
    if (( nrec == 0 )); then
      rm "${contig_vcf}"
    fi
  >>>

  output {
    File? subset_vcf = contig_vcf
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
  String vcf_index_bn = basename(vcf_index)

  command <<<
    set -euo pipefail

    bcftools index --stats '~{vcf}' | cut -f 1 > contigs_in_vcf.list
    read -r vcf_contigs_count _ < <(wc -l contigs_in_vcf.list)
    unset -v other
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
    mv '~{vcf_index}' '~{vcf_index_bn}'
    if [[ ! -s matched_contig.list ]]; then
      printf '\n' > matched_contig.list
      rm '~{vcf_bn}'
      rm '~{vcf_index_bn}'
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
    File? matched_vcf_index = vcf_index_bn
    String? matched_contig = if read_string("matched_contig.list") == "" then null else read_string("matched_contig.list")
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
    set -euo pipefail

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
  >>>

  output {
    File ped_subset = "subset.ped"
    File sample_subset = "sample_subset.list"
    File batch_subset = "batch_subset.list"
  }
}

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
    String linux_docker
    RuntimeAttr? runtime_attr_override
  }

  Float input_size = size(batch_sample_lists, "GB")
  RuntimeAttr default_attr = object {
    mem_gb: 1,
    cpu_cores: 1,
    disk_gb: ceil(input_size) + 16,
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
    docker: linux_docker
  }

  command <<<
    set -euo pipefail

    batch_names='~{write_lines(batch_name_list)}'

    paste "${batch_names}" '~{write_lines(clustered_depth_vcf)}' > 'depth_manifest.tsv'

    manta='~{if length(manta_vcfs) > 0 then write_lines(manta_vcfs) else ""}'
    melt='~{if length(melt_vcfs) > 0 then write_lines(melt_vcfs) else ""}'
    wham='~{if length(wham_vcfs) > 0 then write_lines(wham_vcfs) else ""}'
    scramble='~{if length(scramble_vcfs) > 0 then write_lines(scramble_vcfs) else ""}'

    pesr_manifest='pesr_manifest.tsv'
    : > "${pesr_manifest}"

    if [[ "${manta}" ]]; then
        paste "${batch_names}" "${manta}" >> "${pesr_manifest}"
    fi
    if [[ "${melt}" ]]; then
        paste "${batch_names}" "${melt}" >> "${pesr_manifest}"
    fi
    if [[ "${wham}" ]]; then
        paste "${batch_names}" "${wham}" >> "${pesr_manifest}"
    fi
    if [[ "${scramble}" ]]; then
        paste "${batch_names}" "${scramble}" >> "${pesr_manifest}"
    fi

    if [[ ! -s "${pesr_manifest}" ]]; then
        printf 'at least one non-empty list of PESR evidence VCF lists should be provided\n' >&2
        exit 1
    fi

    paste "${batch_names}" '~{write_lines(batch_sample_lists)}' \
      | awk -F'\t' '{while((getline line < $2) > 0) {print $1 "\t" line}}' > 'sample_manifest.tsv'

    paste "${batch_names}" '~{write_lines(batch_bincov_matrix)}' '~{write_lines(batch_bincov_matrix_index)}' > 'bincov_manifest.tsv'
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
  }
}

# Subset the PESR, depth, sample, and bincov manifests to the given batches.
# The subset samples manifest will only contain sample IDs, not batch IDs.
task SubsetManifestsByBatches {
  input {
    File pesr_manifest
    File depth_manifest
    File sample_manifest
    File bincov_manifest
    File batches
    String linux_docker
    RuntimeAttr? runtime_attr_override
  }

  Float input_size = size([pesr_manifest, depth_manifest, sample_manifest, bincov_manifest,
   batches], "GB")

  RuntimeAttr default_attr = object {
    mem_gb: 1,
    cpu_cores: 1,
    disk_gb: ceil(input_size) + 16,
    boot_disk_gb: 8,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  runtime {
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GB"
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    docker: linux_docker
  }

  command <<<
    set -euo pipefail

    awk -F'\t' 'NR==FNR{a[$1]} NR>FNR && ($1 in a)' '~{batches}' '~{sample_manifest}' > subset_sample_manifest.tsv
    awk -F'\t' 'NR==FNR{a[$1]} NR>FNR && ($1 in a)' '~{batches}' '~{pesr_manifest}' > subset_pesr_manifest.tsv
    awk -F'\t' 'NR==FNR{a[$1]} NR>FNR && ($1 in a)' '~{batches}' '~{depth_manifest}' > subset_depth_manifest.tsv
    awk -F'\t' 'NR==FNR{a[$1]} NR>FNR && ($1 in a)' '~{batches}' '~{bincov_manifest}' > subset_bincov_manifest.tsv
  >>>

  output {
    File subset_sample_manifest = "subset_sample_manifest.tsv"
    File subset_pesr_manifest = "subset_pesr_manifest.tsv"
    File subset_depth_manifest = "subset_depth_manifest.tsv"
    File subset_bincov_manifest = "subset_bincov_manifest.tsv"
  }
}

task BatchVcfToBed {
  input {
    File vcf
    File samples
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  Float input_size = size([vcf, samples], "GB")
  RuntimeAttr default_attr = object {
    mem_gb: 1,
    cpu_cores: 1,
    disk_gb: ceil(input_size * 2) + 16,
    boot_disk_gb: 8,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  runtime {
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GB"
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    docker: sv_pipeline_docker
  }

  File vcf_basename = basename(vcf, ".vcf.gz")
  command <<<
    set -euo pipefail

    bcftools view --output-type z --output temp.vcf.gz --no-update \
      --force-samples --samples-files '~{samples}' '~{vcf}'
    svtk vcf2bed temp.vcf.gz --info SVTYPE '~{vcf_basename}.bed'
    bgzip '~{vcf_basename}.bed'
  >>>

  output {
    File bed = vcf_basename + ".bed.gz"
  }
}

# Merge batch BED files into per-contig BED files.
task MergeBatchBedsToContigs {
  input {
    Array[File]+ beds
    Array[String]+ contigs
    String sv_base_mini_docker
    RuntimeAttr? runtime_attr_override
  }

  Float input_size = size(beds, "GB")
  RuntimeAttr default_attr = object {
    mem_gb: 4,
    cpu_cores: 2,
    disk_gb: ceil(input_size * 2) + 16,
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

  command <<<
    set -euo pipefail

    mkdir splits
    awk '
    BEGIN { FS = "\t"; OFS = "\t" }
    NR == FNR {
      a[$0] = sprintf("splits/%03d.bed.gz", FNR)
      b[$0] = "bgzip -c > " a[$0]
    }
    NR > FNR && FNR == 1 {
      in_cmd = "bgzip -cd " quote($0)
      in_cmd | getline line
      if (line !~ /^#chrom/) {
        print "BED file missing header lines" > "/dev/stderr"
        exit 1
      }
      for (f in b) {
        print line | b[f]
      }
      close(in_cmd)
    }
    NR > FNR {
      in_cmd = "bgzip -cd " quote($0)
      while ((in_cmd | getline line) > 0) {
        if (line ~ /^#/) {
          continue
        }

        split(line, fields)
        if (fields[1] in a) {
          print line | b[fields[1]]
        }
      }
      close(in_cmd)
    }
    END {
      for (contig in a) {
        print contig, a[contig] > "split_beds.tsv"
      }
    }

    function quote(x) {
      gsub(/\047/, "\047\\\047\047", x)
      return x
    }' '~{write_lines(contigs)}' '~{write_lines(beds)}'
  >>>

  output {
    Map[String, File] contig_beds = read_map("split_beds.tsv")
  }
}

# Reformat a single contig BED file produced by converting the GATK-SV
# clustered VCF into a BED file with svtk vcf2bed then extracting a
# single contig. The BED file is split into proband SVs and parental SVs
# reformatted according to the requirements of the denovo_svs.py script.
task ReformatContigBed {
  input {
    File bed
    String contig
    String type
    File pedigree
    String sv_base_mini_docker
    RuntimeAttr? runtime_attr_override
  }

  Float disk_size = size(bed, "GB") * 2 + size(pedigree, "GB") + 16
  RuntimeAttr default_attr = object {
    mem_gb: 4,
    cpu_cores: 1,
    disk_gb: ceil(disk_size),
    boot_disk_gb: 10,
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

  String type_str = if type == "" then "" else ".${type}"
  String proband_bed = "${contig}.proband${type_str}.reformatted.sorted.bed.gz"
  String parents_bed = "${contig}.parents${type_str}.reformatted.sorted.bed.gz"
  command <<<
    set -euo pipefail

    # FamID ParentalID
    awk -F'\t' '!/^#/ && $1 !~ /FamID/' '~{pedigree}' \
      | awk -F'\t' '$3 && $3!=0{print $1,$3} $4 && $4!=0{print $1,$4}' OFS='\t' \
      | sort -u > parents.tsv
    # ChildID
    awk -F'\t' '!/^#/ && $1 !~ /FamID/' '~{pedigree}' \
      | awk -F'\t' 'NR==FNR{a[$2]} NR>FNR && !($2 in a){print $2}' parents.tsv - \
      | sort -u > children.list

    # CHROM start end SVTYPE samples
    bgzip -cd '~{bed}' \
      | grep -v '^#' \
      | awk -F'\t' '{split($6, a, /,/); for(i in a){print $1,$2,$3,$7,a[i]}}' \
      | awk 'BEGIN{OFS="\t"
                   out_c="sort -k1,1 -k2,2n | bgzip -c > ~{proband_bed}"
                   out_p="sort -k1,1 -k2,2n | bgzip -c > ~{parents_bed}"}
             FILENAME == ARGV[1]{c[$1]}
             FILENAME == ARGV[2]{p[$2]=$1}
             FILENAME == ARGV[3] && ($5 in p){print $1"_"$4"_"p[$5],$2,$3,$4,$5 | out_p}
             FILENAME == ARGV[3] && ($5 in c){print $1"_"$4"_"$5,$2,$3,$4,$5 | out_c}' children.list parents.tsv -
  >>>

  # Proband output should be a tab-delimited file with columns:
  # CHROM_SVTYPE_sample start end SVTYPE sample

  # Parents output is a tab-delimited file with columns:
  # CHROM_SVTYPE_FamID start end SVTYPE sample
  output {
    File reformatted_proband_bed = proband_bed
    File reformatted_parents_bed = parents_bed
  }
}

# Merge the BED files output by DeNovoSVsScatter.
task MergeDenovoBedFiles {
  input {
    Array[File] bed_files
    String variant_interpretation_docker
    RuntimeAttr? runtime_attr_override
  }

  Float bed_files_size = size(bed_files, "GB")

  RuntimeAttr default_attr = object {
    mem_gb: 3.75,
    disk_gb: ceil(10 + (bed_files_size) * 2.0),
    cpu_cores: 1,
    preemptible_tries: 2,
    max_retries: 1,
    boot_disk_gb: 8
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    docker: variant_interpretation_docker
  }

  command <<<
    set -exuo pipefail

    zcat ~{bed_files[0]} | awk 'NR<=1' > denovo.merged.bed
    zcat ~{sep=" " bed_files} | grep -v ^chrom >> denovo.merged.bed
    bgzip denovo.merged.bed
  >>>

  output {
    File merged_denovo_output = "denovo.merged.bed.gz"
  }
}

# Check for outliers in the de novo callset.
task CallOutliers {
  input {
    File bed_file
    Int denovo_outlier_factor
    String variant_interpretation_docker
    RuntimeAttr? runtime_attr_override
  }

  Float bed_files_size = size(bed_file, "GB")

  RuntimeAttr default_attr = object {
    mem_gb: 16, # 3.75
    disk_gb: ceil(10 + (bed_files_size) * 2.0),
    cpu_cores: 1,
    preemptible_tries: 2,
    max_retries: 1,
    boot_disk_gb: 8
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    docker: variant_interpretation_docker
  }

  command <<<
    set -exuo pipefail

    python /src/denovo/denovo_outliers.py --bed ~{bed_file} \
      --denovo_outlier_factor ~{denovo_outlier_factor}

    bgzip final.denovo.merged.bed
    bgzip final.denovo.merged.outliers.bed
    bgzip final.denovo.merged.allSamples.bed
  >>>

  output {
    File final_denovo_nonOutliers_output = "final.denovo.merged.bed.gz"
    File final_denovo_outliers_output = "final.denovo.merged.outliers.bed.gz"
    File final_denovo_allSamples_output = "final.denovo.merged.allSamples.bed.gz"
  }
}

# Make plots for the de novo calls.
task CreatePlots {
  input {
    File bed_file
    File ped_file
    String variant_interpretation_docker
    RuntimeAttr? runtime_attr_override
  }

  Float input_size = size(select_all([bed_file, ped_file]), "GB")

  RuntimeAttr default_attr = object {
    mem_gb: 16, # 3.75
    disk_gb: ceil(10 + input_size * 1.2),
    cpu_cores: 1,
    preemptible_tries: 2,
    max_retries: 1,
    boot_disk_gb: 8
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  command <<<
    set -exuo pipefail

    Rscript /src/denovo/denovo_sv_plots.R ~{bed_file} ~{ped_file} output_plots.pdf
  >>>

  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    docker: variant_interpretation_docker
  }

  output {
    File output_plots = "output_plots.pdf"
  }
}
