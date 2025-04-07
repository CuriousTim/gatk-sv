version 1.0

###########################
# IMPORT TOOLS
###########################

import "Structs.wdl"
import "TasksMakeCohortVcf.wdl" as miniTasks
import "DeNovoSVsScatter.wdl" as runDeNovo
import "Utils.wdl" as util

###########################
# MAIN WORKFLOW DEFINITION
###########################

workflow DeNovoSV {
  input {
    # Core input files
    File pedigree
    # One family ID per line to call de novo in subset of families
    File? family_ids

    # Either a single VCF or an array of VCFs with each one containing a single
    # contig. In the case of a single VCF, it is expected that all the contigs
    # in the input contigs are present.
    Array[File] vcfs
    Array[File] vcf_indicies
    File? genomic_disorder_regions
    File exclude_regions

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
    Float? blacklist_overlap
    Int? nearby_insertion
    Int? coverage_cutoff
    Float? gq_min
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
    RuntimeAttr? runtime_override_match_vcf_to_contig
    RuntimeAttr? runtime_override_subset_ped
    RuntimeAttr? runtime_override_subset_vcf
    RuntimeAttr? runtime_override_samples_from_vcf
    RuntimeAttr? runtime_override_get_sample_batches
    RuntimeAttr? runtime_override_subset_manifests
    RuntimeAttr? runtime_override_clean_ped
    RuntimeAttr? runtime_override_clustered_vcf_to_bed
    RuntimeAttr? runtime_override_merge_clustered_bed
    RuntimeAttr? runtime_override_contig_from_bed
    RuntimeAttr? runtime_override_reformat_bed
    RuntimeAttr? runtime_override_subset_vcf
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
  }

  if (length(vcfs) > 1) {
    scatter (i in range(length(vcfs))) {
      call MatchVcfToContig {
        input:
          vcf = vcfs[i],
          vcf_index = vcf_indicies[i],
          contigs = contigs,
          sv_base_mini_docker = sv_base_mini_docker,
          runtime_attr_override = runtime_override_match_vcf_to_contig
      }
    }
  }

  Array[File] split_vcfs = select_all(select_first([SubsetVcfByContig.subset_vcf, MatchVcfToContig.matched_vcf]))
  Array[File] split_vcf_indicies = select_all(select_first([SubsetVcfByContig.subset_vcf_index, MatchVcfToContig.matched_vcf_index]))
  Array[String] kept_contigs = select_all(if defined(MatchVcfToContig.matched_contig) then select_first([MatchVcfToContig.matched_contig]) else contigs)

  # If family_ids file is given, subset the pedigree and VCFs to those families.
  if (defined(family_ids)) {
    call SubsetPedByFams {
      input:
        ped = pedigree,
        fams = select_first([family_ids]),
        linux_docker = linux_docker,
        runtime_attr_override = runtime_override_subset_ped
    }

    scatter (i in range(length(split_vcfs))) {
      call util.SubsetVcfBySamplesList as SubsetVcfWithPed {
        input:
          vcf = split_vcfs[i],
          vcf_idx = split_vcf_indicies[i],
          list_of_samples = SubsetPedByFams.samples,
          outfile_name = "${kept_contigs[i]}-subset.vcf.gz",
          sv_base_mini_docker = sv_base_mini_docker,
          runtime_attr_override = runtime_override_subset_vcf
      }
    }
  }

  File final_ped = select_first([SubsetPedByFams.subset_ped, pedigree])
  Array[File] final_vcfs = select_first([SubsetVcfWithPed.vcf_subset, split_vcfs])
  Array[File] final_vcf_indicies = select_first([SubsetVcfWithPed.vcf_subset_index, split_vcf_indicies])

  scatter (i in range(length(final_vcfs))) {
    call util.GetSampleIdsFromVcf {
      input:
        vcf = final_vcfs[i],
        sv_base_mini_docker = sv_base_mini_docker,
        runtime_attr_override = runtime_override_samples_from_vcf
    }
  }

  call GetSampleBatches {
    input:
      samples = GetSampleIdsFromVcf.out_file,
      sample_batch_map = MakeManifests.sample_manifest,
      runtime_attr_override = runtime_override_get_sample_batches
  }

  call SubsetManifestsByBatches {
    input:
      pesr_manifest = MakeManifests.pesr_manifest,
      depth_manifest = MakeManifests.depth_manifest,
      sample_manifest = MakeManifests.sample_manifest,
      bincov_manifest = MakeManifests.bincov_manifest,
      batches = GetSampleBatches.batches,
      runtime_attr_override = runtime_override_subset_manifests
  }

  # Makes a ped file of singletons, duos, and trios for input into the de novo script (only including families of interest)
  call CleanPed {
    input:
      ped_file = final_ped,
      variant_interpretation_docker = variant_interpretation_docker,
      runtime_attr_override = runtime_override_clean_ped
  }

  Array[String] all_pesr_vcfs = transpose(read_tsv(SubsetManifestsByBatches.subset_pesr_manifest))[1]
  Array[String] all_depth_vcfs = transpose(read_tsv(SubsetManifestsByBatches.subset_depth_manifest))[1]
  scatter (vcf in all_pesr_vcfs) {
    call util.VcfToBed as PesrVcfToBed {
      input:
        vcf_file = vcf,
        args = "--info SVTYPE",
        variant_interpretation_docker = sv_pipeline_docker,
        runtime_attr_override = runtime_override_clustered_vcf_to_bed
    }
  }

  scatter (vcf in all_depth_vcfs) {
    call util.VcfToBed as DepthVcfToBed {
      input:
        vcf_file = vcf,
        args = "--info SVTYPE",
        variant_interpretation_docker = sv_pipeline_docker,
        runtime_attr_override = runtime_override_clustered_vcf_to_bed
    }
  }

  call miniTasks.CatUncompressedFiles as MergePesrBed {
    input:
      shards = PesrVcfToBed.bed_output,
      sv_base_mini_docker = sv_base_mini_docker,
      runtime_attr_override = runtime_override_merge_clustered_bed
  }

  call miniTasks.CatUncompressedFiles as MergeDepthBed {
    input:
      shards = DepthVcfToBed.bed_output,
      sv_base_mini_docker = sv_base_mini_docker,
      runtime_attr_override = runtime_override_merge_clustered_bed
  }

  scatter (contig in kept_contigs) {
    call GetContigFromBed as GetContigFromPesrBed {
      input:
        bed = MergePesrBed.outfile,
        contig = contig,
        sv_base_mini_docker = sv_base_mini_docker,
        runtime_attr_override = runtime_override_contig_from_bed
    }

    call GetContigFromBed as GetContigFromDepthBed {
      input:
        bed = MergeDepthBed.outfile,
        contig = contig,
        sv_base_mini_docker = sv_base_mini_docker,
        runtime_attr_override = runtime_override_contig_from_bed
    }

    call ReformatContigBed as ReformatPesrBed {
      input:
        bed = GetContigFromPesrBed.contig_bed,
        contig = contig,
        type = "",
        pedigree = pedigree,
        sv_base_mini_docker = sv_base_mini_docker,
        runtime_attr_override = runtime_override_reformat_bed
    }

    call ReformatContigBed as ReformatDepthBed {
      input:
        bed = GetContigFromDepthBed.contig_bed,
        contig = contig,
        type = "depth",
        pedigree = pedigree,
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
        pedigree = CleanPed.cleaned_ped,
        vcfs = ScatterVcf.shards,
        chromosome = kept_contigs[i],
        raw_proband = ReformatPesrBed.reformatted_proband_bed[i],
        raw_parents = ReformatPesrBed.reformatted_parents_bed[i],
        raw_depth_proband = ReformatDepthBed.reformatted_proband_bed[i],
        raw_depth_parents = ReformatDepthBed.reformatted_parents_bed[i],
        genomic_disorder_regions = genomic_disorder_regions,
        exclude_regions = exclude_regions,
        sample_batches = SubsetManifestsByBatches.subset_sample_manifest,
        batch_bincov_index = SubsetManifestsByBatches.subset_bincov_manifest,
        small_cnv_size = small_cnv_size,
        intermediate_cnv_size = intermediate_cnv_size,
        depth_only_size = depth_only_size,
        exclude_parent_cnv_size = exclude_parent_cnv_size,
        gnomad_af = gnomad_af,
        parents_af = parents_af,
        cohort_af = cohort_af,
        large_raw_overlap = large_raw_overlap,
        small_raw_overlap = small_raw_overlap,
        parents_overlap = parents_overlap,
        blacklist_overlap = blacklist_overlap,
        nearby_insertion = nearby_insertion,
        coverage_cutoff = coverage_cutoff,
        gq_min = gq_min,
        gnomad_col = gnomad_col,
        alt_gnomad_col = alt_gnomad_col,
        variant_interpretation_docker = variant_interpretation_docker,
        runtime_override_denovo = runtime_override_denovo,
        runtime_override_vcf_to_bed = runtime_override_vcf_to_bed,
        runtime_override_merge_bed = runtime_override_denovo_merge_bed
    }
  }

  # Merges the per chromosome final de novo SV outputs
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
      variant_interpretation_docker=variant_interpretation_docker,
      runtime_attr_override = runtime_override_call_outliers
  }

  # Generates plots for QC
  call CreatePlots {
    input:
      bed_file = CallOutliers.final_denovo_nonOutliers_output,
      ped_file = pedigree,
      variant_interpretation_docker=variant_interpretation_docker,
      runtime_attr_override = runtime_override_create_plots
  }

  output {
    File cleaned_ped = CleanPed.cleaned_ped
    File final_denovo_nonOutliers = CallOutliers.final_denovo_nonOutliers_output
    File final_denovo_outliers = CallOutliers.final_denovo_outliers_output
    File final_denovo_nonOutliers_plots = CreatePlots.output_plots
    Array [File] denovo_output_annotated = GetDeNovo.per_chromosome_annotation_output_file
  }
}

###########################
# TASK DEFINITIONS
###########################

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
  String contig_vcf_index = "${contig}.vcf.gz.tbi"

  command <<<
    set -euo pipefail

    bcftools view --regions '~{contig}' --output-type z --output '~{contig_vcf}' \
      '~{vcf}'
    bcftools index --tbi '~{contig_vcf}'
  >>>

  output {
    File subset_vcf = contig_vcf
    File subset_vcf_index = contig_vcf_index  
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

  command <<<
    set -euo pipefail

    bcftools index --stats '~{vcf}' | cut -f 1 > contigs_in_vcf.list
    wc -l contigs_in_vcf.list | read -r vcf_contigs_count other
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
    if [[ ! -s matched_contig.list ]]; then
      rm matched_contig.list
      rm '~{vcf}'
      rm '~{vcf_index}'
    fi
  >>>

  output {
    File? matched_vcf = vcf
    File? matched_vcf_index = vcf_index
    String? matched_contig = read_string("matched_contig.list")
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

# Subset pedigree by family IDs.
task SubsetPedByFams {
  input {
    File ped
    File fams
    String linux_docker
    RuntimeAttr? runtime_attr_override
  }

  Float input_size = size([ped, fams], "GB")
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

    # The pedigree format required by the pipeline has a header
    awk -F'\t' 'NR==FNR{a[$1]} NR>FNR && FNR==1{print; next} NR>FNR && ($1 in a)' \
      '~{fams}' '~{ped}' > subset.ped
    awk 'FNR>1' | cut -f 2 | LC_ALL=C sort -u > subset_samples.list
  >>>

  output {
    File subset_ped = "subset.ped"
    File samples = "subset_samples.list"
  }
}

# Get the batches associated with the given samples.
task GetSampleBatches {
  input {
    Array[File] samples
    File sample_batch_map
    String linux_docker
    RuntimeAttr? runtime_attr_override
  }

  Float input_size = size(samples, "GB") + size(sample_batch_map, "GB")
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

    cat '~{write_lines(samples)}' \
      | xargs -L -- cat \
      | awk -F'\t' 'NR==FNR{a[$1]} NR>FNR && ($2 in a){print $1}' \
      | LC_ALL=C sort -u > batches.list
  >>>

  output {
    File batches = "batches.list"
  }
}

# Extract a single contig from a BED file.
task GetContigFromBed {
  input {
    File bed
    String contig
    String sv_base_mini_docker
    RuntimeAttr? runtime_attr_override
  }

  Float input_size = size(bed, "GB")
  RuntimeAttr default_attr = object {
    mem_gb: 1,
    cpu_cores: 1,
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
    bgzip -cd '~{bed}' \
      | awk -F'\t' '$1 == "~{contig}"' \
      | bgzip -c > '~{contig}.bed.gz'
  >>>

  output {
    File contig_bed = "${contig}.bed.gz"
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
  command <<<
    # FamID ParentalID
    awk -F'\t' 'BEGIN{OFS="\t"} /^#/ || $1 ~ /FamID/{next} {print $1,$3; print $1,$4}' '~{pedigree}' \
      | sort -u > parents.tsv
    # ChildID
    awk -F'\t' 'BEGIN{OFS="\t"} /^#/ || $1 ~ /FamID/{next} NR==FNR{a[$2]} NR>FNR && !($2 in a){print $2}' parents.tsv '~{pedigree}' \
      | sort -u > children.list

    # CHROM start end SVTYPE samples
    bgzip -cd '~{bed}' \
      | awk -F'\t' 'BEGIN{OFS="\t"} /^#/{next} {split($6, a, /,/); for(i in a){print $1,$2,$3,$7,a[i]}}' \
      | awk 'BEGIN{OFS="\t"
                   out_c="sort -k1,1 -k2,2n | bgzip -c > ~{contig}.proband~{type_str}.reformatted.bed.gz"
                   out_p="sort -k1,1 -k2,2n | bgzip -c > ~{contig}.parents~{type_str}.reformatted.bed.gz"}
             FILENAME == ARGV[1]{c[$1]}
             FILENAME == ARGV[2]{p[$2]=$1}
             FILENAME == ARGV[3] && ($5 in p){print $1"_"$4"_"p[$5],$2,$3,$4,$5 | out_p}
             FILENAME == ARGV[3] && ($5 in c}{print $1"_"$4"_"$5,$2,$3,$4,$5 | out_c}' children.list parents.tsv -
  >>>

  # Proband output should be a tab-delimited file with columns:
  # CHROM_SVTYPE_sample start end SVTYPE sample

  # Parents output is a tab-delimited file with columns:
  # CHROM_SVTYPE_FamID start end SVTYPE sample
  output {
    File reformatted_proband_bed = "${contig}.proband${type_str}.reformatted.sorted.bed.gz"
    File reformatted_parents_bed = "${contig}.proband${type_str}.reformatted.sorted.bed.gz"
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

    python /src/denovo/denovo_outliers.py --bed ~{bed_file}

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

# Organize the pedigree for DeNovoSVsScatter.
task CleanPed {
  input {
    File ped_file
    String variant_interpretation_docker
    RuntimeAttr? runtime_attr_override
  }

  Float ped_size = size(ped_file, "GB")

  RuntimeAttr default_attr = object {
    mem_gb: 3.75,
    disk_gb: ceil(10 + ped_size * 1.5),
    cpu_cores: 1,
    preemptible_tries: 3,
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

    # TODO: this script should get the name of the output file it generates as an input argument.
    # The output filename is currently hardcoded to be 'clean_ped.txt'.
    Rscript /src/denovo/clean_ped.R ~{ped_file}
  >>>

  output {
    File cleaned_ped = "clean_ped.txt"
  }
}
