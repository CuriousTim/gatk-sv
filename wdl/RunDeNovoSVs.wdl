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
    File vcf
    File vcf_index
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
    RuntimeAttr? runtime_override_subset_manifests
    RuntimeAttr? runtime_override_subset_vcf
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

  # Create text files with the paths of the raw data used by the workflow
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

  # If family_ids file is given, subset all other input files to only include the necessary batches.
  if (defined(family_ids)) {
    call SubsetManifestsByFamilies {
      input:
        pesr_manifest = MakeManifests.pesr_manifest,
        depth_manifest = MakeManifests.depth_manifest,
        sample_manifest = MakeManifests.sample_manifest,
        bincov_manifest = MakeManifests.bincov_manifest,
        family_ids = select_first([family_ids]),
        pedigree = pedigree,
        linux_docker = linux_docker,
        runtime_attr_override = runtime_override_subset_manifests
    }

    call util.SubsetVcfBySamplesList {
      input:
        vcf = vcf,
        vcf_idx = vcf_index,
        list_of_samples = SubsetManifestsByFamilies.subset_samples,
        outfile_name = output_prefix,
        keep_af = true,
        remove_private_sites = false,
        sv_base_mini_docker = sv_base_mini_docker,
        runtime_attr_override = runtime_override_subset_vcf
    }
  }

  # Makes a ped file of singletons, duos, and trios for input into the de novo script (only including families of interest)
  call CleanPed {
    input:
      ped_file = pedigree,
      vcf_input = select_first([SubsetVcfBySamplesList.vcf_subset, vcf]),
      variant_interpretation_docker = variant_interpretation_docker,
      runtime_attr_override = runtime_override_clean_ped
  }

  Array[String] all_pesr_vcfs = transpose(read_tsv(select_first([SubsetManifestsByFamilies.subset_pesr_manifest, MakeManifests.pesr_manifest])))[1]
  Array[String] all_depth_vcfs = transpose(read_tsv(select_first([SubsetManifestsByFamilies.subset_depth_manifest, MakeManifests.depth_manifest])))[1]
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

  scatter (contig in contigs) {
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

  # Scatter the following tasks across chromosomes: SubsetVcf, miniTasks.ScatterVCf,
  # and runDeNovo.DeNovoSVsScatter.
  scatter (i in range(length(contigs))) {
    # Splits vcf by chromosome
    call SubsetVcf {
      input:
        vcf = select_first([SubsetVcfBySamplesList.vcf_subset, vcf]),
        vcf_index = vcf_index,
        chromosome = contigs[i],
        sv_base_docker = sv_base_mini_docker,
        runtime_attr_override = runtime_override_subset_vcf
    }

    # Shards vcf
    call miniTasks.ScatterVcf as ScatterVcf {
      input:
        vcf = SubsetVcf.vcf_output,
        vcf_index = SubsetVcf.vcf_index_output,
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
        chromosome = contigs[i],
        raw_proband = ReformatPesrBed.reformatted_proband_bed[i],
        raw_parents = ReformatPesrBed.reformatted_parents_bed[i],
        raw_depth_proband = ReformatDepthBed.reformatted_proband_bed[i],
        raw_depth_parents = ReformatDepthBed.reformatted_parents_bed[i],
        genomic_disorder_regions = genomic_disorder_regions,
        exclude_regions = exclude_regions,
        sample_batches = MakeManifests.sample_manifest,
        batch_bincov_index = select_first([SubsetManifestsByFamilies.subset_bincov_manifest, MakeManifests.bincov_manifest]),
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
    boot_disk_gb: 10,
    preemptible_tries: 3,
    max_retries: 1,
  }
  RuntimeAttr runtime_override = select_first([runtime_attr_override, default_attr])

  Array[String] manta_vcfs = select_first([clustered_manta_vcf, []])
  Array[String] melt_vcfs = select_first([clustered_melt_vcf, []])
  Array[String] wham_vcfs = select_first([clustered_wham_vcf, []])
  Array[String] scramble_vcfs = select_first([clustered_scramble_vcf, []])

  runtime {
    memory: "${select_first([runtime_override.mem_gb, default_attr.mem_gb])} GB"
    cpu: select_first([runtime_override.cpu_cores, default_attr.cpu_cores])
    disks: "local-disk ${select_first([runtime_override.disk_gb, default_attr.disk_gb])} HDD"
    bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, default_attr.boot_disk_gb])
    preemptible: select_first([runtime_override.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_override.max_retries, default_attr.max_retries])
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

# Subset the PESR, depth, sample, and bincov manifests to batches containing the
# given families. The subset samples manifest will only contain sample IDs, not
# batch IDs.
task SubsetManifestsByFamilies {
  input {
    File pesr_manifest
    File depth_manifest
    File sample_manifest
    File bincov_manifest
    File family_ids
    File pedigree
    String linux_docker
    RuntimeAttr? runtime_attr_override
  }

  Float input_size = size([pesr_manifest, depth_manifest, sample_manifest, bincov_manifest,
   family_ids, pedigree], "GB")

  RuntimeAttr default_attr = object {
    mem_gb: 3.75,
    cpu_cores: 1,
    disk_gb: ceil(input_size) + 10,
    boot_disk_gb: 8,
    preemptible_tries: 2,
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
    set -exuo pipefail

    awk -F'\t' 'NR==FNR && NF==1{a[$1]} NR>FNR && !/^#/ && ($1 in a){print $2}' \
      '~{family_ids}' '~{pedigree}' \
      | sort -u > samples.list
    if [[ $(wc -l samples.list | awk '{print $1}') -eq 0 ]]; then
      echo 'No matching family IDs found in pedigree.' >&2
      exit 1
    fi

    awk -F'\t' 'NR==FNR && NF==1{a[$1]} NR>FNR && ($2 in a){print $1}' \
      samples.list '~{sample_manifest}' \
      | sort -u > batches.list
    if [[ $(wc -l batches.list | awk '{print $1}') -eq 0 ]]; then
      echo 'No matching batch IDs found in sample manifest.' >&2
      exit 1
    fi

    awk -F'\t' 'NR==FNR{a[$1]} NR>FNR && ($1 in a)' batches.list '~{pesr_manifest}' > subset_pesr_manifest.tsv
    awk -F'\t' 'NR==FNR{a[$1]} NR>FNR && ($1 in a)' batches.list '~{depth_manifest}' > subset_depth_manifest.tsv
    awk -F'\t' 'NR==FNR{a[$1]} NR>FNR && ($1 in a)' batches.list '~{bincov_manifest}' > subset_bincov_manifest.tsv
  >>>

  output {
    File subset_samples = "samples.list"
    File subset_pesr_manifest = "subset_pesr_manifest.tsv"
    File subset_depth_manifest = "subset_depth_manifest.tsv"
    File subset_bincov_manifest = "subset_bincov_manifest.tsv"
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
    boot_disk_gb: 10,
    preemptible_tries: 3,
    max_retries: 1,
  }
  RuntimeAttr runtime_override = select_first([runtime_attr_override, default_attr])

  runtime {
    memory: "${select_first([runtime_override.mem_gb, default_attr.mem_gb])} GB"
    cpu: select_first([runtime_override.cpu_cores, default_attr.cpu_cores])
    disks: "local-disk ${select_first([runtime_override.disk_gb, default_attr.disk_gb])} HDD"
    bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, default_attr.boot_disk_gb])
    preemptible: select_first([runtime_override.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_override.max_retries, default_attr.max_retries])
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
  RuntimeAttr runtime_override = select_first([runtime_attr_override, default_attr])

  runtime {
    memory: "${select_first([runtime_override.mem_gb, default_attr.mem_gb])} GB"
    cpu: select_first([runtime_override.cpu_cores, default_attr.cpu_cores])
    disks: "local-disk ${select_first([runtime_override.disk_gb, default_attr.disk_gb])} HDD"
    bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, default_attr.boot_disk_gb])
    preemptible: select_first([runtime_override.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_override.max_retries, default_attr.max_retries])
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

# Extract a single contig from a VCF.
task SubsetVcf {
  input {
    File vcf
    File vcf_index
    String chromosome
    String sv_base_docker
    RuntimeAttr? runtime_attr_override
  }

  Float input_size = size(vcf, "GB")

  RuntimeAttr default_attr = object {
    mem_gb: 3.75,
    disk_gb: ceil(10 + input_size * 1.5),
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
    docker: sv_base_docker
  }

  command <<<
    set -euxo pipefail

    bcftools view '~{vcf}' --regions '~{chromosome}' -O z -o '~{chromosome}.vcf.gz'
    bcftools index --tbi '~{chromosome}.vcf.gz'
  >>>

  output {
    File vcf_output = "${chromosome}.vcf.gz"
    File vcf_index_output = "${chromosome}.vcf.gz.tbi"
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
    File vcf_input
    String variant_interpretation_docker
    RuntimeAttr? runtime_attr_override
  }

  Float ped_size = size(ped_file, "GB")
  Float vcf_size = size(vcf_input, "GB")

  RuntimeAttr default_attr = object {
    mem_gb: 3.75,
    disk_gb: ceil(10 + vcf_size + ped_size * 1.5),
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

    # TODO: this script should get the name of the output file it generates as an input argument.
    # The output filename is currently hardcoded to be 'clean_ped.txt'.
    Rscript /src/denovo/clean_ped.R ~{ped_file}
    cut -f2 cleaned_ped.txt | awk 'NR > 1' > all_samples.txt
    bcftools query -l ~{vcf_input} > samples_to_include_in_ped.txt

    # Note: grep returns a non-success exit code (i.e., other than `0`) when it cannot find the
    # match in the following scripts. We do not expect it to find a match for every entry.
    # Hence, to avoid exit with unsuccessful code, we can either drop pipefail from above or use `|| true`.

    grep -w -v -f samples_to_include_in_ped.txt all_samples.txt > excluded_samples.txt || true
    grep -w -f excluded_samples.txt cleaned_ped.txt | cut -f1 | sort -u > excluded_families.txt || true
    grep -w -v -f excluded_families.txt cleaned_ped.txt > subset_cleaned_ped.txt || true
  >>>

  output {
    File cleaned_ped = "subset_cleaned_ped.txt"
  }
}
