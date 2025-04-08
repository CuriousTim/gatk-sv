version 1.0
    
import "Structs.wdl"
import "Utils.wdl" as utils

workflow DeNovoSVsScatter {
  input {
    File pedigree
    Array[File] vcfs
    String chromosome
    File raw_proband
    File raw_parents
    File raw_depth_proband
    File raw_depth_parents
    File exclude_regions
    File sample_batches
    File batch_bincov_index
    File? genomic_disorder_regions
    File exclude_regions

    # Parameters for denovo_svs.py with default values
    Int small_cnv_size = 1000
    Int intermediate_cnv_size = 5000
    Int depth_only_size = 10000
    Int exclude_parent_cnv_size = 10000000
    # Allele frequency
    Float gnomad_af = 0.01
    Float parents_af = 0.05
    Float cohort_af = 0.05
    # Overlap parameters
    Float large_raw_overlap = 0.5
    Float small_raw_overlap = 0.5
    Float parents_overlap = 0.5
    Float blacklist_overlap = 0.5
    Int nearby_insertion = 100
    # SV quality (parents)
    Int coverage_cutoff = 10
    Float gq_min = 0
    # Other
    String gnomad_col = "gnomad_v4.1_sv_AF"
    String alt_gnomad_col = "gnomad_v4.1_sv_AF"

    String variant_interpretation_docker
    String sv_pipeline_docker

    RuntimeAttr? runtime_override_denovo
    RuntimeAttr? runtime_override_vcf_to_bed
    RuntimeAttr? runtime_override_merge_bed
  }
    
  Array[String] coverage_index_files = transpose(read_tsv(batch_bincov_index))[2]

  # Scatter genotyping over shards
  scatter (shard in vcfs) {
    call utils.VcfToBed as VcfToBed {
      input:
        vcf_file = shard,
        args = "--info ALL --include-filters",
        variant_interpretation_docker = sv_pipeline_docker,
        runtime_attr_override = runtime_override_vcf_to_bed
    }

    call RunDeNovo {
      input:
        bed_input = VcfToBed.bed_output,
        pedigree = pedigree,
        vcf = shard,
        chromosome = chromosome,
        raw_proband = raw_proband,
        raw_parents = raw_parents,
        raw_depth_proband = raw_depth_proband,
        raw_depth_parents = raw_depth_parents,
        genomic_disorder_regions = genomic_disorder_regions,
        exclude_regions = exclude_regions,
        coverage_indices = coverage_index_files,
        sample_batches = sample_batches,
        batch_bincov_index = batch_bincov_index,
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
        runtime_attr_override = runtime_override_denovo
    }
  }

  call MergeBedFiles as MergeBedFilesAnnotated {
    input:
      bed_files = RunDeNovo.annotation_output,
      chromosome = chromosome,
      variant_interpretation_docker = variant_interpretation_docker,
      runtime_attr_override = runtime_override_merge_bed
  }

  call MergeBedFiles as MergeBedFilesFinal {
    input:
      bed_files = RunDeNovo.denovo_output,
      chromosome = chromosome,
      variant_interpretation_docker = variant_interpretation_docker,
      runtime_attr_override = runtime_override_merge_bed
  }

  output {
    File per_chromosome_annotation_output_file = MergeBedFilesAnnotated.per_chromosome_denovo_output
    File per_chromosome_final_output_file = MergeBedFilesFinal.per_chromosome_denovo_output
  }
}

task RunDeNovo {
  input {
    File bed_input
    File pedigree
    File vcf
    String chromosome
    File raw_proband
    File raw_parents
    File raw_depth_proband
    File raw_depth_parents
    File? genomic_disorder_regions
    File exclude_regions
    Array[File] coverage_indices
    File batch_bincov_index
    File sample_batches

    # config parameters
    Int small_cnv_size
    Int intermediate_cnv_size
    Int depth_only_size
    Int exclude_parent_cnv_size
    Float gnomad_af
    Float parents_af
    Float cohort_af
    Float large_raw_overlap
    Float small_raw_overlap
    Float parents_overlap
    Float blacklist_overlap
    Int nearby_insertion
    Int coverage_cutoff
    Float gq_min
    String gnomad_col
    String alt_gnomad_col

    String variant_interpretation_docker

    RuntimeAttr? runtime_attr_override
  }

  Float vcf_size = size(vcf, "GB")
  Float bed_size = size(bed_input, "GB")

  RuntimeAttr default_attr = object {
    mem_gb: 16,
    disk_gb: ceil(16 + vcf_size + bed_size * 1.5),
    cpu_cores: 1,
    preemptible_tries: 2,
    max_retries: 1,
    boot_disk_gb: 8
  }
  RuntimeAttr runtime_override = select_first([runtime_attr_override, default_attr])

  String basename = basename(vcf, ".vcf.gz")
#TODO handle optional genomic disorder
  command <<<
    set -exuo pipefail

    printf "small_cnv_size: '%d'\n" '~{small_cnv_size}' > config.py
    printf "intermediate_cnv_size: '%d'\n" '~{intermediate_cnv_size}' >> config.py
    printf "depth_only_size: '%d'\n" '~{depth_only_size}' >> config.py
    printf "exclude_parent_cnv_size: '%d'\n" '~{exclude_parent_cnv_size}' >> config.py
    printf "gnomad_AF: '%0.2f'\n" '~{gnomad_af}' >> config.py
    printf "parents_AF: '%0.2f'\n" '~{parents_af}' >> config.py
    printf "cohort_AF: '%0.2f'\n" '~{cohort_af}' >> config.py
    printf "large_raw_overlap: '%0.1f'\n" '~{large_raw_overlap}' >> config.py
    printf "small_raw_overlap: '%0.1f'\n" '~{small_raw_overlap}' >> config.py
    printf "parents_overlap: '%0.1f'\n" '~{parents_overlap}' >> config.py
    printf "blacklist_overlap: '%0.1f'\n" '~{blacklist_overlap}' >> config.py
    printf "nearby_insertion: '%d'\n" '~{nearby_insertion}' >> config.py
    printf "coverage_cutoff: '%d'\n" '~{coverage_cutoff}' >> config.py
    printf "gq_min: '%0.2f'\n" '~{gq_min}' >> config.py
    printf "gnomad_col: '%s'\n" '~{gnomad_col}' >> config.py
    printf "alt_gnomad_col: '%s'\n" '~{alt_gnomad_col}' >> config.py

    # The MakeManifests tasks uses a batchID, sampleID column order while
    # denovo_svs.py expects the columns flipped.
    awk -F'\t' '{print $2"\t"$1}' '~{sample_batches}' > sample_batches_swapped.tsv

    export GCS_OAUTH_TOKEN=`gcloud auth application-default print-access-token`

    bcftools view ~{vcf} | grep -v ^## | bgzip -c > ~{basename}.noheader.vcf.gz
    python /src/denovo/denovo_svs.py \
        --bed ~{bed_input} \
        --ped ~{pedigree} \
        --vcf ~{basename}.noheader.vcf.gz \
        --out ~{basename}.annotation.bed \
        --out_de_novo ~{basename}.denovo.bed \
        --raw_proband ~{raw_proband} \
        --raw_parents ~{raw_parents} \
        --raw_depth_proband ~{raw_depth_proband} \
        --raw_depth_parents ~{raw_depth_parents} \
        --config config.py \
        --exclude_regions ~{exclude_regions} \
        --coverage ~{batch_bincov_index} \
        --sample_batches sample_batches_swapped.tsv \
        --verbose True

    bgzip ~{basename}.denovo.bed
    bgzip ~{basename}.annotation.bed
  >>>

  runtime {
    cpu: select_first([runtime_override.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_override.mem_gb, default_attr.mem_gb]) + " GB"
    disks: "local-disk " + select_first([runtime_override.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, default_attr.boot_disk_gb])
    preemptible: select_first([runtime_override.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_override.max_retries, default_attr.max_retries])
    docker: variant_interpretation_docker
  }

  output {
    File denovo_output = "${basename}.denovo.bed.gz"
    File annotation_output = "${basename}.annotation.bed.gz"
  }

}

task MergeBedFiles {
  input {
    Array[File] bed_files
    String chromosome
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
  RuntimeAttr runtime_override = select_first([runtime_attr_override, default_attr])

  command <<<
    set -exuo pipefail

    zcat ~{bed_files[0]} | awk 'NR==1' > ~{chromosome}.denovo.merged.bed
    zcat ~{sep=" " bed_files} | grep -v ^chrom >> ~{chromosome}.denovo.merged.bed
    bgzip ~{chromosome}.denovo.merged.bed
  >>>

  runtime {
    cpu: select_first([runtime_override.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_override.mem_gb, default_attr.mem_gb]) + " GB"
    disks: "local-disk " + select_first([runtime_override.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, default_attr.boot_disk_gb])
    preemptible: select_first([runtime_override.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_override.max_retries, default_attr.max_retries])
    docker: variant_interpretation_docker
  }

  output {
    File per_chromosome_denovo_output = "~{chromosome}.denovo.merged.bed.gz"
  }
}
