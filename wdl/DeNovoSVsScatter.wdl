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
    File sample_batches
    File batch_bincov_index

    # Parameters for denovo_svs.py with default values
    Int small_cnv_size = 1000
    Int intermediate_cnv_size = 5000
    Int depth_only_size = 10000
    Int exclude_parent_cnv_size = 10000000
    # Allele frequency
    Float parents_af = 0.05
    # Overlap parameters
    Float large_raw_overlap = 0.5
    Float small_raw_overlap = 0.5
    Float parents_overlap = 0.5
    Int nearby_insertion = 100
    # SV quality (parents)
    Int coverage_cutoff = 10
    Int gq_min = 0

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
        coverage_indices = coverage_index_files,
        sample_batches = sample_batches,
        batch_bincov_index = batch_bincov_index,
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
    Array[File] coverage_indices
    File batch_bincov_index
    File sample_batches

    # config parameters
    Int small_cnv_size
    Int intermediate_cnv_size
    Int depth_only_size
    Int exclude_parent_cnv_size
    Float parents_af
    Float large_raw_overlap
    Float small_raw_overlap
    Float parents_overlap
    Int nearby_insertion
    Int coverage_cutoff
    Int gq_min

    String variant_interpretation_docker

    RuntimeAttr? runtime_attr_override
  }

  Float vcf_size = size(vcf, "GB")
  Float bed_size = size(bed_input, "GB")
  Float raw_files_size = size([raw_proband, raw_parents, raw_depth_proband, raw_depth_parents], "GB")

  RuntimeAttr default_attr = object {
    mem_gb: 16,
    disk_gb: ceil(16 + vcf_size + raw_files_size + bed_size * 3),
    cpu_cores: 1,
    preemptible_tries: 2,
    max_retries: 1,
    boot_disk_gb: 8
  }
  RuntimeAttr runtime_override = select_first([runtime_attr_override, default_attr])

  String basename = basename(vcf, ".vcf.gz")
  command <<<
    set -exuo pipefail

    # The MakeManifests task uses a batchID, sampleID column order while
    # denovo_svs.py expects the columns flipped.
    awk -F'\t' '{print $2"\t"$1}' '~{sample_batches}' > sample_batches_swapped.tsv

    export GCS_OAUTH_TOKEN=`gcloud auth application-default print-access-token`

    bcftools view ~{vcf} | grep -v ^## | bgzip -c > ~{basename}.noheader.vcf.gz
    python /src/denovo/denovo_svs.py \
        ~{bed_input} \
        ~{pedigree} \
        ~{basename}.noheader.vcf.gz \
        ~{basename}.annotation.bed \
        ~{basename}.denovo.bed \
        ~{raw_proband} \
        ~{raw_parents} \
        ~{raw_depth_proband} \
        ~{raw_depth_parents} \
        ~{batch_bincov_index} \
        sample_batches_swapped.tsv \
        --small-cnv-size ~{small_cnv_size} \
        --med-cnv-size ~{intermediate_cnv_size} \
        --depth-only-size ~{depth_only_size} \
        --max-parent-cnv-size ~{exclude_parent_cnv_size} \
        --parents-af ~{parents_af} \
        --large-raw-overlap ~{large_raw_overlap} \
        --small-raw-overlap ~{small_raw_overlap} \
        --parents-overlap ~{parents_overlap} \
        --nearby-insertion ~{nearby_insertion} \
        --coverage-cutoff ~{coverage_cutoff} \
        --gq-min ~{gq_min} \
        --verbose

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
