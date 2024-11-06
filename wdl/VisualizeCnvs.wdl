version 1.0

import "Structs.wdl"

# Plots CNV depth profiles across batches

workflow VisualizeCnvs {
  input{
    # Note vcf will be faster
    File vcf_or_bed  # bed columns: chrom,start,end,name,svtype,samples

    String prefix
    Array[File] median_files
    Array[File] rd_files
    File ped_file
    Int min_size
    Int? variants_per_shard
    String flags
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_shard_variants
    RuntimeAttr? runtime_attr_merge_medians
    RuntimeAttr? runtime_attr_rdtest
    RuntimeAttr? runtime_attr_merge_plot_tars
  }

  scatter (file in rd_files) {
    File rd_file_indexes = file + ".tbi"
  }

  call ShardVariants {
    input:
      vcf_or_bed = vcf_or_bed,
      min_size = min_size,
      variants_per_shard = variants_per_shard,
      sv_pipeline_docker = sv_pipeline_docker,
      runtime_attr_override = runtime_attr_shard_variants
  }

  call MergeMedians {
    input:
      median_files = median_files,
      sv_pipeline_docker = sv_pipeline_docker,
      runtime_attr_override = runtime_attr_merge_medians
  }

  scatter (shard in ShardVariants.shards) {
    call RdTestPlot {
      input:
        variants = shard,
        rd_files = rd_files,
        rd_file_indexes = rd_file_indexes,
        medians = MergeMedians.medians,
        sample_ids = ShardVariants.sample_ids,
        ped_file = ped_file,
        prefix = prefix,
        sv_pipeline_docker = sv_pipeline_docker,
        flags = flags,
        runtime_attr_override = runtime_attr_rdtest
    }
  }

  call MergePlotTars {
    input:
      plot_tars = RdTestPlot.plots,
      prefix = prefix,
      sv_pipeline_docker = sv_pipeline_docker,
      runtime_attr_override = runtime_attr_merge_plot_tars
  }

  output{
    File rdtest_plots = MergePlotTars.plots
  }
}

task MergeMedians {
  input {
    Array[File] median_files
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: 1,
    disk_gb: ceil(size(median_files, "GB") * 2) + 16,
    boot_disk_gb: 16,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: sv_pipeline_docker
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }

  command <<<
    set -o errexit
    set -o nounset
    set -o pipefail

    paste ~{sep = " " median_files} > median_file.txt
  >>>

  output {
    File medians = "median_file.txt"
  }
}

task ShardVariants {
  input {
    File vcf_or_bed
    Int min_size
    Int variants_per_shard = 50
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: 4,
    disk_gb: ceil(size(vcf_or_bed, "GB") * 10) + 16,
    boot_disk_gb: 16,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: sv_pipeline_docker
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }

  command <<<
    set -o errexit
    set -o nounset
    set -o pipefail

    cat2() {
      if [[ "$1" == *.gz ]]; then
        zcat "$1"
      else
        cat "$1"
      fi
    }
    if [[ '~{vcf_or_bed}' == *.vcf.gz ]]; then
      bcftools view --include '(SVTYPE == "DEL" || SVTYPE == "DUP") && SVLEN >= ~{min_size}' '~{vcf_or_bed}' \
        | svtk vcf2bed stdin raw.bed
      awk -F'\t' -v OFS='\t' '{print $1,$2,$3,$4,$6,$5}' raw.bed \
        | LC_ALL=C sort -k1,1 -k2,2n > cnvs.bed
      rm raw.bed
      bcftools query --list-samples '~{vcf_or_bed}' > samples.txt
    elif [[ '~{vcf_or_bed}' == *.bed || '~{vcf_or_bed}' == *.bed.gz ]]; then
      cat2 '~{vcf_or_bed}' \
        | awk -F'\t' -v OFS='\t' '
          BEGIN {cmd="sort -u > samples.txt"}
          !/^#/ && $3 - $2 >= ~{min_size} && ($5 == "DEL" || $5 == "DUP") {print $1,$2,$3,$4,$6,$5; a=$6; gsub(/,/, "\n", a); print a | cmd}' \
        | LC_ALL=C sort -k1,1 -k2,2n > cnvs.bed
    else
      echo "Invalid extension for input calls. Must be .vcf.gz, .bed.gz, or .bed"
      exit 1
    fi

    mkdir variant_splits
    split -l ~{variants_per_shard} cnvs.bed variant_splits/cnvs_
  >>>

  output {
    Array[File] shards = glob("variant_splits/cnvs_*")
    File sample_ids = "samples.txt"
  }
}

task RdTestPlot {
  input {
    File variants
    Array[File] rd_files
    Array[File] rd_file_indexes
    File medians
    File sample_ids
    File ped_file
    String prefix
    String sv_pipeline_docker
    String flags
    RuntimeAttr? runtime_attr_override
  }

  Int variant_count = length(read_lines(variants))
  Int sample_count = length(read_lines(sample_ids))

  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: ceil(0.001 * sample_count),
    disk_gb: ceil((0.01 * variant_count) + size(variants, "GB") + size(rd_files, "GB")) + 32,
    boot_disk_gb: 16,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: sv_pipeline_docker
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }

  command <<<
    set -o errexit
    set -o nounset
    set -o pipefail

    bedtools merge -i '~{variants}' | cut -f1-3 > merged.bed

    mkdir rd_matrices
    while read -r rd; do
      mv "${rd}" "${rd}.tbi" rd_matrices
    done < '~{write_lines(rd_files)}'

    Rscript /opt/RdTest/RdTestV2.R \
      -b '~{variants}' \
      -n ~{prefix} \
      -x rd_matrices \
      -m '~{medians}' \
      -f '~{ped_file}' \
      -p TRUE \
      -w '~{sample_ids}' \
      ~{flags}

    mkdir ~{prefix}_rd_plots
    mv *.jpg ~{prefix}_rd_plots
    tar -czvf ~{prefix}_rd_plots.tar.gz ~{prefix}_rd_plots/
  >>>

  output {
    File plots = "${prefix}_rd_plots.tar.gz"
  }
}

task MergePlotTars {
  input {
    Array[File] plot_tars
    String prefix
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: 1,
    disk_gb: ceil(size(plot_tars, "GB") * 3) + 16,
    boot_disk_gb: 16,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: sv_pipeline_docker
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }

  command <<<
    set -o errexit
    set -o nounset
    set -o pipefail

    while read -r pt; do
      tar -xzf "${pt}"
    done < '~{write_lines(plot_tars)}'

    tar -czf ~{prefix}_rd_plots.tar.gz ~{prefix}_rd_plots
  >>>

  output {
    File plots = "${prefix}_rd_plots.tar.gz"
  }
}
