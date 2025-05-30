version 1.0

import "Structs.wdl"

# Plots CNV depth profiles across batches

workflow VisualizeCnvs {
  input {
    # Note vcf will be faster
    Array[File] vcfs_or_beds  # bed columns: chrom,start,end,name,svtype,samples
    Boolean pass_only = true
    String plot_prefix
    Int min_size
    Int? variants_per_shard
    String rdtest_flags

    File sample_table # TSV with sample_id, sample_set_id
    Array[String] sample_set_ids
    Array[File] median_files
    Array[File] rd_files
    Array[File] rd_file_indicies

    String linux_docker
    String sv_base_mini_docker
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_make_manifests
    RuntimeAttr? runtime_attr_format_vcf_or_bed
    RuntimeAttr? runtime_attr_merge_intervals
    RuntimeAttr? runtime_attr_merge_intervals
    RuntimeAttr? runtime_attr_subset_rd_matrix
    RuntimeAttr? runtime_attr_gather_uris
    RuntimeAttr? runtime_attr_rdtest
    RuntimeAttr? runtime_attr_merge_plot_tars
  }

  call MakeManifests {
    input:
      sample_set_ids = sample_set_ids,
      median_files = median_files,
      rd_files = rd_files,
      rd_file_indicies = rd_file_indicies,
      linux_docker = linux_docker,
      runtime_attr_override = runtime_attr_make_manifests
  }

  scatter (file in vcfs_or_beds) {
    call FormatVcfOrBed {
      input:
        vcf_or_bed = file,
        min_size = min_size,
        sample_table = sample_table,
        pass_only = pass_only,
        variants_per_shard = variants_per_shard,
        sv_pipeline_docker = sv_pipeline_docker,
        runtime_attr_override = runtime_attr_format_vcf_or_bed
    }
  }

  call MergeIntervals {
    input:
      interval_tars = FormatVcfOrBed.intervals,
      sv_base_mini_docker = sv_base_mini_docker,
      runtime_attr_override = runtime_attr_merge_intervals
  }

  scatter (intervals in MergeIntervals.merged_intervals) {
    String batch_id = basename(intervals, ".intervals")
    call SubsetRdMatrix {
      input:
        batch_id = batch_id,
        intervals = intervals,
        rd_file = MakeManifests.rd_manifest[batch_id],
        rd_file_index = MakeManifests.rd_indicies_manifest[batch_id],
        sv_base_mini_docker = sv_base_mini_docker,
        runtime_attr_override = runtime_attr_subset_rd_matrix
    }
  }

  Array[File] variant_shards = flatten(FormatVcfOrBed.variants)
  Array[File] variant_batch_shards = flatten(FormatVcfOrBed.variant_batches)
  scatter (i in range(length(variant_shards))) {
    call GatherShardUris {
      input:
        batches = variant_batch_shards[i],
        rd_batch_ids = SubsetRdMatrix.batch_id_out,
        rd_subsets = SubsetRdMatrix.rd_subset,
        rd_subset_indices = SubsetRdMatrix.rd_subset_index,
        medians_manifest = MakeManifests.medians_manifest,
        linux_docker = linux_docker,
        runtime_attr_override = runtime_attr_gather_uris
    }

    call RdTestPlot {
      input:
        variants = variant_shards[i],
        rd_files = GatherShardUris.rd_uris,
        rd_file_indicies = GatherShardUris.rd_index_uris,
        median_files = GatherShardUris.medians_uris,
        sample_table = sample_table,
        plot_prefix = plot_prefix,
        sv_pipeline_docker = sv_pipeline_docker,
        flags = rdtest_flags,
        runtime_attr_override = runtime_attr_rdtest
    }
  }

  call MergePlotTars {
    input:
      plot_tars = RdTestPlot.plots,
      plot_prefix = plot_prefix,
      sv_pipeline_docker = sv_pipeline_docker,
      runtime_attr_override = runtime_attr_merge_plot_tars
  }

  output {
    File rdtest_plots = MergePlotTars.plots
  }
}

task MakeManifests {
  input {
    Array[String] sample_set_ids
    Array[String] median_files
    Array[String] rd_files
    Array[String] rd_file_indicies
    String linux_docker
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: 1,
    disk_gb: 16,
    boot_disk_gb: 8,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: linux_docker
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }

  command <<<
    set -o errexit
    set -o nounset
    set -o pipefail

    sample_set_ids='~{write_lines(sample_set_ids)}'
    paste "${sample_set_ids}" '~{write_lines(median_files)}' > medians_manifest.tsv
    paste "${sample_set_ids}" '~{write_lines(rd_files)}' > rd_manifest.tsv
    paste "${sample_set_ids}" '~{write_lines(rd_file_indicies)}' > rd_indicies_manifest.tsv
  >>>

  output {
    Map[String, String] medians_manifest = read_map("medians_manifest.tsv")
    Map[String, String] rd_manifest = read_map("rd_manifest.tsv")
    Map[String, String] rd_indicies_manifest = read_map("rd_indicies_manifest.tsv")
  }
}

task FormatVcfOrBed {
  input {
    File vcf_or_bed
    Int min_size
    File sample_table
    Boolean pass_only
    Int variants_per_shard = 40
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: 4,
    disk_gb: ceil(size(vcf_or_bed, "GB") * 3 + size(sample_table, "GB")) + 16,
    boot_disk_gb: 16,
    preemptible_tries: 2,
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
    shopt -s failglob

    cat2() {
      if [[ "$1" == *.gz ]]; then
        zcat "$1"
      else
        cat "$1"
      fi
    }

    if (( ~{variants_per_shard} <= 0 )); then
      printf 'variants per shard must be > 0\n' >&2
      exit 1
    fi

    if [[ '~{vcf_or_bed}' == *.vcf.gz ]]; then
      bcftools view --output-type u --include '(SVTYPE == "DEL" || SVTYPE == "DUP") && SVLEN >= ~{min_size}' '~{vcf_or_bed}' \
      ~{if pass_only then "| bcftools view --output-type u --include 'FILTER ~ \"PASS\"'" else ""} \
        | svtk vcf2bed stdin raw.bed
      awk -F'\t' -v OFS='\t' '!/^#/{print $1,$2,$3,$4,$6,$5}' raw.bed > cnvs.bed
      rm raw.bed
    elif [[ '~{vcf_or_bed}' == *.bed || '~{vcf_or_bed}' == *.bed.gz ]]; then
      cat2 '~{vcf_or_bed}' \
        | awk -F'\t' -v OFS='\t' '
          !/^#/ && $3 - $2 >= ~{min_size} && ($5 == "DEL" || $5 == "DUP") {print $1,$2,$3,$4,$6,$5}' > cnvs.bed
    else
      echo "Invalid extension for input calls. Must be .vcf.gz, .bed.gz, or .bed"
      exit 1
    fi

    mkdir intervals variants
    # Split the variants into BED files of no more than {variants_per_shard}
    # variants. For each split, create a file with the list of batches in that
    # split. Create per batch files listing all the intervals needed from that
    # batch's RD matrix.
    awk -F'\t' '
      BEGIN {
        OFS="\t"
        Batches_out = sprintf("variants/%09d.bed.batches", Shard_n)
        Variants_out = sprintf("variants/%09d.bed", Shard_n)
      }
      NR==FNR {samples[$1]=$2}
      NR>FNR && ++Variants_n > Max_n {
        for (batch in Batches) {
          print batch > Batches_out
        }
        ++Shard_n
        Variants_n = 0
        delete Batches
        Batches_out = sprintf("variants/%09d.bed.batches", Shard_n)
        Variants_out = sprintf("variants/%09d.bed", Shard_n)
      }
      NR>FNR {
        print > Variants_out
        split($5, a, /,/)
        for(i in a) {
          Batches[samples[a[i]]]
          print $1,$2,$3 > ("intervals/"samples[a[i]]".intervals.tmp")
        }
      }
      END {
        if (length(Batches) > 0) {
          for (batch in Batches) {
            print batch > Batches_out
          }
        }
      }' Max_n=~{variants_per_shard} '~{sample_table}' cnvs.bed

    for f in intervals/*.intervals.tmp; do
      # Merge intervals that are close to prevent duplicates when retrieving
      # with tabix
      sort -k1,1 -k2,2n "${f}" \
        | bedtools merge -i stdin -d 101 > "${f%.tmp}"
      rm "${f}"
    done

    tar -zcf intervals.tar.gz intervals
  >>>

  output {
    # The variants files for a single split will have the same prefix so as long as the glob
    # expansion is consistently lexicographical, the output arrays should be parallel.
    Array[File] variants = glob("variants/*.bed")
    Array[File] variant_batches = glob("variants/*.bed.batches")
    File intervals = "intervals.tar.gz"
  }
}

task MergeIntervals {
  input {
    Array[File] interval_tars
    String sv_base_mini_docker
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: 2,
    disk_gb: ceil(size(interval_tars, "GB") * 4) + 16,
    boot_disk_gb: 8,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: sv_base_mini_docker
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }

  command <<<
    set -o errexit
    set -o nounset
    set -o pipefail
    shopt -s failglob

    mkdir intervals
    while read -r f; do
      tar -C "$(dirname "${f}")" -zxf "${f}"
      for g in "${f%.tar.gz}/"*.intervals; do
        cat "${g}" >> "intervals/$(basename "${g}").tmp"
      done
    done < '~{write_lines(interval_tars)}'

    for f in intervals/*.intervals.tmp; do
      sort -k1,1 -k2,2n "${f}" \
        | bedtools merge -i "stdin" -d 101 > "${f%.tmp}"
    done
  >>>

  output {
    Array[File] merged_intervals = glob("intervals/*.intervals")
  }
}

task SubsetRdMatrix {
  input {
    String batch_id
    File intervals
    File rd_file
    File rd_file_index
    String sv_base_mini_docker
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: 1,
    disk_gb: ceil(size(rd_file, "GB") * 1.1 + size(rd_file_index, "GB") + size(intervals, "GB")) + 16,
    boot_disk_gb: 8,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: sv_base_mini_docker
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }

  String subset_name = batch_id + ".RD.txt.gz"
  command <<<
    set -o errexit
    set -o nounset
    set -o pipefail

    # tabix query regions are always 1-based, inclusive but RD matrices are
    # 0-based 
    intervals_tmp="$(mktemp -p "${PWD}" tmp_XXXXXXXXX)"
    awk -F'\t' '{print $1 "\t" ($2 + 1) "\t" $3}' '~{intervals}' > "${intervals_tmp}"
    tabix --print-header --regions "${intervals_tmp}" '~{rd_file}' \
      | bgzip -c > '~{subset_name}'
    tabix --preset bed '~{subset_name}'
    rm "${intervals_tmp}"
  >>>

  output {
    String batch_id_out = batch_id
    File rd_subset = subset_name
    File rd_subset_index = subset_name + ".tbi"
  }
}

task GatherShardUris {
  input {
    File batches
    Array[String] rd_batch_ids
    Array[String] rd_subsets
    Array[String] rd_subset_indices
    Map[String, String] medians_manifest
    String linux_docker
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: 1,
    disk_gb: 16,
    boot_disk_gb: 8,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: linux_docker
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }

  command <<<
    set -o errexit
    set -o nounset
    set -o pipefail

    rd_batch_ids='~{write_lines(rd_batch_ids)}'
    rd_subsets='~{write_lines(rd_subsets)}'
    rd_subset_indices='~{write_lines(rd_subset_indices)}'
    medians_manifest='~{write_map(medians_manifest)}'

    paste "${rd_batch_ids}" "${rd_subsets}" > matrix_manifest
    paste "${rd_batch_ids}" "${rd_subset_indices}" > index_manifest
    awk -F'\t' 'NR==FNR{a[$1]} NR>FNR && ($1 in a){print $2}' '~{batches}' matrix_manifest > rd_uris
    awk -F'\t' 'NR==FNR{a[$1]} NR>FNR && ($1 in a){print $2}' '~{batches}' index_manifest > rd_index_uris
    awk -F'\t' 'NR==FNR{a[$1]} NR>FNR && ($1 in a){print $2}' '~{batches}' "${medians_manifest}" > medians_uris
  >>>

  output {
    Array[String] rd_uris = read_lines("rd_uris")
    Array[String] rd_index_uris = read_lines("rd_index_uris")
    Array[String] medians_uris = read_lines("medians_uris")
  }
}

task RdTestPlot {
  input {
    File variants
    Array[File] rd_files
    Array[File] rd_file_indicies
    Array[File] median_files
    File sample_table
    String plot_prefix
    String sv_pipeline_docker
    String flags
    RuntimeAttr? runtime_attr_override
  }

  Int variant_count = length(read_lines(variants))
  Float input_size = size(rd_files, "GB")
    + size(rd_file_indicies, "GB")
    + size(median_files, "GB")
    + size(sample_table, "GB")
    + size(variants, "GB")

  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: 8,
    disk_gb: ceil(input_size + variant_count * 0.01) + 16,
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

    mkdir median_files
    xargs -I '{}' mv '{}' median_files < '~{write_lines(median_files)}'

    mkdir rd_files
    xargs -I '{}' mv '{}' rd_files < '~{write_lines(rd_files)}'
    xargs -I '{}' mv '{}' rd_files < '~{write_lines(rd_file_indicies)}'

    mkdir rd_links
    while IFS=$'\t' read -r chr start end vid samples svtype; do
      : > variant.bed
      : > batch_ids.list
      : > merged_medians.bed

      printf '%s\t%s\t%s\t%s\t%s\t%s\n' "${chr}" "${start}" "${end}" "${vid}" "${samples}" "${svtype}" > variant.bed
      tr ',' '\n' <<< "${samples}" \
        | awk -F'\t' 'NR==FNR{a[$1]} NR>FNR && ($1 in a){print $2}' - '~{sample_table}' \
        | sort -u > batch_ids.list
      awk '{
        print ENVIRON["PWD"] "/rd_files/" $0 ".RD.txt.gz"
        print "rd_links/" $0 ".RD.txt.gz"
        print ENVIRON["PWD"] "/rd_files/" $0 ".RD.txt.gz.tbi"
        print "rd_links/" $0 ".RD.txt.gz.tbi"
      }' batch_ids.list \
        | xargs -L 2 ln -s
      awk '{print "median_files/" $0 "_medianCov.transposed.bed"}' batch_ids.list \
        | xargs paste > merged_medians.bed

      Rscript /opt/RdTest/RdTest.R \
        -b variant.bed \
        -n '~{plot_prefix}' \
        -c rd_links \
        -m merged_medians.bed \
        -p TRUE \
        ~{flags}

      find rd_links -type l -delete
    done < '~{variants}'

    mkdir '~{plot_prefix}_rd_plots'
    mv *.jpg '~{plot_prefix}_rd_plots'
    tar -czvf '~{plot_prefix}_rd_plots.tar.gz' '~{plot_prefix}_rd_plots/'
  >>>

  output {
    File plots = "${plot_prefix}_rd_plots.tar.gz"
  }
}

task MergePlotTars {
  input {
    Array[File] plot_tars
    String plot_prefix
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

    tar -czf '~{plot_prefix}_rd_plots.tar.gz' '~{plot_prefix}_rd_plots'
  >>>

  output {
    File plots = "${plot_prefix}_rd_plots.tar.gz"
  }
}
