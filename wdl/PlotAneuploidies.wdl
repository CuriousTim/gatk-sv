version 1.0

import "Structs.wdl"
import "PloidyEstimation.wdl" as pe

workflow PlotAneuploidies {
  input {
    File samples
    
    # Batch info
    File batch_map
    Array[String] batch_ids
    Array[String] bincov_matrices

    String? chr_x
    String? chr_y

    Int bin_size = 1000000
    Int background_size = 100

    # Runtime parameters
    String linux_docker
    String sv_base_mini_docker
    # String sv_pipeline_qc_docker

    RuntimeAttr? runtime_override_group_samples
    RuntimeAttr? runtime_override_ploidy_matrix
    RuntimeAttr? runtime_override_ploidy
  }

  call GroupSamples {
    input:
      samples = samples,
      batch_map = batch_map,
      batch_ids = batch_ids,
      bincov_matrices = bincov_matrices,
      linux_docker = linux_docker,
      runtime_override = runtime_override_group_samples
  }

  scatter (i in range(length(GroupSamples.groups))) {
    String batch_name = "batch_" + i
    call MakePloidyMatrix {
      input:
        group = GroupSamples.groups[i],
        bin_size = bin_size,
        background_size = background_size,
        output_prefix = batch_name,
        sv_base_mini_docker = sv_base_mini_docker,
        runtime_override = runtime_override_ploidy_matrix
    }

    # call PlotPloidy {
    #   input:
    #     ploidy_matrix = PloidyMatrix.ploidy_matrix
    #     batch = batch_name, 
    #     chr_x = chr_x,
    #     chr_y = chr_y,
    #     sv_pipeline_qc_docker = sv_pipeline_qc_docker,
    #     runtime_override = runtime_override_ploidy
    # }
  }

  # call MergePlots {
  #   input:
  #     plots = PlotPloidy.plots,
  #     linux_docker = linux_docker,
  #     runtime_override = runtime_override_ploidy
  # }

  output {
    Array[File] ploidy_matricies = MakePloidyMatrix.ploidy_matrix
  }
}

task GroupSamples {
  input {
    File samples
    File batch_map
    Array[String] batch_ids
    Array[String] bincov_matrices
    String linux_docker

    RuntimeAttr? runtime_override
  }

  Float input_size = size([samples, batch_map], "GB")
  RuntimeAttr runtime_default = object {
    cpu_cores: 1,
    mem_gb: 2,
    boot_disk_gb: 8,
    preemptible_tries: 3,
    max_retries: 1,
    disk_gb: ceil(input_size * 2) + 16
  }
  RuntimeAttr runtime_attr = select_first([runtime_override, runtime_default])

  runtime {
    docker: linux_docker
    cpu: select_first([runtime_attr.cpu_cores, runtime_default.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, runtime_default.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, runtime_default.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, runtime_default.boot_disk_gb])
    preemptible: select_first([runtime_attr.preemptible_tries, runtime_default.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, runtime_default.max_retries])
  }

  command <<<
    set -euo pipefail

    mkdir groups
    batch_ids='~{write_lines(batch_ids)}'
    bincov_matrices='~{write_lines(bincov_matrices)}'
    batch_n=$(cat "${batch_ids}" | wc -l)
    bincov_n=$(cat "${bincov_matrices}" | wc -l)
    if (( batch_n != bincov_n )); then
      printf 'unequal number of batches and bincov matrices (%d vs %d)\n' "${batch_n}" "${bincov_n}" >&2
      exit 1
    fi
    paste "${batch_ids}" "${bincov_matrices}" > bincovs.tsv

    awk -F'\t' '
      FILENAME == ARGV[1] {bincovs[$1] = $2}
      FILENAME == ARGV[2] {batches[$2] = $1}
      FILENAME == ARGV[3] {
        if (!($1 in batches)) {
          printf "sample \047%s\047 does not have a batch\n", $1 > "/dev/stderr"
          exit 1
        }
        q = batches[$1]
        if (!(q in bincovs)) {
          printf "batch \047s\047 does not have a bincov matrix\n", q > "/dev/stderr"
          exit 1
        }

        if (!(q in e)) {
          e[q] = sprintf("groups/%06d", i++)
          print bincovs[q] > e[q]
        }
        print $1 > e[q]
      }' bincovs.tsv '~{batch_map}' <(LC_ALL=C sort -u '~{samples}')
  >>>

  output {
    Array[File] groups = glob("groups/*")
  }
}

task MakePloidyMatrix {
  input {
    # A file in which the first line is the path to the bincov matrix for a
    # batch and the remaining lines are samples to plot (those with suspected
    # aneuploidies).
    File group
    Int bin_size
    Int background_size
    String output_prefix
    String sv_base_mini_docker
    RuntimeAttr? runtime_override
  }

  File bincov = read_lines(group)[0]
  Float input_size = size(bincov, "GB")
  RuntimeAttr runtime_default = object {
    cpu_cores: 1,
    mem_gb: 4,
    boot_disk_gb: 8,
    preemptible_tries: 3,
    max_retries: 1,
    disk_gb: ceil(input_size * 1.4) + 16
  }
  RuntimeAttr runtime_attr = select_first([runtime_override, runtime_default])

  runtime {
    docker: sv_base_mini_docker
    cpu: select_first([runtime_attr.cpu_cores, runtime_default.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, runtime_default.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, runtime_default.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, runtime_default.boot_disk_gb])
    preemptible: select_first([runtime_attr.preemptible_tries, runtime_default.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, runtime_default.max_retries])
  }

  String output_name = output_prefix + "_ploidy_matrix.bed.gz"
  command <<<
    set -euo pipefail

    bgzip -cd '~{bincov}' \
      | head -n 1 \
      | awk -F'\t' '{for (i=4; i<=NF; ++i){print $i}}' \
      | LC_ALL=C sort -u > bincov_samples.list
    awk -F'\t' 'NR>1' '~{group}' | LC_ALL=C sort -u > samples_to_plot.list
    LC_ALL=C comm -13 samples_to_plot.list bincov_samples.list \
      | shuf --head-count ~{background_size} - > bg_samples.list 
    cat samples_to_plot.list bg_samples.list > keep.list
    ls >&2

    bgzip -cd '~{bincov}' \
      | awk -F'\t' '
      NR == FNR { a[$1] }
      NR > FNR && FNR == 1 {
        sub(/^#/, "", $1)
        printf "%s\t%s\t%s", $1, $2, $3
        for (i = 4; i <= NF; ++i) {
          if ($i in a) {
            index[++j] = i
            printf "\t%s", $i
          }
        }
        print
        next
      }
      NR > FNR && FNR == 2 {
        chr = $1
        start = $2
        end = start + bin_size
      }
      NR > FNR && (chr && chr != $1) || ($2 >= end) {
        printf "%s\t%s\t%s", $1, start, end
        for (i = 1; i <= j; ++i) {
          printf "\t%s", cov[i]
          cov[i] = 0
        }
        print
        chr = $1
        start = $2
        end = start + bin_size
      }
      NR > FNR {
        for (i = 1; i <= j; ++i) {
          cov[i] += $(index[i])
        }
      }' bin_size=~{bin_size} keep.list - \
      | bgzip -c - > '~{output_name}' 
  >>>

  output {
    File ploidy_matrix = output_name
  }
}
#
# task PlotPloidy {
#   input {
#     File ploidy_matrix
#     String chr_x = "chrX"
#     String chr_y = "chrY"
#     String sv_pipeline_qc_docker
#     RuntimeAttr? runtime_override    
#   }
# 
#   Float input_size = size(ploidy_matrix, "GB")
#   RuntimeAttr runtime_default = object {
#     cpu_cores: 1,
#     mem_gb: 2,
#     boot_disk_gb: 8,
#     preemptible_tries: 3,
#     max_retries: 1,
#     disk_gb: ceil(input_size * 2) + 16
#   }
#   RuntimeAttr runtime_attr = select_first([runtime_override, runtime_default])
# 
#   runtime {
#     docker: sv_pipeline_qc_docker
#     cpu: select_first([runtime_attr.cpu_cores, runtime_default.cpu_cores])
#     memory: select_first([runtime_attr.mem_gb, runtime_default.mem_gb]) + " GiB"
#     disks: "local-disk " + select_first([runtime_attr.disk_gb, runtime_default.disk_gb]) + " HDD"
#     bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, runtime_default.boot_disk_gb])
#     preemptible: select_first([runtime_attr.preemptible_tries, runtime_default.preemptible_tries])
#     maxRetries: select_first([runtime_attr.max_retries, runtime_default.max_retries])
#   }
# 
#   command <<<
#     set -euo pipefail
#   >>>
#}
