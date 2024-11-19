version 1.0

import "Structs.wdl"
import "MakeBincovMatrix.wdl" as mbm
import "PloidyEstimation.wdl" as pe

workflow CohortEvidenceQC {
  input {
    String sample_set_set_id
    Array[File] qc_tables
    Array[File] wgd_matrices
    Array[File] bincov_matrices
    RuntimeAttr? runtime_attr_merge_wgd
    RuntimeAttr? runtime_attr_merge_bincov
  }

  call MergeMatrices as merge_wgd {
    input:
      matrices = wgd_matrices,
      linux_docker = linux_docker,
      runtime_attr_override = runtime_attr_merge_wgd
  }

  call MergeMatrices as merge_bincov {
    input:
      matrices = bincov_matrices,
      linux_docker = linux_docker,
      runtime_attr_override = runtime_attr_merge_bincov
  }
}

task MergeMatrices {
  input {
    Array[File] matrices
    String linux_docker
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1, 
    mem_gb: 2,
    disk_gb: ceil(size(matrices, "GB"), * 2) + 16,
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
    docker: linux_docker
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }

  command <<<
    set -o errexit
    set -o nounset
    set -o pipefail

    awk '
      {
        ++i
        cmds[i] = "gzip -cd " $0
        files[i] = $0
      }
      END {
        merge_success = 1
        while (1) {
          for (j = 1; j <= i; ++j) {
            ret_code = cmds[j] | getline line   
            if (ret_code < 0) {
              print "error reading file: " files[j] > "/dev/stderr"
              merge_success = 0
              break
            }

            if (ret_code == 0) {
              if (j > 1) {
                print "unexpected end-of-file: " files[j] > "/dev/stderr"
                merge_success = 0
              }
              break
            }

            if (j == 1) {
              printf "%s", line
            } else {
              sub(/([^\t]+\t){3}/, "", line)
              printf "\t%s", line
            }
            did_print = 1
          }

          if (did_print == 1) {
            printf "\n"
          }

          if (ret_code <= 0) {
            break
          }
        }

        for (j = 1; j <= i; ++j) {
          close(cmds[j])
        }

        if (!merge_success) {
          exit 1
        }
      }' '~{write_lines(matrices)}' | gzip > merged_matrix.tsv.gz
  >>>

  output {
    File merged_matrix = "merged_matrix.tsv.gz"
  }
}
