version 1.0

import "Structs.wdl"
import "PloidyEstimation.wdl" as pe
import "WGD.wdl" as wgd

workflow CohortEvidenceQC {
  input {
    String sample_set_set_id
    Array[File] wgd_matrices
    Array[File] ploidy_matrices
    File wgd_scoring_mask
    String linux_docker
    String sv_pipeline_qc_docker
    String? chr_x
    String? chr_y
    RuntimeAttr? runtime_attr_merge_wgd
    RuntimeAttr? runtime_attr_merge_ploidy
    RuntimeAttr? runtime_attr_score_wgd
    RuntimeAttr? runtime_attr_score_ploidy
  }

  call MergeMatrices as merge_wgd {
    input:
      matrices = wgd_matrices,
      linux_docker = linux_docker,
      runtime_attr_override = runtime_attr_merge_wgd
  }

  call MergeMatrices as merge_ploidy {
    input:
      matrices = ploidy_matrices,
      linux_docker = linux_docker,
      runtime_attr_override = runtime_attr_merge_ploidy
  }

  RuntimeAttr score_wgd_default_attr = object { mem_gb: ceil(length(wgd_matrices) * 0.2) }
  call wgd.WGDScore as score_wgd {
    input:
      wgd_scoring_mask = wgd_scoring_mask,
      WGD_matrix = merge_wgd.merged_matrix,
      batch = sample_set_set_id,
      sv_pipeline_qc_docker = sv_pipeline_qc_docker,
      runtime_attr_override = select_first([runtime_attr_score_wgd, score_wgd_default_attr])
  }

  RuntimeAttr score_ploidy_default_attr = object { mem_gb: ceil(length(ploidy_matrices) * 0.2) }
  call pe.PloidyScore as score_ploidy {
    input:
      ploidy_matrix = merge_ploidy.merged_matrix,
      batch = sample_set_set_id,
      chr_x = chr_x,
      chr_y = chr_y,
      sv_pipeline_qc_docker = sv_pipeline_qc_docker,
      runtime_attr_override = select_first([runtime_attr_score_ploidy, score_ploidy_default_attr])
  }

  output {
    File cohort_wgd_dist = score_wgd.WGD_dist
    File cohort_ploidy_plots = score_ploidy.ploidy_plots
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
    disk_gb: ceil(size(matrices, "GB") * 2) + 16,
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
