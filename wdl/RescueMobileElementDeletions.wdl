version 1.0

import "Structs.wdl"
import "CleanVcfChromosome.wdl"

workflow RescueMobileElementDeletions {
  input {
    File vcf
    File vcf_index
    File contig_list
    String prefix
    File LINE1_reference
    File HERVK_reference
    String sv_pipeline_docker
    String sv_base_mini_docker
    RuntimeAttr? runtime_attr_override
  }

  Array[String] contigs = transpose(read_tsv(contig_list))[0]

  call SplitVCF {
    input:
      vcf = vcf,
      vcf_index = vcf_index,
      contigs = contigs,
      sv_base_mini_docker = sv_base_mini_docker
  }

  scatter (i in range(length(SplitVCF.splits))) {
    call CleanVcfChromosome.RescueMobileElementDeletions as rescue_med {
      input:
        vcf = SplitVCF.splits[i],
        prefix = i,
        LINE1 = LINE1_reference,
        HERVK = HERVK_reference,
        sv_pipeline_docker = sv_pipeline_docker,
        runtime_attr_override = runtime_attr_override
    }
  }

  call MergeVCFs {
    input:
      vcfs = rescue_med.out,
      prefix = prefix,
      sv_base_mini_docker = sv_base_mini_docker
  }

  output {
    File out = MergeVCFs.merged_vcf
    File out_index = MergeVCFs.merged_vcf_index
  }
}

task SplitVCF {
  input {
    File vcf
    File vcf_index
    Array[String] contigs
    String sv_base_mini_docker
  }

  Int disk_size_gb = ceil(size(vcf, "GB") * 2.5) + 16

  runtime {
    memory: "1GiB"
    disks: "local-disk ${disk_size_gb} HDD"
    cpu: 1
    preemptible: 3
    maxRetries: 1
    docker: sv_base_mini_docker
    bootDiskSizeGb: 16
  }

  command <<<
    set -o errexit
    set -o nounset
    set -o pipefail

    mkdir splits
    while read -r co; do
      bcftools view --regions "${co}" --output-type z '~{vcf}' > "splits/${co}.vcf.gz"
    done < '~{write_lines(contigs)}'
  >>>

  output {
    Array[File] splits = glob("splits/*.vcf.gz")
  }
}

task MergeVCFs {
  input {
    Array[File] vcfs
    String prefix
    String sv_base_mini_docker
  }

  Int disk_size_gb = ceil(size(vcfs, "GB") * 2.5) + 16

  runtime {
    memory: "1GiB"
    disks: "local-disk ${disk_size_gb} HDD"
    cpu: 1
    preemptible: 3
    maxRetries: 1
    docker: sv_base_mini_docker
    bootDiskSizeGb: 16
  }

  command <<<
    set -o errexit
    set -o nounset
    set -o pipefail

    bcftools concat --file-list '~{write_lines(vcfs)}' \
      --naive --output-type z > '~{prefix}-me_dels.vcf.gz'
  >>>

  output {
    File merged_vcf = "${prefix}-me_dels.vcf.gz"
    File merged_vcf_index = "${prefix}-me_dels.vcf.gz.tbi"
  }
}
