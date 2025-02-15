version 1.0

import "Structs.wdl"
import "CleanVcfChromosome.wdl"

workflow RescueMobileElementDeletions {
  input {
    File vcf
    File vcf_index
    File contig_list
    String prefix
    Boolean? has_end2
    File LINE1_reference
    File HERVK_reference
    String sv_pipeline_docker
    String sv_base_mini_docker
    RuntimeAttr? runtime_attr_override
  }

  Array[String] contigs = read_lines(contig_list)

  scatter (co in contigs) {
    call GetVcfContig {
      input:
        vcf = vcf,
        vcf_index = vcf_index,
        contig = co,
        sv_base_mini_docker = sv_base_mini_docker,
    }

    call CleanVcfChromosome.RescueMobileElementDeletions as rescue_med {
      input:
        vcf = GetVcfContig.vcf_contig,
        prefix = co,
        has_end2 = has_end2,
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

task GetVcfContig {
  input {
    File vcf
    File vcf_index
    String contig
    String sv_base_mini_docker
  }

  Int disk_size_gb = ceil(size(vcf, "GB") * 2.5 + size(vcf, "GB")) + 16

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

    bcftools view --regions '~{contig}' --output-type z '~{vcf}' \
      --output '~{contig}.vcf.gz'
  >>>

  output {
    File vcf_contig = "${contig}.vcf.gz"
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
    bcftools index --tbi '~{prefix}-me_dels.vcf.gz'
  >>>

  output {
    File merged_vcf = "${prefix}-me_dels.vcf.gz"
    File merged_vcf_index = "${prefix}-me_dels.vcf.gz.tbi"
  }
}
