version 1.0

import "Structs.wdl"
import "CleanVcfChromosome.wdl"

workflow RescueMobileElementDeletions {
  input {
    File vcf
    String prefix
    File LINE1_reference
    File HERVK_reference
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  call CleanVcfChromosome.RescueMobileElementDeletions as rescue_med {
    input:
      vcf = vcf,
      prefix = "~{prefix}.rescue_me_dels",
      LINE1 = LINE1_reference,
      HERVK = HERVK_reference,
      sv_pipeline_docker = sv_pipeline_docker,
      runtime_attr_override = runtime_attr_override
  }

  call IndexVCF {
    input:
      vcf = rescue_med.out,
      sv_pipeline_docker = sv_pipeline_docker
  }

  output {
    File out = rescue_med.out
    File out_index = IndexVCF.index
  }
}

task IndexVCF {
  input {
    File vcf
    String sv_pipeline_docker
  }

  Int disk_size_gb = ceil(size(vcf, "GB")) + 16
  
  runtime {
    memory: "1GiB"
    disks: "local-disk ${disk_size_gb} HDD"
    cpu: 1
    preemptible: 3
    maxRetries: 1
    docker: sv_pipeline_docker
    bootDiskSizeGb: 16
  }

  command <<<
    tabix -p vcf '~{vcf}'
  >>>

  output {
    File index = "~{vcf}.tbi"
  }
}
