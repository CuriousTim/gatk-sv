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

  output {
    File out = rescue_med.out
  }
}
