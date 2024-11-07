version 1.0

import "Structs.wdl"

workflow VCF2BED {
  input {
    File vcf
    File vcf_index
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_convert
  }

  call Convert {
    input:
      vcf = vcf,
      vcf_index = vcf_index,
      sv_pipeline_docker = sv_pipeline_docker,
      runtime_attr_override = runtime_attr_convert
  }

  output {
    File bed = Convert.bed
  }
}

task Convert {
  input {
    File vcf
    File vcf_index
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: 4,
    disk_gb: ceil(size(vcf, "GB") * 20) + 16,
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

  String bed_name = sub(sub(basename(vcf), "\\.gz$", ""), "\\.vcf$", "") + ".bed"
  command <<<
    set -o errexit
    set -o nounset
    set -o pipefail

    svtk vcf2bed -i ALL --include-filters '~{vcf}' '~{bed_name}'
    gzip "~{bed_name}"
  >>>

  output {
    File bed = bed_name + ".gz"
  }
}
