version 1.0

import "Structs.wdl"

# Group the offspring BCF by offspring batch, father batch, and mother batch.
workflow DeNovoSvsGroupOffspringBcf {
  input {
    File bcf
    Array[File] by_offspring_batch
    Array[File] by_father_batch
    Array[File] by_mother_batch
    Array[String] batch_ids

    String sv_base_mini_docker

    RuntimeAttr? runtime_override_group_bcf_by_family_batch
  }

  parameter_meta {
    bcf: "BCF in which all samples are offspring."
    by_offspring_batch: "Offspring samples grouped by offspring batch."
    by_father_batch: "Offspring samples grouped by father batch."
    by_mother_batch: "Offspring samples grouped by mother batch."
    batch_ids: "Batch IDs, parallel to offspring sample groups."
    sv_base_mini_docker: "The corresponding Docker image from GATK-SV."
    runtime_override_group_bcf_by_family_batch: "Runtime attribute overrides for GroupBcfByFamilyBatch."
  }

  output {
    Array[File] offspring_batch_grouped_bcfs = GroupBcfByFamilyBatch.by_offspring_bcf
    Array[File] father_batch_grouped_bcfs = GroupBcfByFamilyBatch.by_father_bcf
    Array[File] mother_batch_grouped_bcfs = GroupBcfByFamilyBatch.by_mother_bcf
  }

  scatter (i in range(length(batch_ids))) {
    call GroupBcfByFamilyBatch {
      input:
        bcf = bcf,
        by_offspring = by_offspring_batch[i],
        by_father = by_father_batch[i],
        by_mother = by_mother_batch[i],
        group_prefix = batch_ids[i],
        sv_base_mini_docker = sv_base_mini_docker,
        runtime_attr_override = runtime_override_group_bcf_by_family_batch
    }
  }
}

# Subset a VCF or BCF offspring samples into three different groups. The
# task should be scattered over the batches and the samples in each
# group are the offspring samples in the current batch, the samples whose
# fathers are in the current batch, and the samples whose mothers are in
# the current batch.
task GroupBcfByFamilyBatch {
  input {
    File bcf
    File by_offspring
    File by_father
    File by_mother
    String group_prefix
    String sv_base_mini_docker
    RuntimeAttr? runtime_attr_override
  }

  parameter_meta {
    bcf: "VCF or BCF file to group."
    by_offspring: "Samples, one per line, in the offspring group."
    by_father: "Samples, one per line, in the father group."
    by_mother: "Samples, one per line, in the mother group."
    group_prefix: "Prefix to use for the outputs."
    sv_base_mini_docker: "The corresponding Docker image from GATK-SV."
    runtime_attr_override: "Runtime attribute overrides."
  }

  output {
    File by_offspring_bcf = by_offspring_name
    File by_father_bcf = by_father_name
    File by_mother_bcf = by_mother_name
  }

  Float bcf_size = size(bcf, "GB")
  Float samples_size = size([by_offspring, by_father, by_mother], "GB")
  RuntimeAttr default_attr = object {
    mem_gb: 4,
    cpu_cores: 1,
    disk_gb: ceil(bcf_size * 4 + samples_size) + 16,
    boot_disk_gb: 8,
    preemptible_tries: 3,
    max_retries: 1,
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  runtime {
    memory: "${select_first([runtime_attr.mem_gb, default_attr.mem_gb])} GB"
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    disks: "local-disk ${select_first([runtime_attr.disk_gb, default_attr.disk_gb])} HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    docker: sv_base_mini_docker
  }

  String by_offspring_name = "${group_prefix}-by_offspring.bcf"
  String by_father_name = "${group_prefix}-by_father.bcf"
  String by_mother_name = "${group_prefix}-by_mother.bcf"

  command <<<
    set -euxo pipefail

    : > groups.tsv
    if [[ -s '~{by_offspring}' ]]; then
      awk '{print $1"\t-\t"f}' f='~{by_offspring_name}' '~{by_offspring}' >> groups.tsv
    else
      touch '~{by_offspring_name}'
    fi
    if [[ -s '~{by_father}' ]]; then
      awk '{print $1"\t-\t"f}' f='~{by_father_name}' '~{by_father}' >> groups.tsv
    else
      touch '~{by_father_name}'
    fi
    if [[ -s '~{by_mother}' ]]; then
      awk '{print $1"\t-\t"f}' f='~{by_mother_name}' '~{by_mother}' >> groups.tsv
    else
      touch '~{by_mother_name}'
    fi
    bcftools +split --groups-file groups.tsv --output-type b '~{bcf}'
  >>>
}
