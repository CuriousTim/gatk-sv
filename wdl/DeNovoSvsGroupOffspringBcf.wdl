version 1.0

import "Structs.wdl"

# Group the offspring BCF by offspring batch, father batch, and mother batch.
workflow DeNovoSvsGroupOffspringBcf {
  input {
    File bcf
    File offspring

    File batches
    File pedigree
    File sample_manifest

    String linux_docker
    String sv_base_mini_docker

    RuntimeAttr? runtime_override_group_offspring_by_batch
    RuntimeAttr? runtime_override_group_bcf_by_family_batch
  }

  parameter_meta {
    bcf: "BCF in which all samples are offspring."
    offspring: "Sample IDs, one per line, of offspring."
    batches: "Batches, one per line, of all samples in the VCF input."
    pedigree: "Pedigree."
    sample_manifest: "TSV with batches in the first column and samples in the second."
    linux_docker: "A Linux Docker image."
    sv_base_mini_docker: "The corresponding Docker image from GATK-SV."
    runtime_override_group_offspring_by_batch: "Runtime attribute overrides for GroupOffspringByBatch."
    runtime_override_group_bcf_by_family_batch: "Runtime attribute overrides for GroupBcfByFamilyBatch."
  }

  output {
    Array[File] offspring_batch_grouped_bcfs = GroupBcfByFamilyBatch.by_offspring_bcf
    Array[File] father_batch_grouped_bcfs = GroupBcfByFamilyBatch.by_father_bcf
    Array[File] mother_batch_grouped_bcfs = GroupBcfByFamilyBatch.by_mother_bcf
    Array[File] offspring_batch_offspring_ids = GroupOffspringByBatch.by_offspring
    Array[File] father_batch_father_ids = GroupOffspringByBatch.fathers
    Array[File] mother_batch_mother_ids = GroupOffspringByBatch.mothers
  }

  call GroupOffspringByBatch {
    input:
      offspring = offspring,
      batches = batches,
      pedigree = pedigree,
      sample_manifest = sample_manifest,
      linux_docker = linux_docker,
      runtime_attr_override = runtime_override_group_offspring_by_batch
  }

  Array[String] batch_ids = read_lines(batches)
  scatter (i in range(length(batch_ids))) {
    call GroupBcfByFamilyBatch {
      input:
        bcf = bcf,
        by_offspring = GroupOffspringByBatch.by_offspring[i],
        by_father = GroupOffspringByBatch.by_father[i],
        by_mother = GroupOffspringByBatch.by_mother[i],
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
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GB"
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
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

    if [[ -s '~{by_offspring}' ]]; then
      bcftools view --no-update --samples-file '~{by_offspring}' --output-type b \
        --output '~{by_offspring_name}' '~{bcf}'
    else
      touch '~{by_offspring_name}'
    fi
    if [[ -s '~{by_father}' ]]; then
      bcftools view --no-update --samples-file '~{by_father}' --output-type b \
        --output '~{by_father_name}' '~{bcf}'
    else
      touch '~{by_father_name}'
    fi
    if [[ -s '~{by_mother}' ]]; then
      bcftools view --no-update --samples-file '~{by_mother}' --output-type b \
        --output '~{by_mother_name}' '~{bcf}'
    else
      touch '~{by_mother_name}'
    fi
  >>>
}

# Group offspring IDs by different batches.
# The offspring IDs are grouped once by offspring batch, once by father batch and once by mother
# batch. The parental IDs are also grouped by their batches.
task GroupOffspringByBatch {
  input {
    File offspring
    File batches
    File pedigree
    File sample_manifest
    String linux_docker
    RuntimeAttr? runtime_attr_override
  }

  parameter_meta {
    offspring: "Offspring sample IDs, one per line."
    batches: "IDs of batches being processed, one per line."
    pedigree: "Pedigree."
    sample_manifest: "TSV of all samples (including parents) and batches being processed. First column is batch ID, second is sample ID."
    linux_docker: "A Linux Docker image."
    runtime_attr_override: "Runtime attribute overrides."
  }

  output {
    Array[File] by_offspring = glob("by_offspring/*")
    Array[File] by_father = glob("by_father/*")
    Array[File] fathers = glob("fathers/*")
    Array[File] by_mother = glob("by_mother/*")
    Array[File] mothers = glob("mothers/*")
  }

  Float input_size = size([offspring, batches, pedigree, sample_manifest], "GB")

  RuntimeAttr default_attr = object {
    mem_gb: 2,
    cpu_cores: 1,
    disk_gb: ceil(input_size * 2) + 16,
    boot_disk_gb: 8,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  runtime {
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GB"
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    docker: linux_docker
  }

  command <<<
    set -euxo pipefail

    mkdir by_offspring by_father fathers by_mother mothers
    # We need to unconditionally create a file for each batch for transpose to work
    while read -r b; do
      touch "by_offspring/${b}" "by_father/${b}" "fathers/${b}" "by_mother/${b}" "mothers/${b}"
    done < <(sort -u '~{batches}')

    awk -F'\t' 'FILENAME == ARGV[1] {a[$2]=$1}
                FILENAME == ARGV[2] {b[$2]=$3;c[$2]=$4}
                FILENAME == ARGV[3] {
                  print $1 > ("by_offspring/" a[$1])
                  print $1 > ("by_father/" a[b[$1]])
                  print b[$1] > ("fathers/" a[b[$1]])
                  print $1 > ("by_mother/" a[c[$1]])
                  print c[$1] > ("mothers/" a[c[$1]])
                }' '~{sample_manifest}' '~{pedigree}' '~{offspring}'
  >>>
}
