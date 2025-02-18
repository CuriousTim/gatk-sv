version 1.0

import "Structs.wdl"

workflow NullContigGenotypes {
  input {
    File vcf
    File vcf_index
    File samples_list
    Array[String]? contigs

    String sv_base_mini_docker
    RuntimeAttr? runtime_attr_override_regenotype
  }

  call SetGenotypesToNull {
    input:
      vcf = vcf,
      vcf_index = vcf_index,
      samples_list = samples_list,
      contigs = contigs,
      runtime_docker = sv_base_mini_docker,
      runtime_attr_override = runtime_attr_override_regenotype
  }

  output {
    File genotypes_nulled_vcf = SetGenotypesToNull.output_vcf
    File genotypes_nulled_vcf_index = SetGenotypesToNull.output_vcf_index
  }
}

task SetGenotypesToNull {
  input {
    File vcf
    File vcf_index
    File samples_list
    Array[String] contigs = ["chrX", "chrY"]
    String runtime_docker 
    RuntimeAttr? runtime_attr_override
  }

  Float input_size = size([vcf, vcf_index, samples_list], "GB") 
  RuntimeAttr runtime_default = object {
    mem_gb: 4,
    disk_gb: ceil(input_size * 2 + 16),
    cpu_cores: 1,
    preemptible_tries: 3,
    max_retries: 1,
    boot_disk_gb: 16
  }

  RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])

  runtime {
    memory: "${select_first([runtime_override.mem_gb, runtime_default.mem_gb])} GB"
    disks: "local-disk ${select_first([runtime_override.disk_gb, runtime_default.disk_gb])} HDD"
    cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
    preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
    maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
    docker: runtime_docker
    bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
  }

  command <<<
    set -o errexit
    set -o nounset
    set -o pipefail

    contigs_list='~{write_lines(contigs)}'
    bcftools query --regions-file "${contigs_list}" \
      --format '%CHROM\t%POS\n' '~{vcf}' \
      | awk -F'\t' '{print $0 "\t./."}' \
      | bgzip -c > 'annotations.tsv.gz'
    tabix --begin 2 --end 2 --sequence 1 'annotations.tsv.gz'
    bcftools annotate --annotations 'annotations.tsv.gz' \
      --columns 'CHROM,POS,.FORMAT/GT' --regions-file "${contigs_list}" \
      --output 'genotypes_nulled.vcf.gz' --output-type z \
      --samples-file '~{samples_list}'
    tabix --preset vcf 'genotypes_nulled.vcf.gz'
  >>>

  output {
    File output_vcf = "genotypes_nulled.vcf.gz"
    File output_vcf_index = "genotypes_nulled.vcf.gz.tbi"
  }
}
