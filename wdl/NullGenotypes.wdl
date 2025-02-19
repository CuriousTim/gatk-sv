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

    if [[ '~{vcf_index}' =~ .*\.tbi$ ]]; then
      vcf_index_ext='tbi'
    elif [[ '~{vcf_index}' =~ .*\.csi$ ]]; then
      vcf_index_ext='csi'
    else
      vcf_index_ext='index'
    fi

    bcftools query --regions '~{sep="," contigs}' --format '%CHROM\t%POS\n' '~{vcf}' > 'positions.tsv'
    if awk 'END{exit NR > 0}' 'positions.tsv'; then
      printf 'no positions found in the VCF for the given contigs\n' >&2
      printf 'no genotypes will be changed\n' >&2
      cp '~{vcf}' 'genotypes_nulled.vcf.gz'
      cp '~{vcf_index}' "genotypes_nulled.vcf.gz.${vcf_index_ext}"
      exit 0
    fi

    contigs_list='~{write_lines(contigs)}'
    printf '##fileformat=VCFv4.5\n' > 'annotations.vcf'
    awk '
    $1 !~ /^[0-9A-Za-z!#$%&+.\/:;?@\^_|~-][0-9A-Za-z!#$%&*+.\/:;=?@\^_|~-]*$/ {
      print $1 " does not conform to VCF spec for contig IDs" > "/dev/stderr"
      had_err = 1
      exit 87
    }
    {
      print "##contig=<ID=" $1 ">"
      ++contig_count
    }
    END {
      if (contig_count == 0) {
        print "at least one contig must given" > "/dev/stderr"
        exit 88
      }
    }' "${contigs_list}" >> 'annotations.vcf'
    printf '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n' >> 'annotations.vcf'

    filtered_samples='samples_to_annotate.list'
    bcftools query --list-samples '~{vcf}' | LC_ALL=C sort > 'samples_in_vcf.list'
    LC_ALL=C sort -u '~{samples_list}' > 'requested_samples.list'
    LC_ALL=C comm -12 'samples_in_vcf.list' 'requested_samples.list' > "${filtered_samples}"
    awk '
    {a[$1]}
    END{
      if (NR == 0) {
        print "no samples in the input list were found in the VCF" > "/dev/stderr"
        exit 88
      }

      printf "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT"
      for (s in a) {
        printf "\t%s", s
      }
      printf "\n"
    }' "${filtered_samples}" >> 'annotations.vcf'

    awk -F '\t' '
    NR == FNR { ++count }
    NR > FNR {
      printf "%s\t%d\t.\tA\t.\t.\t.\t.\tGT", $1, $2
      for (i = 1; i <= count; ++i) {
        printf "\t./."
      }
      printf "\n"
    }' "${filtered_samples}" 'positions.tsv' >> 'annotations.vcf'
    bgzip 'annotations.vcf'
    bcftools index 'annotations.vcf.gz'

    bcftools annotate --annotations 'annotations.vcf.gz' \
      --columns 'CHROM,POS,.FORMAT/GT' --output 'genotypes_nulled.vcf.gz' \
      --output-type z --samples-file "${filtered_samples}" '~{vcf}' \
      --pair-logic all
    bcftools index 'genotypes_nulled.vcf.gz'
  >>>

  output {
    File output_vcf = "genotypes_nulled.vcf.gz"
    File output_vcf_index = "genotypes_nulled.vcf.gz.csi"
  }
}
