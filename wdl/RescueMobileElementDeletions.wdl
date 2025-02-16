version 1.0

import "Structs.wdl"
import "CleanVcfChromosome.wdl"

workflow RescueMobileElementDeletions {
  input {
    File vcf
    File vcf_index
    File contig_list
    String prefix
    Boolean? add_end2
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
        add_end2 = add_end2,
        contig = co,
        sv_base_mini_docker = sv_base_mini_docker,
    }

    call CleanVcfChromosome.RescueMobileElementDeletions as rescue_med {
      input:
        vcf = GetVcfContig.vcf_contig,
        prefix = co,
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
    Boolean add_end2 = false
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

    if [[ ~{true="true" false="false" add_end2} == 'true' ]]; then
      bcftools view --regions '~{contig}' --output-type v '~{vcf} \
        | awk -F'\t' '$0 ~ /^##/ {print; if ($0 ~ /##INFO=<ID=END2/) { a = 1 } next}
            $0 ~ /^#CHROM/ {if (!a) { print "##INFO=<ID=END2,Number=1,Type=Integer,Description=\"Position of breakpoint on CHR2\">" } print; next}
            $8 ~ /SVTYPE=BND/ && $8 ~ /STRANDS=\+-/ {
              chr2 = ""
              svlen = 0
              match($8, /CHR2=[^;]+/)
              if (RSTART) {
                chr2 = substr($8, RSTART + 5, RLENGTH - 5)
              }
              match($8, /SVLEN=[^;]+/)
              if (RSTART) {
                svlen = substr($8, RSTART + 6, RLENGTH - 6)
              }

              if ($1 == chr2 && svlen) {
                sub(/END2=[^;]+;?/, "", $8)
                $8 = $8 ";END2=" ($2 + svlen)
              }
          } 1' OFS='\t' - \
        | bgzip -c > '~{contig}.vcf.gz'
    else
      bcftools view --regions '~{contig}' --output-type z --output '~{contig}.vcf.gz' '~{vcf}'
    fi
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
