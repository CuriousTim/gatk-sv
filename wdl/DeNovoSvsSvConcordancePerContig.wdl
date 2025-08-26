version 1.0

import "Structs.wdl"

workflow DeNovoSvsSvConcordancePerContig {
  input {
    Array[File]+ truth_vcfs
    File truth_samples
    Array[File]+ eval_bcfs
    Array[String]+ contigs
    File reference_dict
    String batch

    String sv_base_mini_docker
    String svconcordance_keep_all_docker
    String linux_docker

    RuntimeAttr? runtime_override_concat_raw_evidence
    RuntimeAttr? runtime_override_svconcordance
  }

  parameter_meta {
    truth_vcfs: "Clustered raw evidence VCFs (e.g. Manta, Wham, depth) for a single batch."
    truth_samples: "Samples that should be used to subset the truth VCFs."
    eval_bcfs: "BCFs to annotate with variants from the truth VCFs, split by contig."
    contigs: "Contigs present in the eval BCFs."
    reference_dict: "Sequence dictionary in the form of a '.dict' file."
    batch: "Batch ID."
    sv_base_mini_docker: "The corresponding Docker image from GATK-SV."
    svconcordance_keep_all_docker: "Docker with a build of GATK that supports the `--keep-all` option of SVConcordance."
    linux_docker: "A Linux Docker image."
    runtime_override_concat_raw_evidence: "Runtime attribute overrides for ConcatRawEvidence."
    runtime_override_svconcordance: "Runtime attribute overrides for SVConcordance."
  }

  output {
    Array[File] concordance_vcfs = SVConcordance.concordance_vcf
  }

  call ConcatRawEvidence {
    input:
      vcfs = truth_vcfs,
      contigs = contigs,
      samples = truth_samples,
      sv_base_mini_docker = sv_base_mini_docker,
      runtime_attr_override = runtime_override_concat_raw_evidence
  }

  call MakeRawEvidenceMap {
    input:
      bcf_paths = ConcatRawEvidence.concat_bcfs,
      linux_docker = linux_docker
  }

  scatter (i in range(length(contigs))) {
    String current_contig = contigs[i]
    call SVConcordance {
      input:
        truth_bcf = MakeRawEvidenceMap.contig_map[current_contig],
        eval_bcf = eval_bcfs[i],
        concordance_prefix = "${batch}-${current_contig}",
        reference_dict = reference_dict,
        svconcordance_keep_all_docker = svconcordance_keep_all_docker,
        runtime_attr_override = runtime_override_svconcordance
    }
  }
}

task ConcatRawEvidence {
  input {
    Array[File]+ vcfs
    Array[String] contigs
    File samples
    String sv_base_mini_docker
    RuntimeAttr? runtime_attr_override
  }

  parameter_meta {
    vcfs: "VCFs to concatenate."
    contigs: "Contigs into which the concatenated BCF files should be split."
    samples: "Samples, one per line, to use to subset the concatenated files."
    sv_base_mini_docker: "The corresponding Docker image from GATK-SV."
    runtime_attr_override: "Runtime attribute overrides."
  }

  output {
    Array[File] concat_bcfs = glob("merged/*.bcf")
    Array[File] concat_bcf_indexes = glob("merged/*.csi")
  }

  Float input_size = size(vcfs, "GB")
  RuntimeAttr default_attr = object {
    mem_gb: 4,
    cpu_cores: 1,
    disk_gb: ceil(input_size * 4) + 16,
    boot_disk_gb: 8,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  Float mem = select_first([runtime_attr.mem_gb, default_attr.mem_gb])

  runtime {
    memory: "${mem} GB"
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    docker: sv_base_mini_docker
  }

  Float max_sort_mem = mem * 0.8

  command <<<
    set -euxo pipefail

    LC_ALL=C sort -u '~{samples}' > samples_uniq
    contigs='~{write_lines(contigs)}'
    mkdir splits
    declare -i i=0
    while read -r src; do
      bcftools index "${src}"
      while read -r contig; do
        dest="splits/${contig}"
        if [[ ! -d "${dest}" ]]; then
          mkdir "${dest}"
        fi
        bcftools view --samples-file 'samples_uniq' --regions "${contig}" --no-update \
          --exclude 'INFO/SVTYPE == "BND"' --output-type u "${src}" \
          | bcftools view --include 'COUNT(GT="alt")' --output "${dest}/${i}.bcf" --output-type b
        bcftools index "${dest}/${i}.bcf"
      done < "${contigs}"
      i+=1
    done < '~{write_lines(vcfs)}'

    mkdir merged
    while read -r contig; do
      bcftools concat --allow-overlaps --file-list <(find "splits/${contig}" -type f -name '*.bcf') \
        --output-type u \
        | bcftools sort --max-mem '~{max_sort_mem}G' --output "merged/${contig}.bcf" --output-type b
      bcftools index "merged/${contig}.bcf"
    done < "${contigs}"
  >>>
}

task MakeRawEvidenceMap {
  input {
    Array[String] bcf_paths
    String linux_docker
  }

  parameter_meta {
    bcf_paths: "Paths to the merged raw evidence files, split by contig."
    linux_docker: "A Linux Docker image."
  }

  output {
    Map[String, File] contig_map = read_map("raw.tsv")
  }

  runtime {
    memory: "1 GB"
    cpus: 1
    disks: "local-disk 16 HDD"
    bootDiskSizeGb: 8
    preemptible: 3
    maxRetries: 1
    docker: linux_docker
  }

  command <<<
    set -euxo pipefail

    while read -r p; do
      bn="$(basename "${p}")"
      printf '%s\t%s\n' "${bn%.bcf}" "${p}" > raw.tsv
    done < '~{write_lines(bcf_paths)}'
  >>>
}

task SVConcordance {
  input {
    File truth_bcf
    File eval_bcf
    String concordance_prefix
    File reference_dict
    String svconcordance_keep_all_docker
    RuntimeAttr? runtime_attr_override
  }

  parameter_meta {
    truth_bcf: "BCF against which to match variants."
    eval_bcf: "BCF to annotate with variants matched in the truth BCF."
    concordance_prefix: "Prefix to use for the output BCF."
    reference_dict: "Sequence dictionary in the form of a '.dict' file."
    svconcordance_keep_all_docker: "Docker with a build of GATK that supports the `--keep-all` option of SVConcordance."
    runtime_attr_override: "Runtime attribute overrides."
  }

  output {
    File concordance_vcf = concordance_name
  }

  Float input_size = size([truth_bcf, eval_bcf], "GB")
  RuntimeAttr default_attr = object {
    mem_gb: 4,
    cpu_cores: 1,
    disk_gb: ceil(input_size * 3) + 16,
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
    docker: svconcordance_keep_all_docker
  }

  String concordance_name = "${concordance_prefix}.vcf.gz"

  command <<<
    set -euxo pipefail

    gatk --java-options "-Xmx3400M" SVConcordance \
      --sequence-dictionary '~{reference_dict}' \
      --eval '~{eval_bcf}' \
      --truth '~{truth_bcf}' \
      --output '~{concordance_name}'
  >>>
}
