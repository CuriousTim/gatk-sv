version 1.0

import "Structs.wdl"

# VCF processing workflow for de novo SV calling.
# 1. Subset VCF to offspring and convert to BCF
# 2. Remove uncalled SV types
# 3. Apply site level filters
# 4. Group offspring BCF by batches
workflow DeNovoSvsProcessOffspringVcf {
  input {
    File vcf
    File offspring

    Float max_cohort_af = 0.02
    Float max_gnomad_af = 0.01
    Int large_cnv_size = 1000
    Int depth_only_size = 5000
    Array[File]? exclude_regions
    Float exclude_regions_ovp = 0.5
    File gd_regions
    Float gd_regions_ovp = 0.5

    File batches
    File pedigree
    File sample_manifest

    String linux_docker
    String sv_base_mini_docker

    RuntimeAttr? runtime_override_subset_bcf_by_samples
    RuntimeAttr? runtime_override_remove_uncalled_svtypes
    RuntimeAttr? runtime_override_filter_offspring_sites
    RuntimeAttr? runtime_override_group_offspring_by_batch
  }

  parameter_meta {
    vcf: "VCF in which some subset of samples are offspring from trios."
    offspring: "Sample IDs, one per line, of trio offspring."
    max_cohort_af: "Maximum cohort allele frequency allowed."
    max_gnomad_af: "Maximum gnomAD allele frequency allowed. The value should be in the INFO field with the key 'gnomad_v4.1_sv_AF'."
    large_cnv_size: "Minimum size, in bases, of a large CNV."
    depth_only_size: "Minimum size, in bases, of a DUP that is allowed to a depth-only call."
    exclude_regions: "BED3 files of genomic regions to exclude. The files are concatenated before testing for coverage."
    exclude_regions_ovp: "Fraction of SV that must be covered by exclude regions to be dropped."
    gd_regions: "BED3 files of genomic disorder regions."
    gd_regions_ovp: "Fraction of SV that must be covered by genomic disorder regions to be dropped."
    batches: "Batches, one per line, of all samples in the VCF input."
    pedigree: "Pedigree."
    sample_manifest: "TSV with batches in the first column and samples in the second."
    linux_docker: "A Linux Docker image."
    sv_base_mini_docker: "The corresponding Docker image from GATK-SV."
    runtime_override_subset_bcf_by_samples: "Runtime attribute overrides for SubsetBcfBySamples."
    runtime_override_remove_uncalled_svtypes: "Runtime attribute overrides for RemoveUncalledSvtypes."
    runtime_override_filter_offspring_sites: "Runtime attribute overrides for FilterOffspringSites."
    runtime_override_group_offspring_by_batch: "Runtime attribute overrides for GroupOffspringByBatch."
  }

  output {
    Array[File] offspring_batch_grouped_bcfs = group_bcf_by_offspring_batch.subset_bcf
    Array[File] father_batch_grouped_bcfs = group_bcf_by_father_batch.subset_bcf
    Array[File] mother_batch_grouped_bcfs = group_bcf_by_mother_batch.subset_bcf
    Array[File] offspring_batch_offspring_ids = GroupOffspringByBatch.by_offspring
    Array[File] father_batch_father_ids = GroupOffspringByBatch.fathers
    Array[File] mother_batch_mother_ids = GroupOffspringByBatch.mothers
  }

  call SubsetBcfBySamples as make_offspring_bcf {
    input:
      bcf = vcf,
      samples = offspring,
      subset_prefix = "offspring",
      sv_base_mini_docker = sv_base_mini_docker,
      runtime_attr_override = runtime_override_subset_bcf_by_samples
  }

  call RemoveUncalledSvtypes {
    input:
      bcf = make_offspring_bcf.subset_bcf,
      sv_base_mini_docker = sv_base_mini_docker,
      runtime_attr_override = runtime_override_remove_uncalled_svtypes
  }

  call FilterOffspringSites {
    input:
      bcf = RemoveUncalledSvtypes.filtered_bcf,
      max_cohort_af = max_cohort_af,
      max_gnomad_af = max_gnomad_af,
      large_cnv_size = large_cnv_size,
      depth_only_size = depth_only_size,
      exclude_regions = exclude_regions,
      exclude_regions_ovp = exclude_regions_ovp,
      gd_regions = gd_regions,
      gd_regions_ovp = gd_regions_ovp,
      sv_base_mini_docker = sv_base_mini_docker,
      runtime_attr_override = runtime_override_filter_offspring_sites
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

  scatter (os in GroupOffspringByBatch.by_offspring) {
    call SubsetBcfBySamples as group_bcf_by_offspring_batch {
      input:
        bcf = FilterOffspringSites.filtered_bcf,
        samples = os,
        subset_prefix = basename(os),
        sv_base_mini_docker = sv_base_mini_docker,
        runtime_attr_override = runtime_override_subset_bcf_by_samples
    }
  }

  scatter (os in GroupOffspringByBatch.by_father) {
    call SubsetBcfBySamples as group_bcf_by_father_batch {
      input:
        bcf = FilterOffspringSites.filtered_bcf,
        samples = os,
        subset_prefix = basename(os),
        sv_base_mini_docker = sv_base_mini_docker,
        runtime_attr_override = runtime_override_subset_bcf_by_samples
    }
  }

  scatter (os in GroupOffspringByBatch.by_mother) {
    call SubsetBcfBySamples as group_bcf_by_mother_batch {
      input:
        bcf = FilterOffspringSites.filtered_bcf,
        samples = os,
        subset_prefix = basename(os),
        sv_base_mini_docker = sv_base_mini_docker,
        runtime_attr_override = runtime_override_subset_bcf_by_samples
    }
  }
}

# Subset a VCF/BCF to a set of samples. Output will be a BCF or an empty file.
task SubsetBcfBySamples {
  input {
    File bcf
    File samples
    String subset_prefix
    String sv_base_mini_docker
    RuntimeAttr? runtime_attr_override
  }

  parameter_meta {
    bcf: "VCF or BCF file to subset."
    samples: "Samples, one per line, to extract. If empty, the output will be an empty file."
    subset_prefix: "Prefix to use for the output."
    sv_base_mini_docker: "The corresponding Docker image from GATK-SV."
    runtime_attr_override: "Runtime attribute overrides."
  }

  output {
    File subset_bcf = subset_bcf_name
  }

  Float bcf_size = size(bcf, "GB")
  RuntimeAttr default_attr = object {
    mem_gb: 4,
    cpu_cores: 1,
    disk_gb: ceil(bcf_size * 2 + size(samples, "GB")) + 16,
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

  String subset_bcf_name = "${subset_prefix}.bcf"

  command <<<
    set -euxo pipefail

    if [[ -s '~{samples}' ]]; then
      bcftools view --no-update --samples-file '~{samples}' --output-type b \
        --output '~{subset_bcf_name}' '~{bcf}'
    else
      touch '~{subset_bcf_name}'
    fi
  >>>
}

# Remove SV types from a BCF of offspring sites that are not handled in the de novo pipeline. All
# BND and CNV sites are removed. Then the remaining sites are split into two: those that are CPX or
# CTX and those that are not. The file with CPX and CTX events will be a VCF and the file without
# those SV types will be BCF.
task RemoveUncalledSvtypes {
  input {
    File bcf
    String sv_base_mini_docker
    RuntimeAttr? runtime_attr_override
  }

  parameter_meta {
    bcf: "BCF file to filter."
    sv_base_mini_docker: "The corresponding Docker image from GATK-SV."
    runtime_attr_override: "Runtime attribute overrides."
  }

  output {
    File filtered_bcf = filtered_bcf_name
    File cpx_vcf = cpx_vcf_name
  }

  RuntimeAttr default_attr = object {
    mem_gb: 4,
    cpu_cores: 1,
    disk_gb: ceil(size(bcf, "GB") * 3)  + 16,
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

  String filtered_bcf_name = "svtypes_filtered-" + basename(bcf)
  String cpx_vcf_name = "cpx_ctx-" + basename(bcf, ".bcf") + ".vcf.gz"

  command <<<
    set -euxo pipefail

    bcftools view --exclude 'INFO/SVTYPE = "BND" || INFO/SVTYPE = "CNV"' \
      --output-type b --output tmp.bcf '~{bcf}'

    bcftools view --exclude 'INFO/SVTYPE = "CPX" || INFO/SVTYPE = "CTX"' \
      --output-type b --output '~{filtered_bcf_name}' tmp.bcf
    bcftools view --include 'INFO/SVTYPE = "CPX" || INFO/SVTYPE = "CTX"' \
      --output-type z --output '~{cpx_vcf_name}' tmp.bcf
  >>>
}

# Filter sites in an offspring BCF for potential de novos
# 1. Remove all sites that:
#    a. have an cohort or gnomAD allele frequency greater than the input
#       thresholds
#    b. are overlapped by exclude regions by a minimum of
#       `exclude_regions_ovp` fraction of the SV
#    c. small CNVs that are SR-only and don't have BOTHSIDES_SUPPORT
#    d. are depth-only DUPs and are smaller than the depth-only size threshold
#    e. are not covered by genomic disorder regions by a minimum of
#       `gd_regions_ovp` fraction of the SV (any site meeting this criteria will be
#       kept, even if it would otherwise excluded by the previous criteria)
task FilterOffspringSites {
  input {
    File bcf
    Float max_cohort_af
    Float max_gnomad_af
    Int large_cnv_size
    Int depth_only_size
    Array[File]? exclude_regions
    Float exclude_regions_ovp
    File gd_regions
    Float gd_regions_ovp
    String sv_base_mini_docker
    RuntimeAttr? runtime_attr_override
  }

  parameter_meta {
    bcf: "BCF with offspring samples."
    max_cohort_af: "Maximum cohort allele frequency allowed."
    max_gnomad_af: "Maximum gnomAD allele frequency allowed. The value should be in the INFO field with the key 'gnomad_v4.1_sv_AF'."
    large_cnv_size: "Minimum size, in bases, of a large CNV."
    depth_only_size: "Minimum size, in bases, of a DUP that is allowed to a depth-only call."
    exclude_regions: "BED3 files of genomic regions to exclude. The files are concatenated before testing for coverage."
    exclude_regions_ovp: "Fraction of SV that must be covered by exclude regions to be dropped."
    gd_regions: "BED3 files of genomic disorder regions."
    gd_regions_ovp: "Fraction of SV that must be covered by genomic disorder regions to be dropped."
    sv_base_mini_docker: "The corresponding Docker image from GATK-SV."
    runtime_attr_override: "Runtime attribute overrides."
  }

  output {
    File filtered_bcf = filtered_bcf_name
  }

  Float bcf_size = size(bcf, "GB")
  Float other_size = size(gd_regions, "GB") + (if defined(exclude_regions) then size(select_first([exclude_regions]), "GB") else 0)

  RuntimeAttr default_attr = object {
    mem_gb: 4,
    cpu_cores: 2,
    disk_gb: ceil(bcf_size * 4 + other_size) + 16,
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

  String filtered_bcf_name = "sites_filtered-${basename(bcf)}"

  command <<<
    set -euxo pipefail

    cat2() {
      if [[ "$1" = *.gz ]]; then
        zcat "$1"
      else
        cat "$1"
      fi
    }

    bcftools view --drop-genotypes --output-type b --output sites_only.bcf '~{bcf}'
    bcftools query \
      --include 'AF > ~{max_cohort_af} || (gnomad_v4.1_sv_AF != "." && gnomad_v4.1_sv_AF > ~{max_gnomad_af})' \
      --format '%ID\n' \
      sites_only.bcf > af_fail

    bcftools head sites_only.bcf | grep '^##' > headers.txt

    # Older GATK-SV VCFs have BOTHSIDES_SUPPORT in the FILTER field while
    # newer ones have it in the INFO field
    if grep -qF '##INFO=<ID=BOTHSIDES_SUPPORT,' headers.txt; then
      bothsides_filter='INFO/BOTHSIDES_SUPPORT = 1'
    elif grep -qF '##FILTER=<ID=BOTHSIDES_SUPPORT,' headers.txt; then
      bothsides_filter='FILTER ~ "BOTHSIDES_SUPPORT"'
    else
      printf 'BOTHSIDES_SUPPORT not found in BCF\n' >&2
      exit 1
    fi
    bcftools view \
      --include '(SVTYPE = "DEL" || SVTYPE = "DUP") && (EVIDENCE ~ "^RD,SR$" || EVIDENCE = "SR") && SVLEN < ~{large_cnv_size}' \
      --output-type u \
      sites_only.bcf \
      | bcftools query --exclude "${bothsides_filter}" --format '%ID\n' > bothsides_fail
    bcftools query \
      --include 'SVTYPE = "DUP" && ALGORITHMS = "depth" && SVLEN < ~{depth_only_size}' \
      --format '%ID\n' \
      sites_only.bcf > depth_only_fail

    bcftools query --format '%CHROM\t%POS0\t%END\t%ID\n' sites_only.bcf > sites.bed
    : > exclude_regions_fail
    er_paths='~{if defined(exclude_regions) then write_lines(select_first([exclude_regions])) else ""}'
    if [[ -n "${er_paths:-}" ]]; then
      while read -r f; do cat2 "${f}"; done < "${er_paths}" \
        | LC_ALL=C sort -k1,1 -k2,2n > er_merged.bed

        bedtools coverage -a sites.bed -b er_merged.bed -sorted \
          | awk -F'\t' '$8 >= ovp {print $4}' ovp=~{exclude_regions_ovp} >> exclude_regions_fail
    fi

    bedtools coverage -a sites.bed -b '~{gd_regions}' \
      | awk -F'\t' '$8 >= ~{gd_regions_ovp} {print $4}' > gd_pass

    sort -u gd_pass > whitelist
    cat af_fail bothsides_fail depth_only_fail exclude_regions_fail | sort -u > blacklist
    comm -13 whitelist blacklist > blacklist_clean

    bcftools view --exclude 'ID = "@blacklist_clean"' --output-type u \
      --output '~{filtered_bcf_name}' '~{bcf}'
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
