version 1.0

import "Structs.wdl"
import "FilterGenotypes.wdl" as filter
import "MainVcfQc.wdl" as qc

# Wrapper for FilterGenotypes for using sharded VCFs

workflow FilterGenotypesSharded {
  input {
    # VCFs sharded by position
    Array[File] vcfs
    String output_prefix
    File ploidy_table

    File gq_recalibrator_model_file
    Array[String] recalibrate_gq_args = []
    Array[File] genome_tracks = []
    Float no_call_rate_cutoff = 0.05  # Set to 1 to disable NCR filtering
    Float fmax_beta = 0.4  # Recommended range [0.3, 0.5] (use lower values for higher specificity)

    # One of the following must be provided
    File? truth_json  # If given, SL cutoffs will be automatically optimized. Overrides sl_filter_args.
    String? sl_filter_args  # Explicitly set SL cutoffs. See apply_sl_filter.py for arguments.

    Int optimize_vcf_records_per_shard = 50000
    Int filter_vcf_records_per_shard = 20000

    # For MainVcfQc
    Boolean run_qc = true
    Boolean do_per_sample_qc = true
    File primary_contigs_fai
    File? ped_file
    Array[Array[String]]? site_level_comparison_datasets    # Array of two-element arrays, one per dataset, each of format [prefix, gs:// path to directory with one BED per population]
    Array[Array[String]]? sample_level_comparison_datasets  # Array of two-element arrays, one per dataset, each of format [prefix, gs:// path to per-sample tarballs]
    File? sample_renaming_tsv # File with mapping to rename per-sample benchmark sample IDs for compatibility with cohort
    String qc_bcftools_preprocessing_options = "-i 'FILTER=\"PASS\" || FILTER=\"MULTIALLELIC\"'"
    Int qc_sv_per_shard = 2500
    Int qc_samples_per_shard = 600
    RuntimeAttr? runtime_override_plot_qc_per_family
    RuntimeAttr? runtime_override_per_sample_benchmark_plot

    String linux_docker
    String gatk_docker
    String sv_base_mini_docker
    String sv_pipeline_docker
  }

  scatter (vcf in vcfs) {
    call filter.FilterGenotypes {
      input:
        vcf = vcf,
        output_prefix = output_prefix,
        ploidy_table = ploidy_table,
        gq_recalibrator_model_file = gq_recalibrator_model_file,
        recalibrate_gq_args = recalibrate_gq_args,
        genome_tracks = genome_tracks,
        no_call_rate_cutoff = no_call_rate_cutoff,
        fmax_beta = fmax_beta,
        truth_json = truth_json,
        sl_filter_args = sl_filter_args,
        optimize_vcf_records_per_shard = optimize_vcf_records_per_shard,
        filter_vcf_records_per_shard = filter_vcf_records_per_shard,
        run_qc = false,
        linux_docker = linux_docker,
        gatk_docker = gatk_docker,
        sv_base_mini_docker = sv_base_mini_docker,
        sv_pipeline_docker = sv_pipeline_docker
    }
  }

  if (run_qc) {
    call qc.MainVcfQc {
      input:
        vcfs = FilterGenotypes.filtered_vcf,
        prefix = "~{output_prefix}.filter_genotypes",
        bcftools_preprocessing_options = qc_bcftools_preprocessing_options,
        ped_file = ped_file,
        do_per_sample_qc = do_per_sample_qc,
        primary_contigs_fai = primary_contigs_fai,
        sample_renaming_tsv = sample_renaming_tsv,
        sv_per_shard = qc_sv_per_shard,
        runtime_override_per_sample_benchmark_plot = runtime_override_per_sample_benchmark_plot,
        runtime_override_plot_qc_per_family = runtime_override_plot_qc_per_family,
        samples_per_shard = qc_samples_per_shard,
        sv_pipeline_qc_docker = sv_pipeline_docker,
        sv_base_mini_docker = sv_base_mini_docker,
        sv_pipeline_docker = sv_pipeline_docker
    }
  }

  output {
    Array[File] filtered_vcf = FilterGenotypes.filtered_vcf
    Array[File] filtered_vcf_index = FilterGenotypes.filtered_vcf_index
    File? main_vcf_qc_tarball = MainVcfQc.sv_vcf_qc_output

    # For optional analysis
    Array[File?] vcf_optimization_table = FilterGenotypes.vcf_optimization_table
    Array[File?] sl_cutoff_qc_tarball = FilterGenotypes.sl_cutoff_qc_tarball
    Array[File] unfiltered_recalibrated_vcf = FilterGenotypes.unfiltered_recalibrated_vcf
    Array[File] unfiltered_recalibrated_vcf_index = FilterGenotypes.unfiltered_recalibrated_vcf_index
  }
}
