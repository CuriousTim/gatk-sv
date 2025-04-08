version 1.0

import "CollectQcVcfWide.wdl" as vcfwideqc
import "CollectQcPerSample.wdl" as persample
import "TasksMakeCohortVcf.wdl" as MiniTasks
import "Utils.wdl" as Utils

# Main workflow to perform comprehensive quality control (QC) on
# an SV VCF output by GATK-SV
workflow PlotQcPerFamily {
  input {
    Array[File] vcfs  # Option to provide a single GATK-SV VCF or an array of position-sharded SV VCFs. Must be indexed
    Boolean vcf_format_has_cn = true
    String? bcftools_preprocessing_options
    File ped_file
    Int max_trios = 1000
    String prefix
    Int sv_per_shard
    Int samples_per_shard
    Boolean do_per_sample_qc = true
    File primary_contigs_fai
    Int? random_seed
    Int? max_gq  # Max GQ for plotting. Default = 99, ie. GQ is on a scale of [0,99]. Prior to CleanVcf, use 999

    String linux_docker
    String sv_base_mini_docker
    String sv_pipeline_docker
    String sv_pipeline_qc_docker

    # overrides for local tasks
    RuntimeAttr? runtime_override_plot_qc_per_family

    # overrides for MiniTasks or Utils
    RuntimeAttr? runtime_override_subset_vcf
    RuntimeAttr? runtime_override_merge_vcfwide_stat_shards
    RuntimeAttr? runtime_override_merge_vcf_2_bed

    # overrides for CollectQcVcfWide
    RuntimeAttr? runtime_override_preprocess_vcf
    RuntimeAttr? runtime_override_collect_sharded_vcf_stats
    RuntimeAttr? runtime_override_svtk_vcf_2_bed
    RuntimeAttr? runtime_override_scatter_vcf
    RuntimeAttr? runtime_override_merge_subvcf_stat_shards

    # overrides for CollectQcPerSample
    RuntimeAttr? runtime_override_collect_vids_per_sample
    RuntimeAttr? runtime_override_split_samples_list
    RuntimeAttr? runtime_override_tar_shard_vid_lists
    RuntimeAttr? runtime_override_merge_sharded_per_sample_vid_lists
  }

  Array[String] contigs = transpose(read_tsv(primary_contigs_fai))[0]

  # Scatter raw variant data collection per chromosome
  scatter (i in range(length(contigs))) {
    # Collect VCF-wide summary stats
    call vcfwideqc.CollectQcVcfWide {
      input:
        vcfs=vcfs,
        contig=contigs[i],
        sv_per_shard=sv_per_shard,
        bcftools_preprocessing_options=bcftools_preprocessing_options,
        prefix="~{prefix}.~{contigs[i]}",
        sv_base_mini_docker=sv_base_mini_docker,
        sv_pipeline_docker=sv_pipeline_docker,
        runtime_override_preprocess_vcf=runtime_override_preprocess_vcf,
        runtime_override_collect_sharded_vcf_stats=runtime_override_collect_sharded_vcf_stats,
        runtime_override_svtk_vcf_2_bed=runtime_override_svtk_vcf_2_bed,
        runtime_override_scatter_vcf=runtime_override_scatter_vcf,
        runtime_override_merge_subvcf_stat_shards=runtime_override_merge_subvcf_stat_shards,
        runtime_override_merge_svtk_vcf_2_bed=runtime_override_merge_vcf_2_bed
    }
  }

  # Merge shards into single VCF stats file
  call MiniTasks.ConcatBeds as MergeVcfWideStats {
    input:
      shard_bed_files=CollectQcVcfWide.vcf_stats,
      prefix="~{prefix}.VCF_sites.stats",
      index_output=true,
      sv_base_mini_docker=sv_base_mini_docker,
      runtime_attr_override=runtime_override_merge_vcfwide_stat_shards
  }

  # Shard sample list
  call MiniTasks.SplitUncompressed as SplitSamplesList {
    input:
      whole_file=CollectQcVcfWide.samples_list[0],
      lines_per_shard=samples_per_shard,
      shard_prefix="~{prefix}.list_shard.",
      sv_pipeline_docker=sv_pipeline_docker,
      runtime_attr_override=runtime_override_split_samples_list
  }

  # Collect per-sample VID lists for each sample shard
  scatter ( shard in SplitSamplesList.shards ) {
    call persample.CollectQcPerSample {
      input:
        vcfs=vcfs,
        vcf_format_has_cn=vcf_format_has_cn,
        samples_list=shard,
        prefix=prefix,
        sv_base_mini_docker=sv_base_mini_docker,
        sv_pipeline_docker=sv_pipeline_docker,
        runtime_override_collect_vids_per_sample=runtime_override_collect_vids_per_sample,
        runtime_override_merge_sharded_per_sample_vid_lists=runtime_override_merge_sharded_per_sample_vid_lists
    }
  }

  # Merge all VID lists into single output directory and tar it
  call TarShardVidLists {
    input:
      in_tarballs=CollectQcPerSample.vid_lists,
      folder_name="~{prefix}_perSample_VIDs_merged",
      tarball_prefix="~{prefix}_perSample_VIDs",
      sv_base_mini_docker=sv_base_mini_docker,
      runtime_attr_override=runtime_override_tar_shard_vid_lists
  }

  Int max_gq_ = select_first([max_gq, 99])

  # Plot per-family stats if .ped file provided as input
  call PlotQcPerFamily {
    input:
      vcf_stats=MergeVcfWideStats.merged_bed_file,
      samples_list=CollectQcVcfWide.samples_list[0],
      ped_file=ped_file,
      max_trios=max_trios,
      per_sample_tarball=TarShardVidLists.vid_lists,
      prefix=prefix,
      max_gq=max_gq_,
      sv_pipeline_qc_docker=sv_pipeline_qc_docker,
      runtime_attr_override=runtime_override_plot_qc_per_family
  }

  # Final output
  output {
    File qc_per_family_plots = PlotQcPerFamily.perFamily_plots_tarball
    File denovos = PlotQcPerFamily.denovo_dump
  }
}

# Task to merge VID lists across shards
task TarShardVidLists {
  input {
    Array[File] in_tarballs
    String? folder_name
    String? tarball_prefix
    String sv_base_mini_docker
    RuntimeAttr? runtime_attr_override
  }

  String tar_folder_name = select_first([folder_name, "merged"])
  String outfile_name = select_first([tarball_prefix, tar_folder_name]) + ".tar.gz"

  # Since the input files are often/always compressed themselves, assume compression factor for tarring is 1.0
  Float input_size = size(in_tarballs, "GB")
  RuntimeAttr runtime_default = object {
    mem_gb: 2.0,
    disk_gb: ceil(10.0 + input_size * 2.0),
    cpu_cores: 1,
    preemptible_tries: 1,
    max_retries: 1,
    boot_disk_gb: 10
  }
  RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])
  runtime {
    memory: "~{select_first([runtime_override.mem_gb, runtime_default.mem_gb])} GB"
    disks: "local-disk ~{select_first([runtime_override.disk_gb, runtime_default.disk_gb])} HDD"
    cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
    preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
    maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
    docker: sv_base_mini_docker
    bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
  }

  command <<<
    # Create final output directory
    mkdir "~{tar_folder_name}"

    while read tarball_path; do
      tar -xzvf "$tarball_path" --directory ~{tar_folder_name}/
    done < ~{write_lines(in_tarballs)}

    # Compress final output directory
    tar -czvf "~{outfile_name}" "~{tar_folder_name}"
  >>>

  output {
    File vid_lists = outfile_name
  }
}


# Plot per-family stats
task PlotQcPerFamily {
  input {
    File vcf_stats
    File samples_list
    File ped_file
    File per_sample_tarball
    Int max_trios
    Int? random_seed
    String prefix
    Int max_gq
    String sv_pipeline_qc_docker
    RuntimeAttr? runtime_attr_override
  }
  Int random_seed_ = select_first([random_seed, 0])
  RuntimeAttr runtime_default = object {
    mem_gb: 15,
    disk_gb: 100,
    cpu_cores: 1,
    preemptible_tries: 1,
    max_retries: 1,
    boot_disk_gb: 10
  }
  RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])
  runtime {
    memory: "~{select_first([runtime_override.mem_gb, runtime_default.mem_gb])} GiB"
    disks: "local-disk ~{select_first([runtime_override.disk_gb, runtime_default.disk_gb])} HDD"
    cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
    preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
    maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
    docker: sv_pipeline_qc_docker
    bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
  }

  command <<<
    set -eu -o pipefail
    
    # Clean fam file
    /opt/sv-pipeline/scripts/vcf_qc/cleanFamFile.sh \
      ~{samples_list} \
      ~{ped_file} \
      ~{prefix}.cleaned.fam
    rm ~{ped_file} ~{samples_list}

    # Only run if any families remain after cleaning
    n_fams=$( grep -Ev "^#" ~{prefix}.cleaned.fam | wc -l ) || true
    echo -e "DETECTED $n_fams FAMILIES"
    if [ $n_fams -gt 0 ]; then

      # Make per-sample directory
      mkdir ~{prefix}_perSample/

      # Untar per-sample VID lists
      mkdir tmp_untar/
      tar -xvzf ~{per_sample_tarball} \
        --directory tmp_untar/
      for FILE in $( find tmp_untar/ -name "*.VIDs_genotypes.txt.gz" ); do
        mv -v $FILE ~{prefix}_perSample/
      done

      # Subset fam file, if optioned
      n_trios=$( grep -Ev "^#" ~{prefix}.cleaned.fam \
                 | awk '{ if ($2 != "0" && $2 != "." && \
                              $3 != "0" && $3 != "." && \
                              $4 != "0" && $4 != ".") print $0 }' \
                 | wc -l )
      echo -e "DETECTED $n_trios COMPLETE TRIOS"
      if [ $n_trios -gt ~{max_trios} ]; then
        grep -E '^#' ~{prefix}.cleaned.fam > fam_header.txt
        grep -Ev "^#" ~{prefix}.cleaned.fam \
        | awk '{ if ($2 != "0" && $2 != "." && \
                     $3 != "0" && $3 != "." && \
                     $4 != "0" && $4 != ".") print $0 }' \
        | sort -R --random-source <( yes ~{random_seed_} ) \
        > cleaned.shuffled.fam
        awk -v max_trios="~{max_trios}" 'NR <= max_trios' cleaned.shuffled.fam \
        | cat fam_header.txt - \
        > cleaned.subset.fam
        echo -e "SUBSETTED TO $( cat cleaned.subset.fam | wc -l | awk '{ print $1-1 }' ) RANDOM FAMILIES"
      else
        echo -e "NUMBER OF TRIOS DETECTED ( $n_trios ) LESS THAN MAX_TRIOS ( ~{max_trios} ); PROCEEDING WITHOUT DOWNSAMPLING"
        cp ~{prefix}.cleaned.fam cleaned.subset.fam
      fi
      
      # Run family analysis
      echo -e "STARTING FAMILY-BASED ANALYSIS"
      /opt/sv-pipeline/scripts/vcf_qc/analyze_fams.R \
        -S /opt/sv-pipeline/scripts/vcf_qc/SV_colors.txt \
        ~{vcf_stats} \
        cleaned.subset.fam \
        ~{prefix}_perSample/ \
        ~{prefix}_perFamily_plots/ \
        denovo_dump.rdx \
        --maxgq ~{max_gq}

    else

      mkdir ~{prefix}_perFamily_plots/

    fi

    # Prepare output
    echo -e "COMPRESSING RESULTS AS A TARBALL"
    tar -czvf ~{prefix}.plotQC_perFamily.tar.gz \
      ~{prefix}_perFamily_plots
  >>>

  output {
    File perFamily_plots_tarball = "~{prefix}.plotQC_perFamily.tar.gz"
    File cleaned_fam_file = "~{prefix}.cleaned.fam"
    File denovo_dump = "denovo_dump.rdx"
  }
}
