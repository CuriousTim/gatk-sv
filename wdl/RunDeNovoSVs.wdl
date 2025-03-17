version 1.0

###########################
# IMPORT TOOLS
###########################

import "Structs.wdl"
import "ReformatRawFiles.wdl" as raw
import "TasksMakeCohortVcf.wdl" as miniTasks
import "DeNovoSVsScatter.wdl" as runDeNovo
import "Utils.wdl" as util

###########################
# MAIN WORKFLOW DEFINITION
###########################

workflow DeNovoSV {

    # Define the input values for the workflow
    input {

        # Core input files
        File ped_input
        File? family_ids_txt               # Kept to run tests on family subsets
        File vcf_file
        File vcf_file_index                    
        File? genomic_disorder_input       # Kept in case we want to use the file for filtering within the python script.
        File exclude_regions

        # Running parameters
        String prefix
        Int records_per_shard
        Array[String] contigs = ["chr1", "chr2", "chr3", "chr4", "chr5", 
                                "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14",
                                "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX"]

        # Fetch raw data
        Array[String] batch_name_list
        Array[File] batch_sample_lists
        Array[String] batch_bincov_matrix
        Array[String] batch_bincov_matrix_index
        Array[String]? clustered_manta_vcf
        Array[String]? clustered_melt_vcf
        Array[String]? clustered_wham_vcf
        Array[String]? clustered_scramble_vcf
        Array[String] clustered_depth_vcf

        # Parameters for denovo_svs.py with default values
        # Size parameters
        Int small_cnv_size = 1000                # Renamed
        Int intermediate_cnv_size = 5000         # New
        Int depth_only_size = 10000 
        Int exclude_parent_cnv_size = 10000000   # New
        # Allele frequency
        Float gnomad_af = 0.01
        Float parents_af = 0.05
        Float cohort_af = 0.05
        # Overlap parameters
        Float large_raw_overlap = 0.5
        Float small_raw_overlap = 0.5
        Float parents_overlap = 0.5
        Float blacklist_overlap = 0.5            # New
        Int nearby_insertion = 100               # New
        # SV quality (parents)
        Int coverage_cutoff = 10
        Float gq_min = 0
        # Other
        String gnomad_col = gnomAD_V2_AF
        String alt_gnomad_col = gnomad_v2.1_sv_AF     # Check if that one is necessary
        String? af_column_name                        # Not sure what this one is? 

        # Parameters for denovo_outliers.py with default values
        #Int denovo_outlier_factor = 3     # New! To implement

        # Dockers
        String variant_interpretation_docker
        String sv_pipeline_updates_docker
        String linux_docker
        String python_docker                                # Kept this one in in case we'd like to make use of the possibility to run on a subset of the cohort
        String sv_base_mini_docker                          # I saw that you removed that one and also removed it in the code; Ok to do that on this version as well? 

        # Runtime attributes (ordered by task calling)
        RuntimeAttr? runtime_attr_make_manifests
        RuntimeAttr? runtime_attr_get_batched_files
        RuntimeAttr? runtime_attr_subset_by_samples
        RuntimeAttr? runtime_attr_clean_ped
        RuntimeAttr? runtime_attr_vcf_to_bed
        RuntimeAttr? runtime_attr_raw_vcf_to_bed
        RuntimeAttr? runtime_attr_raw_merge_bed
        RuntimeAttr? runtime_attr_raw_divide_by_chrom
        RuntimeAttr? runtime_attr_raw_reformat_bed
        RuntimeAttr? runtime_attr_subset_vcf
        RuntimeAttr? runtime_attr_shard_vcf
        RuntimeAttr? runtime_attr_denovo
        RuntimeAttr? runtime_attr_merge_final_bed_files
        RuntimeAttr? runtime_attr_call_outliers
        RuntimeAttr? runtime_attr_create_plots
        # Removed attributes: not used anywhere: runtime_attr_merge, runtime_attr_batch_vcf (this one was kept in the denovo-dev branch from variant-interpretation but I couldn't find it being used), runtime_attr_gd & runtime_attr_merge_gd (tasks were deleted)

    }

    # Create text files with the paths of the raw data used by the workflow
    call makeManifests {
        input:
            batch_name_list = batch_name_list,
            batch_sample_lists = batch_sample_lists,
            batch_bincov_matrix = batch_bincov_matrix,
            batch_bincov_matrix_index = batch_bincov_matrix_index,
            clustered_manta_vcf = clustered_manta_vcf,
            clustered_melt_vcf = clustered_melt_vcf,
            clustered_wham_vcf = clustered_wham_vcf,
            clustered_scramble_vcf = clustered_scramble_vcf,
            clustered_depth_vcf = clustered_depth_vcf,
            runtime_docker = linux_docker,
            runtime_attr_override = runtime_attr_make_manifests
    }

    # If family_ids_txt is provided (text file with one family_id/line for all families to be analyzed), subset all other input files to only include the necessary batches.
    if (defined(family_ids_txt)) {
        File family_ids_txt_ = select_first([family_ids_txt])
        call GetBatchedFiles {
            input:
                batch_raw_file = batch_raw_file,
                batch_depth_raw_file = batch_depth_raw_file,
                ped_input = ped_input,
                family_ids_txt = family_ids_txt_,
                sample_batches = sample_batches,
                batch_bincov_index = batch_bincov_index,
                python_docker=python_docker,
                runtime_attr_override = runtime_attr_get_batched_files
        }

        call util.SubsetVcfBySamplesList {
            input:
                vcf = vcf_file,
                vcf_idx = vcf_index,
                list_of_samples = GetBatchedFiles.samples,
                outfile_name = prefix,
                keep_af = true,
                remove_private_sites = false,
                sv_base_mini_docker = sv_base_mini_docker,
                runtime_attr_override = runtime_attr_subset_by_samples
        }
    }

    # Makes a ped file of singletons, duos, and trios for input into the de novo script (only including families of interest)
    call CleanPed {
        input:
            ped_input = ped_input,
            vcf_input = select_first([SubsetVcfBySamplesList.vcf_subset, vcf_file]),
            variant_interpretation_docker=variant_interpretation_docker,
            runtime_attr_override = runtime_attr_clean_ped
    }

    # Splits raw files into probands and parents and reformats to have chrom_svtype_sample as the first column for probands and chrom_svtype_famid as the first column for parents
    call raw.ReformatRawFiles as ReformatRawFiles {
        input:
            contigs = contigs,
            raw_files_list = select_first([GetBatchedFiles.batch_raw_files_list, batch_raw_file]),
            ped_input = CleanPed.cleaned_ped,
            depth = false,
            variant_interpretation_docker = variant_interpretation_docker,
            sv_base_mini_docker = sv_base_mini_docker,
            runtime_attr_vcf_to_bed = runtime_attr_raw_vcf_to_bed,
            runtime_attr_merge_bed = runtime_attr_raw_merge_bed,
            runtime_attr_divide_by_chrom = runtime_attr_raw_divide_by_chrom,
            runtime_attr_reformat_bed = runtime_attr_raw_reformat_bed
    }

    # Splits raw files into probands and parents and reformats to have chrom_svtype_sample as the first column for probands and chrom_svtype_famid as the first column for parents
    call raw.ReformatRawFiles as ReformatDepthRawFiles {
        input:
            contigs = contigs,
            raw_files_list = select_first([GetBatchedFiles.batch_depth_raw_files_list, batch_depth_raw_file]),
            ped_input = CleanPed.cleaned_ped,
            depth = true,
            variant_interpretation_docker = variant_interpretation_docker,
            sv_base_mini_docker = sv_base_mini_docker,
            runtime_attr_vcf_to_bed = runtime_attr_raw_vcf_to_bed,
            runtime_attr_merge_bed = runtime_attr_raw_merge_bed,
            runtime_attr_divide_by_chrom = runtime_attr_raw_divide_by_chrom,
            runtime_attr_reformat_bed = runtime_attr_raw_reformat_bed
    }

    # Scatter following tasks across chromosomes: SubsetVcf; ScatterVcf (miniTasks.ScatterVcf); GetDeNovo (runDeNovo.DeNovoSVsScatter)
    scatter (i in range(length(contigs))) {

        # Splits vcf by chromosome
        call SubsetVcf {
            input:
                vcf_file = select_first([SubsetVcfBySamplesList.vcf_subset, vcf_file]),
                chromosome=contigs[i],
                variant_interpretation_docker=variant_interpretation_docker,
                runtime_attr_override = runtime_attr_subset_vcf
        }

        # Shards vcf
        call miniTasks.ScatterVcf as ScatterVcf {
            input:
                vcf=SubsetVcf.vcf_output,
                prefix=prefix,
                records_per_shard=records_per_shard,
                sv_pipeline_docker=sv_pipeline_updates_docker,
                runtime_attr_override=runtime_attr_shard_vcf
        }
    
        # Runs the de novo calling python script on each shard and outputs a per chromosome list of de novo SVs
        call runDeNovo.DeNovoSVsScatter as GetDeNovo {
            input:
                ped_input=CleanPed.cleaned_ped,
                vcf_files=ScatterVcf.shards,
                chromosome=contigs[i],
                raw_proband=ReformatRawFiles.reformatted_proband_raw_files[i],
                raw_parents=ReformatRawFiles.reformatted_parents_raw_files[i],
                raw_depth_proband=ReformatDepthRawFiles.reformatted_proband_raw_files[i],
                raw_depth_parents=ReformatDepthRawFiles.reformatted_parents_raw_files[i],
                exclude_regions = exclude_regions,
                sample_batches = sample_batches,
                batch_bincov_index = select_first([GetBatchedFiles.batch_bincov_index_subset, batch_bincov_index]),
                python_config=python_config,
                variant_interpretation_docker=variant_interpretation_docker,
                runtime_attr_denovo = runtime_attr_denovo,
                runtime_attr_vcf_to_bed = runtime_attr_vcf_to_bed
        }
    }

    # Merges the per chromosome final de novo SV outputs
    call MergeDenovoBedFiles {
        input:
            bed_files = GetDeNovo.per_chromosome_final_output_file,
            variant_interpretation_docker=variant_interpretation_docker,
            runtime_attr_override = runtime_attr_merge_final_bed_files
    }

    # Outputs a final callset of de novo SVs as well as outlier de novo SV calls
    call CallOutliers {
        input:
            bed_file = MergeDenovoBedFiles.merged_denovo_output,
            variant_interpretation_docker=variant_interpretation_docker,
            runtime_attr_override = runtime_attr_call_outliers
    }

    # Generates plots for QC
    call CreatePlots {
        input:
            bed_file = CallOutliers.final_denovo_nonOutliers_output,
            ped_input = ped_input,
            variant_interpretation_docker=variant_interpretation_docker,
            runtime_attr_override = runtime_attr_create_plots
    }

    output {
        File cleaned_ped = CleanPed.cleaned_ped
        File final_denovo_nonOutliers = CallOutliers.final_denovo_nonOutliers_output
        File final_denovo_outliers = CallOutliers.final_denovo_outliers_output
        File final_denovo_nonOutliers_plots = CreatePlots.output_plots
        Array [File] denovo_output_annotated = GetDeNovo.per_chromosome_annotation_output_file
    }
}


###########################
# TASK DEFINITION
###########################

############ makeManifests ############
task makeManifests {
    input {
        Array[String] batch_name_list
        Array[File] batch_sample_lists
        Array[String]? clustered_manta_vcf
        Array[String]? clustered_melt_vcf
        Array[String]? clustered_wham_vcf
        Array[String]? clustered_scramble_vcf
        Array[String] clustered_depth_vcf
        Array[String] batch_bincov_matrix
        Array[String] batch_bincov_matrix_index
        String runtime_docker

        RuntimeAttr? runtime_attr_override
    }

    Float input_size = size(batch_sample_lists, "GB")
    RuntimeAttr runtime_default = object {
                                      mem_gb: 1,
                                      disk_gb: 16 + ceil(input_size),
                                      cpu: 1,
                                      preemptible: 3,
                                      max_retries: 1,
                                      boot_disk_gb: 16
                                    }
    RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])
    Array[String] manta_vcfs = select_first([clustered_manta_vcf, []])
    Array[String] melt_vcfs = select_first([clustered_melt_vcf, []])
    Array[String] wham_vcfs = select_first([clustered_wham_vcf, []])
    Array[String] scramble_vcfs = select_first([clustered_scramble_vcf, []])

    runtime {
        cpu: select_first([runtime_override.cpu, runtime_default.cpu])
        memory: "${select_first([runtime_override.mem_gb, runtime_default.mem_gb])} GB"
        disks: "local-disk ${select_first([runtime_override.disk_gb, runtime_default.disk_gb])} HDD"
        preemptible: select_first([runtime_override.preemptible, runtime_default.preemptible])
        maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
        docker: runtime_docker
        bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
    }

    command <<<
        set -euo pipefail

        batch_names='~{write_lines(batch_name_list)}'

        paste "${batch_names}" '~{write_lines(clustered_depth_vcf)}' > 'depth_manifest.tsv'

        manta='~{if length(manta_vcfs) > 0 then write_lines(manta_vcfs) else ""}'
        melt='~{if length(melt_vcfs) > 0 then write_lines(melt_vcfs) else ""}'
        wham='~{if length(wham_vcfs) > 0 then write_lines(wham_vcfs) else ""}'
        scramble='~{if length(scramble_vcfs) > 0 then write_lines(scramble_vcfs) else ""}'

        pesr_manifest='pesr_manifest.tsv'
        : > "${pesr_manifest}"

        if [[ "${manta}" ]]; then
            paste "${batch_names}" "${manta}" >> "${pesr_manifest}"
        fi
        if [[ "${melt}" ]]; then
            paste "${batch_names}" "${melt}" >> "${pesr_manifest}"
        fi
        if [[ "${wham}" ]]; then
            paste "${batch_names}" "${wham}" >> "${pesr_manifest}"
        fi
        if [[ "${scramble}" ]]; then
            paste "${batch_names}" "${scramble}" >> "${pesr_manifest}"
        fi

        if [[ ! -s "${pesr_manifest}" ]]; then
            printf 'at least one non-empty list of PESR evidence VCF lists should be provided\n' >&2
            exit 1
        fi

        paste "${batch_names}" '~{write_lines(batch_sample_lists)}' \
          | awk -F'\t' '{while((getline line < $2) > 0) {print line "\t" $1}}' > 'sample_manifest.tsv'
        
        paste "${batch_names}" '~{write_lines(batch_bincov_matrix)}' '~{write_lines(batch_bincov_matrix_index)}' > 'bincov_manifest.tsv'
    >>>

    output {
        File depth_manifest = "depth_manifest.tsv"
        File pesr_manifest = "pesr_manifest.tsv"
        File bincov_manifest = "bincov_manifest.tsv"
        File sample_manifest = "sample_manifest.tsv"
    }
}


############ GetBatchedFiles ############
task GetBatchedFiles {

    input {
        File batch_raw_file
        File batch_depth_raw_file
        File family_ids_txt
        File ped_input
        File sample_batches
        File batch_bincov_index
        String python_docker
        RuntimeAttr? runtime_attr_override
    }

    Float input_size = size(select_all([batch_raw_file, batch_depth_raw_file, ped_input, sample_batches, batch_bincov_index, family_ids_txt]), "GB")

    RuntimeAttr default_attr = object {
        mem_gb: 3.75,
        disk_gb: ceil(10 + input_size),
        cpu_cores: 1,
        preemptible_tries: 2,
        max_retries: 1,
        boot_disk_gb: 8
    }
    
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    output {
        File batch_raw_files_list = "batch_raw_files_list.txt"
        File batch_depth_raw_files_list = "batch_depth_raw_files_list.txt"
        File batch_bincov_index_subset = "batch_bincov_index.txt"
        File samples = "samples.txt"
    }

    command {
        set -exuo pipefail

        if grep -q -w -f ~{family_ids_txt} ~{ped_input}; then
            grep -w -f ~{family_ids_txt} ~{ped_input} | cut -f2 | sort -u > samples.txt
        else
            echo "No matching family IDs from family_ids_txt found in ped_input file." >&2
            exit 1
        fi

        if grep -q -w -f samples.txt ~{sample_batches}; then
            grep -w -f samples.txt ~{sample_batches} | cut -f2 | sort -u > batches.txt
        else
            echo "No matching individual IDs found in the sample_batches file." >&2
            exit 1
        fi

        grep -w -f batches.txt ~{batch_bincov_index} > batch_bincov_index.txt
        grep -w -f batches.txt ~{batch_raw_file} > batch_raw_files_list.txt
        grep -w -f batches.txt ~{batch_depth_raw_file} > batch_depth_raw_files_list.txt
    }

    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
        docker: python_docker
    }
}


############ CleanPed ############
task CleanPed {

    input {
        File ped_input
        File vcf_input
        String variant_interpretation_docker
        RuntimeAttr? runtime_attr_override
    }

    Float ped_size = size(ped_input, "GB")
    Float vcf_size = size(vcf_input, "GB")

    RuntimeAttr default_attr = object {
        mem_gb: 3.75,
        disk_gb: ceil(10 + vcf_size + ped_size * 1.5),
        cpu_cores: 1,
        preemptible_tries: 2,
        max_retries: 1,
        boot_disk_gb: 8
    }
    
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    output {
        File cleaned_ped = "subset_cleaned_ped.txt"
    }

    command {
        set -exuo pipefail

        # TODO: this script should get the name of the output file it generates as an input argument.
        # The output filename is currently hardcoded to be 'clean_ped.txt'.
        Rscript /src/denovo/clean_ped.R ~{ped_input}
        cut -f2 cleaned_ped.txt | awk 'NR > 1' > all_samples.txt
        bcftools query -l ~{vcf_input} > samples_to_include_in_ped.txt

        # Note: grep returns a non-success exit code (i.e., other than `0`) when it cannot find the
        # match in the following scripts. We do not expect it to find a match for every entry.
        # Hence, to avoid exit with unsuccessful code, we can either drop pipefail from above or use `|| true`.

        grep -w -v -f samples_to_include_in_ped.txt all_samples.txt > excluded_samples.txt || true
        grep -w -f excluded_samples.txt cleaned_ped.txt | cut -f1 | sort -u > excluded_families.txt || true
        grep -w -v -f excluded_families.txt cleaned_ped.txt > subset_cleaned_ped.txt || true
    }

    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
        docker: variant_interpretation_docker
    }
}


############ SubsetVcf ############
task SubsetVcf {

    input {
        File vcf_file
        String chromosome
        String variant_interpretation_docker
        RuntimeAttr? runtime_attr_override
    }

    Float input_size = size(vcf_file, "GB")

    RuntimeAttr default_attr = object {
        mem_gb: 3.75,
        disk_gb: ceil(10 + input_size * 1.5),
        cpu_cores: 1,
        preemptible_tries: 2,
        max_retries: 1,
        boot_disk_gb: 8
    }
    
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    output {
        File vcf_output = "~{chromosome}.vcf.gz"
    }

    command <<<
        set -euxo pipefail

        bcftools index ~{vcf_file}
        bcftools view ~{vcf_file} --regions ~{chromosome} -O z -o  ~{chromosome}.vcf.gz
    >>>

    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
        docker: variant_interpretation_docker
    }
}


############ MergeDenovoBedFiles ############
task MergeDenovoBedFiles {

    input {
        Array[File] bed_files
        String variant_interpretation_docker
        RuntimeAttr? runtime_attr_override
    }

    Float bed_files_size = size(bed_files, "GB")

    RuntimeAttr default_attr = object {
        mem_gb: 3.75,
        disk_gb: ceil(10 + (bed_files_size) * 2.0),
        cpu_cores: 1,
        preemptible_tries: 2,
        max_retries: 1,
        boot_disk_gb: 8
    }
    
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    output {
        File merged_denovo_output = "denovo.merged.bed.gz"
    }

    command {
        set -exuo pipefail

        zcat ~{bed_files[0]} | awk 'NR<=1' > denovo.merged.bed
        zcat ~{sep=" " bed_files} | grep -v ^chrom >> denovo.merged.bed
        bgzip denovo.merged.bed
    }

    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
        docker: variant_interpretation_docker
    }
}


############ CallOutliers ############
task CallOutliers {

    input {
        File bed_file
        String variant_interpretation_docker
        RuntimeAttr? runtime_attr_override
    }

    Float bed_files_size = size(bed_file, "GB")

    RuntimeAttr default_attr = object {
        mem_gb: 16, # 3.75
        disk_gb: ceil(10 + (bed_files_size) * 2.0),
        cpu_cores: 1,
        preemptible_tries: 2,
        max_retries: 1,
        boot_disk_gb: 8
    }
    
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    output {
        File final_denovo_nonOutliers_output = "final.denovo.merged.bed.gz"
        File final_denovo_outliers_output = "final.denovo.merged.outliers.bed.gz"
        File final_denovo_allSamples_output = "final.denovo.merged.allSamples.bed.gz"
    }

    command {
        set -exuo pipefail

        python /src/denovo/denovo_outliers.py --bed ~{bed_file}

        bgzip final.denovo.merged.bed
        bgzip final.denovo.merged.outliers.bed
        bgzip final.denovo.merged.allSamples.bed
    }

    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
        docker: variant_interpretation_docker
    }
}


############ CreatePlots ############
task CreatePlots {

    input {
        File bed_file
        File ped_input
        String variant_interpretation_docker
        RuntimeAttr? runtime_attr_override
    }

    Float input_size = size(select_all([bed_file, ped_input]), "GB")

    RuntimeAttr default_attr = object {
        mem_gb: 16, # 3.75
        disk_gb: ceil(10 + input_size * 1.2),
        cpu_cores: 1,
        preemptible_tries: 2,
        max_retries: 1,
        boot_disk_gb: 8
    }
    
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    output {
        File output_plots = "output_plots.pdf"
    }

    command {
        set -exuo pipefail

        Rscript /src/denovo/denovo_sv_plots.R ~{bed_file} ~{ped_input} output_plots.pdf
    }

    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
        docker: variant_interpretation_docker
    }
}

