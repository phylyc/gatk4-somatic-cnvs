version development


struct Runtime {
    String gatk_docker
    File? gatk_override
    Int max_retries
    Int preemptible
    Int cpu
    Int machine_mem
    Int command_mem
    Int disk
    Int boot_disk_size
}

workflow callCollectReadCounts {
    input {
        File interval_list

        File ref_fasta
        File ref_fasta_index
        File ref_dict

        File bam
        File bai

        String format = "TSV"

        File? annotated_interval_list
        File? read_count_panel_of_normals

        # runtime
        String gatk_docker = "broadinstitute/gatk"
        File? gatk_override
        Int preemptible = 1
        Int max_retries = 2
        Int emergency_extra_diskGB = 0

        # memory assignments in MB
        Int get_sample_name_mem = 512  # 256
    }

    Int gatk_override_size = if defined(gatk_override) then ceil(size(gatk_override, "GB")) else 0
    Int disk_padGB = 1 + gatk_override_size + emergency_extra_diskGB

    Runtime standard_runtime = {
        "gatk_docker": gatk_docker,
        "gatk_override": gatk_override,
        "max_retries": max_retries,
        "preemptible": preemptible,
        "cpu": 1,
        "machine_mem": 2024,
        "command_mem": 2024,
        "disk": 1 + disk_padGB,
        "boot_disk_size": 12  # needs to be > 10
    }

    call GetSampleName {
        input:
            bam = bam,
            runtime_params = standard_runtime,
            memoryMB = get_sample_name_mem
    }

	call CollectReadCounts {
		input:
			interval_list = interval_list,
            ref_fasta = ref_fasta,
			ref_fasta_index = ref_fasta_index,
			ref_dict = ref_dict,
            bam = bam,
            bai = bai,
            format = format,
            sample_name = GetSampleName.sample_name,
            runtime_params = standard_runtime,
	}

    if (defined(annotated_interval_list) || defined(read_count_panel_of_normals)) {
        call DenoiseReadCounts {
            input:
                read_counts = CollectReadCounts.read_counts,
                sample_name = GetSampleName.sample_name,
                annotated_interval_list = annotated_interval_list,
                count_panel_of_normals = read_count_panel_of_normals,
                runtime_params = standard_runtime,
        }
    }

    output {
        File read_counts = CollectReadCounts.read_counts
        File? denoised_read_counts = DenoiseReadCounts.denoised_read_counts
        File? standardized_copy_ratios = DenoiseReadCounts.standardized_copy_ratios
    }
}

task GetSampleName {
    input {
        File bam

        Runtime runtime_params
        Int? memoryMB = 256
    }

    parameter_meta {
        bam: {localization_optional: true}
    }

    command <<<
        set -e
        export GATK_LOCAL_JAR=~{default="/root/gatk.jar" runtime_params.gatk_override}
        gatk --java-options "-Xmx~{select_first([memoryMB, runtime_params.command_mem])}m" \
            GetSampleName \
            -I ~{bam} \
            -O bam_name.txt \
            -encode
        cat bam_name.txt
    >>>

    output {
        String sample_name = read_string(stdout())
    }

    runtime {
        docker: runtime_params.gatk_docker
        bootDiskSizeGb: runtime_params.boot_disk_size
        memory: select_first([memoryMB, runtime_params.machine_mem]) + " MB"
        disks: "local-disk " + runtime_params.disk + " HDD"
        preemptible: runtime_params.preemptible
        maxRetries: runtime_params.max_retries
        cpu: runtime_params.cpu
    }
}

task CollectReadCounts {
    input {
        File interval_list
        File ref_dict
        File ref_fasta
        File ref_fasta_index
        File bam
        File bai
        String sample_name
        String interval_merging_rule = "OVERLAPPING_ONLY"
        String format = "TSV"

        Runtime runtime_params
        Int? memoryMB = 4096
    }

    parameter_meta {
        interval_list: {localization_optional: true}
        ref_fasta: {localization_optional: true}
        ref_fasta_index: {localization_optional: true}
        ref_dict: {localization_optional: true}
        bam: {localization_optional: true}
        bai: {localization_optional: true}
    }

    String output_name = sample_name + "_counts.tsv"

	command <<<
        set -e
        export GATK_LOCAL_JAR=~{default="/root/gatk.jar" runtime_params.gatk_override}
        gatk --java-options "-Xmx~{select_first([memoryMB, runtime_params.command_mem])}m" \
            CollectReadCounts \
            -I ~{bam} \
            -L ~{interval_list} \
            -R ~{ref_fasta} \
            -O ~{output_name} \
            --interval-merging-rule ~{interval_merging_rule} \
            --format ~{format}
	>>>

	output {
		File read_counts = output_name
	}

    runtime {
        docker: runtime_params.gatk_docker
        bootDiskSizeGb: runtime_params.boot_disk_size
        memory: select_first([memoryMB, runtime_params.machine_mem]) + " MB"
        disks: "local-disk " + runtime_params.disk + " HDD"
        preemptible: runtime_params.preemptible
        maxRetries: runtime_params.max_retries
        cpu: runtime_params.cpu
    }
}

task DenoiseReadCounts {
    input {
        File read_counts
        String sample_name
        File? annotated_interval_list
        File? count_panel_of_normals

        Runtime runtime_params
        Int? memoryMB = 4096
    }

    parameter_meta {
        read_counts: {localization_optional: true}
        annotated_interval_list: {localization_optional: true}
        count_panel_of_normals: {localization_optional: true}
    }

    String base_name = basename(read_counts, ".tsv")
    String output_denoised_read_counts = base_name + "_denoised.tsv"
    String output_standardized_copy_ratios = sample_name + "_standardized_copy_ratios.tsv"

	command <<<
        set -e
        export GATK_LOCAL_JAR=~{default="/root/gatk.jar" runtime_params.gatk_override}
        gatk --java-options "-Xmx~{select_first([memoryMB, runtime_params.command_mem])}m" \
            DenoiseReadCounts \
            -I ~{read_counts} \
            --denoised-copy-ratios ~{output_denoised_read_counts} \
            --standardized-copy-ratios ~{output_standardized_copy_ratios} \
            ~{"--annotated-intervals " + annotated_interval_list} \
            ~{"--count-panel-of-normals " + count_panel_of_normals}
	>>>

	output {
        File denoised_read_counts = output_denoised_read_counts
        File standardized_copy_ratios = output_standardized_copy_ratios
	}

    runtime {
        docker: runtime_params.gatk_docker
        bootDiskSizeGb: runtime_params.boot_disk_size
        memory: select_first([memoryMB, runtime_params.machine_mem]) + " MB"
        disks: "local-disk " + runtime_params.disk + " HDD"
        preemptible: runtime_params.preemptible
        maxRetries: runtime_params.max_retries
        cpu: runtime_params.cpu
    }
}