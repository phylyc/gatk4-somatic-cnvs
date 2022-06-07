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

workflow callAnnotateIntervals {
    input {
        File interval_list
        File ref_fasta
        File ref_fasta_index
        File ref_dict

        File? mappability_track
        File? segmental_duplication_track

        # runtime
        String gatk_docker = "broadinstitute/gatk"
        File? gatk_override
        Int preemptible = 1
        Int max_retries = 2
        Int emergency_extra_diskGB = 0

        # memory assignments in MB
        Int annotate_intervals_mem = 512
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

    call AnnotateIntervals {
        input:
            interval_list = interval_list,
            ref_fasta = ref_fasta,
            ref_fasta_index = ref_fasta_index,
            ref_dict = ref_dict,
            mappability_track = mappability_track,
            segmental_duplication_track = segmental_duplication_track,
            runtime_params = standard_runtime,
            memoryMB = annotate_intervals_mem
    }

    output {
        File annotated_interval_list = AnnotateIntervals.annotated_interval_list
    }
}

task AnnotateIntervals {
    input {
        File interval_list

        File ref_fasta
        File ref_fasta_index
        File ref_dict

        File? mappability_track
        File? segmental_duplication_track

        Runtime runtime_params
        Int? memoryMB = 512
    }

    String output_file = basename(interval_list, ".interval_list") + ".annotated.interval_list"

    parameter_meta {
        interval_list: {localization_optional: true}
        ref_fasta: {localization_optional: true}
        ref_fasta_index: {localization_optional: true}
        ref_dict: {localization_optional: true}
    }

	command <<<
        set -e
        export GATK_LOCAL_JAR=~{default="/root/gatk.jar" runtime_params.gatk_override}
        gatk --java-options "-Xmx~{select_first([memoryMB, runtime_params.command_mem])}m" \
            AnnotateIntervals \
            -R ~{ref_fasta} \
            -L ~{interval_list} \
            -O ~{output_file} \
            --interval-merging-rule OVERLAPPING_ONLY \
            ~{"--mappability-track " + mappability_track} \
            ~{"--segmental-duplication-track " + segmental_duplication_track}
	>>>

	output {
		File annotated_interval_list = output_file
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