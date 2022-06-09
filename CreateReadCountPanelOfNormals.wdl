version development

#import "CollectReadCounts.wdl" as crc
#import "AnnotateIntervals.wdl" as ai
import "https://github.com/phylyc/gatk4-somatic-cnvs/raw/main/CollectReadCounts.wdl" as crc
import "https://github.com/phylyc/gatk4-somatic-cnvs/raw/main/AnnotateIntervals.wdl" as ai


workflow callCreateReadCountPanelOfNormals {
    input {
        File interval_list
        File ref_fasta
        File ref_fasta_index
        File ref_dict

        Array[File]+ normal_bams
        Array[File]+ normal_bais

        String pon_name

        File? annotated_interval_list
        # If annotated_interval_list is specified, those args are ignored:
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
        Int create_panel_mem = 8192
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

    scatter (normal_bam in zip(normal_bams, normal_bais)) {
        call crc.callCollectReadCounts {
            input:
                interval_list = interval_list,
                ref_fasta = ref_fasta,
                ref_fasta_index = ref_fasta_index,
                ref_dict = ref_dict,
                bam = normal_bam.left,
                bai = normal_bam.right,
                format = "HDF5",
                gatk_override = gatk_override,
                gatk_docker = gatk_docker,
                preemptible = preemptible,
                max_retries = max_retries
        }
    }

    if (!defined(annotated_interval_list)) {
        call ai.AnnotateIntervals {
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
    }

    File this_annotated_interval_list = select_first([
        annotated_interval_list, AnnotateIntervals.annotated_interval_list
    ])

	call CreateReadCountPanelOfNormals {
		input:
            input_counts = callCollectReadCounts.read_counts,
            output_name = pon_name,
            annotated_interval_list = this_annotated_interval_list,
            runtime_params = standard_runtime,
            memoryMB = create_panel_mem
	}

    output {
        File pon = CreateReadCountPanelOfNormals.cnv_pon
        File? new_annotated_interval_list = AnnotateIntervals.annotated_interval_list
    }
}

task CreateReadCountPanelOfNormals {
    input {
        Array[File]+ input_counts
        String output_name

        File? annotated_interval_list
        Int number_of_eigensamples = 20

        Runtime runtime_params
        Int? memoryMB = 8192
    }

    String output_pon = output_name + ".hdf5"

    parameter_meta {
        input_counts: {localization_optional: true}
        annotated_interval_list: {localization_optional: true}
    }

	command <<<
        set -e
        export GATK_LOCAL_JAR=~{default="/root/gatk.jar" runtime_params.gatk_override}
        gatk --java-options "-Xmx~{select_first([memoryMB, runtime_params.command_mem])}m" \
            CreateReadCountPanelOfNormals \
            -I ~{sep=' -I ' input_counts} \
            -O ~{output_pon} \
            ~{"--annotated-intervals " + annotated_interval_list} \
            --number-of-eigensamples ~{number_of_eigensamples}
	>>>

	output {
		File cnv_pon = output_pon
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