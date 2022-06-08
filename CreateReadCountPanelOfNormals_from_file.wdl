version development

#import "CreateReadCountPanelOfNormals.wdl" as crcpon
import "https://github.com/phylyc/gatk4-somatic-cnvs/raw/main/CreateReadCountPanelOfNormals.wdl" as crcpon


workflow callCreateReadCountPanelOfNormals_from_File {
    input {
        File interval_list
        File ref_fasta
        File ref_fasta_index
        File ref_dict

        File normal_bams_file
        File normal_bais_file

        String pon_name

        File? annotated_interval_list
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

    call crcpon.callCreateReadCountPanelOfNormals {
        input:
            interval_list = interval_list,
            ref_fasta = ref_fasta,
            ref_fasta_index = ref_fasta_index,
            ref_dict = ref_dict,

            normal_bams = read_lines(normal_bams_file),
            normal_bais = read_lines(normal_bais_file),

            pon_name = pon_name,

            annotated_interval_list = annotated_interval_list,
            mappability_track = mappability_track,
            segmental_duplication_track = segmental_duplication_track,

            gatk_docker = gatk_docker,
            gatk_override = gatk_override,
            preemptible = preemptible,
            max_retries = max_retries,

            emergency_extra_diskGB = emergency_extra_diskGB,
            annotate_intervals_mem = annotate_intervals_mem,
            create_panel_mem = create_panel_mem
    }

    output {
        File pon = callCreateReadCountPanelOfNormals.cnv_pon
        File? annotated_interval_list = callCreateReadCountPanelOfNormals.annotated_interval_list
    }
}