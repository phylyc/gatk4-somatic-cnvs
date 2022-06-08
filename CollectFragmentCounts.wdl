version 1.0

workflow callCollectFragmentCounts {
    input {
        File bam
        File bai
        String sample_name
        File ref_fasta
        File interval_list
        File ref_fasta_index
        File ref_dict
        String? gatk_docker
    }

	call CollectFragmentCounts {
		input:
			interval_list=interval_list,
            ref_fasta=ref_fasta,
			ref_fasta_index=ref_fasta_index,
			ref_dict=ref_dict,
            bam=bam,
            bai=bai,
            sample_name=sample_name,
            gatk_docker=gatk_docker
	}

    output {
        File fragment_counts = CollectFragmentCounts.fragment_counts
    }
}

task CollectFragmentCounts {
    input {
        File interval_list
        File ref_dict
        File ref_fasta
        File ref_fasta_index
        File bam
        File bai
        String sample_name
        String interval_merging_rule = "OVERLAPPING_ONLY"

        # runtime
        String gatk_docker = "broadinstitute/gatk:4.0.0.0"
        File? gatk_override
        Int disk_spaceGB = 2
        Int memoryGB = 8
        Int preemptible = 2
        Int max_retries = 2
        Int cpu = 1
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
        export GATK_LOCAL_JAR=~{default="/root/gatk.jar" gatk_override}
        gatk --java-options "-Xmx~{memoryGB}g" \
            CollectFragmentCounts \
            -I ~{bam} \
            -L ~{interval_list} \
            -R ~{ref_fasta} \
            -O ~{output_name} \
            -imr ~{interval_merging_rule} \
            --format TSV
	>>>

	output {
		File fragment_counts = output_name
	}

    runtime {
        docker: gatk_docker
        bootDiskSizeGb: 12
        memory: memoryGB + " GB"
        disks: "local-disk " + disk_spaceGB + " HDD"
        preemptible: preemptible
        maxRetries: max_retries
        cpu: cpu
    }
}