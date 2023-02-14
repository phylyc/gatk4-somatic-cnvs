version 1.0

workflow callPileup {
    input {
        File bam
        File bai

        String sample_name = sub(bam, "/", "+")
    }

	call Pileup {
		input:
            bam=bam,
            bai=bai,
            sample_name=sample_name,
	}

    output {
        File pileup = Pileup.pileup
    }
}

task Pileup {
    input {
        File bam
        File bai
        String sample_name
        File? intervals
        File? ref_fasta
        File? ref_fasta_index

        Boolean output_insert_length = false

        # runtime
        String docker = "broadinstitute/gatk"
        File? gatk_override
        Int extra_disk_spaceGB = 10
        Int memoryMB = 2048
        Int preemptible = 1
        Int max_retries = 1
    }

    Int ref_size = if defined(ref_fasta) then ceil(size(ref_fasta, "GB") + size(ref_fasta_index, "GB")) else 0
    Int disk_spaceGB = extra_disk_spaceGB # + ceil(size(bam, "GB") + size(bai, "GB")) + ref_size

    parameter_meta {
        ref_fasta: {localization_optional: true}
        ref_fasta_index: {localization_optional: true}
        intervals: {localization_optional: true}
        bam: {localization_optional: true}
        bai: {localization_optional: true}
    }

    String output_name = sample_name + ".txt"

	command <<<
        set -e
        export GATK_LOCAL_JAR=~{default="/root/gatk.jar" gatk_override}
        gatk --java-options "-Xmx~{memoryMB}m" \
            Pileup \
            ~{"-R " + ref_fasta} \
            ~{"-L " + intervals} \
            -I ~{bam} \
            --output-insert-length ~{output_insert_length} \
            -O ~{output_name}
	>>>

	output {
		File pileup = output_name
	}

    runtime {
        docker: docker
        bootDiskSizeGb: 12
        memory: memoryMB + " MB"
        disks: "local-disk " + disk_spaceGB + " HDD"
        preemptible: preemptible
        maxRetries: max_retries
    }
}