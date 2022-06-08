version 1.0

workflow callSamtoolsCoverage {
    input {
        File bam
        File bai
    }

    # Downstream analysis requires sample name to be name of the bam.
    # Since the bam name is not unique, we use the folder structure too:
    # gs://bucket/project/protocol/data_type/sample_id/version/sample_id.bam
#    String bam_name = basename(bam, ".bam")
    String subbed_bam_path = sub(bam, "/", "+")

	call Coverage {
		input:
            bam=bam,
            bai=bai,
            sample_name=subbed_bam_path,
	}

    output {
        File coverage = Coverage.coverage
    }
}

task Coverage {
    input {
        File bam
        File bai
        String sample_name

        # runtime
        String docker = "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.3.3-1513176735"
        Int extra_disk_spaceGB = 10
        Int memoryGB = 4
        Int preemptible = 1
        Int max_retries = 1
    }

    Int disk_spaceGB = ceil(size(bam, "GB") + size(bai, "GB")) + extra_disk_spaceGB

    # samtools is not capable of streaming in files, so has to localize
#    parameter_meta {
#        bam: {localization_optional: true}
#        bai: {localization_optional: true}
#    }

    String output_name = sample_name + ".txt"

	command <<<
        set -e
        samtools coverage \
            ~{bam} \
            --output ~{output_name}
	>>>

	output {
		File coverage = output_name
	}

    runtime {
        docker: docker
        bootDiskSizeGb: 12
        memory: memoryGB + " GB"
        disks: "local-disk " + disk_spaceGB + " HDD"
        preemptible: preemptible
        maxRetries: max_retries
    }
}