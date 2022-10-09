version development

import "https://github.com/phylyc/gatk4-somatic-cnvs/raw/main/gistic2.wdl" as gistic2


workflow Gistic2_per_Sample {
    input {
        File seg_file
        File samples_table
        File? markers_file
        File refgene_file
        File cnv_files

        Float amp_thresh = 0.3  # default 0.1
        Float del_thresh = 0.1
        Float qv_thresh = 0.1  # default 0.25
        Float cap = 1.5  # default 1.5
        Float broad_length_cutoff = 0.5  # default 0.98
        Float conf_level = 0.9  # default 0.75
        Int join_segment_size = 10
        Int max_sample_segs = 20000
        Int max_marker_spacing = 10000  # default 10000
        Int arm_peel = 1  # default 0
        Int do_gene_gistic = 1  # default 0
        Int remove_X = 1  # default 1
        String gene_collapse_method = "extreme"  # default "mean"

        String docker = "broadinstitute/gatk"

        Int memoryMB = 2048
        Int disk_size = 12
        Int preemptible = 1
    }

    call get_samples_to_scatter {
        input:
            samples_table = samples_table,
            docker = docker
    }

    scatter (sample in get_samples_to_scatter.samples) {
        call aggregate_segs_by_patient {
            input:
                seg_file = seg_file,
                samples_table = samples_table,
                sample = sample,
                docker = docker
        }

        call gistic2.Gistic2 {
            input:
                seg_file = aggregate_segs_by_patient.aggregated_seg_file,
                markers_file = markers_file,
                refgene_file = refgene_file,
                cnv_files = cnv_files,

                amp_thresh = amp_thresh,
                del_thresh = del_thresh,
                qv_thresh = qv_thresh,
                cap = cap,
                broad_length_cutoff = broad_length_cutoff,
                conf_level = conf_level,
                join_segment_size = join_segment_size,
                max_sample_segs = max_sample_segs,
                max_marker_spacing = max_marker_spacing,
                arm_peel = arm_peel,
                do_gene_gistic = do_gene_gistic,
                remove_X = remove_X,
                gene_collapse_method = gene_collapse_method,

                memoryMB = memoryMB,
                disk_size = disk_size,
                preemptible = preemptible
        }
    }

    call merge_gistic_output {
        input:
            samples_table = samples_table,
            samples_to_scatter = get_samples_to_scatter.samples,
            all_data_by_genes_array = Gistic2.all_data_by_genes,
            all_thresholded_by_genes_array = Gistic2.all_thresholded_by_genes,
            broad_values_by_arm_array = Gistic2.broad_values_by_arm,
            docker = docker
    }

    output {
        File all_data_by_genes = merge_gistic_output.all_data_by_genes
        File all_thresholded_by_genes = merge_gistic_output.all_thresholded_by_genes
        File broad_values_by_arm = merge_gistic_output.broad_values_by_arm
    }
}

task get_samples_to_scatter {
    input {
        File samples_table
        String docker
    }

    File samples_to_scatter = "samples_to_scatter.txt"

    command <<<
        touch ~{samples_to_scatter}

        python <<CODE
            import pandas as pd

            samples_table = pd.read_csv(~{samples_table}, sep="\t")
            samples_to_scatter = []
            for patient, sample_group in samples_table.groupby(by="Patient"):
                if sample_group["Sample"].nunique() > 1:
                    samples_to_scatter += list(sample_group["Sample"].to_numpy())

            # Ensure that the list is not empty by taking the first sample if it is:
            if not len(samples_to_scatter):
                samples_to_scatter = [samples_table.loc[0, "Sample"]]

            with open("~{samples_to_scatter}", "w") as f:
                for sample in samples_to_scatter:
                    f.write(f"{sample}\n")
        CODE
    >>>

    output {
        Array[String] samples = read_lines(samples_to_scatter)
    }

    runtime {
        docker: docker
        bootDiskSizeGb: 12
        memory: 512 + " MB"
        disks: "local-disk " + 1 + " HDD"
        preemptible: 1
        maxRetries: 1
        cpu: 1
    }
}

task aggregate_segs_by_patient {
    input {
        File seg_file
        File samples_table
        String sample

        String docker
    }

    File out_file = "aggregated_seg_file_" + sample + ".tsv"

    command <<<
        echo "Choosing ~{sample} as representative sample for its patient."
        echo "Aggregating copy ratios per segment over all samples per patient."

        python <<CODE
            import numpy as np
            import pandas as pd

            samples_table = pd.read_csv(~{samples_table}, sep="\t")
            segs = pd.read_csv(~{seg_file}, sep="\t")
            segs = segs.merge(samples_table, on="Sample")
            segs["tCR"] = 2 ** (segs["Seg.CN"] + 1)

            patient = samples_table.set_index("Sample").loc[~{sample}, "Patient"]
            segs_to_aggregate = segs.loc[segs["Patient"] != patient]
            representative_segs = segs.loc[segs["Sample"] == ~{sample}]

            agg_segs = segs_to_aggregate.groupby(
                ["Chromosome", "Start Position", "End Position", "Patient"]
            ).agg(
                **{
                    "tCR": pd.NamedAgg(column="tCR", aggfunc="mean"),
                    # Num Markers were already adjusted during pre-processing.
                    # If your data did not account for that, you need to do more here:
                    "Num Markers": pd.NamedAgg(column="Num Markers", aggfunc="max")
                }
            ).reset_index()
            agg_segs["Seg.CN"] = np.log2(agg_segs["tCR"]) - 1
            agg_segs.drop_duplicates(inplace=True)
            agg_segs.rename(columns={"Patient": "Sample"}, inplace=True)

            gistic_input_columns = ["Sample", "Chromosome", "Start Position", "End Position", "Num Markers", "Seg.CN"]
            agg_segs = pd.concat([agg_segs[gistic_input_columns], representative_segs[gistic_input_columns]])
            agg_segs.sort_values(by=["Sample", "Chromosome", "Start Position"], inplace=True)
            agg_segs.to_csv(~{out_file}, sep="\t", index=False)
        CODE
    >>>

    output {
        File aggregated_seg_file = out_file
    }

    runtime {
        docker: docker
        bootDiskSizeGb: 12
        memory: 1024 + " MB"
        disks: "local-disk " + 1 + " HDD"
        preemptible: 1
        maxRetries: 1
        cpu: 1
    }
}

task merge_gistic_output {
    input {
        File samples_table
        Array[String] samples_to_scatter
        Array[File] all_data_by_genes_array
        Array[File] all_thresholded_by_genes_array
        Array[File] broad_values_by_arm_array

        String docker
    }

    File all_data_by_genes_merged = "all_data_by_genes_merged.txt"
    File all_thresholded_by_genes_merged = "all_thresholded_by_genes_merged.txt"
    File broad_values_by_arm_merged = "broad_values_by_arm_merged.txt"

    command <<<
        python <<CODE
            import pandas as pd

            samples_table = pd.read_csv(~{samples_table}, sep="\t")
            scatter_mask = samples_table["Sample"].isin(samples_to_scatter)
            patient_to_sample = samples_table.loc[~scatter_mask].set_index("Patient")["Sample"].to_dict()
            single_patients = list(patient_to_sample.keys())

            def merge_tables(file_array, index_col, out_file):
                to_concat = []
                to_avg = []
                for file in file_array:
                    table = pd.read_csv(file, sep="\t", index_col=index_col)
                    to_concat.append(table[[c for c in table.columns if c in samples_to_scatter]])
                    to_avg.append(table[[c for c in table.columns if c in single_patients]])
                avg = pd.concat(to_avg).groupby(level=0).mean()
                to_concat.append(avg.rename(columns=patient_to_sample))
                merged_table = pd.concat(to_concat, axis=1)
                merged_table.to_csv(out_file, sep="\t")

            merge_tables(~{all_data_by_genes_array}, ["Gene Symbol", "Gene ID", "Cytoband"], ~{all_data_by_genes_merged})
            merge_tables(~{all_thresholded_by_genes_array}, ["Gene Symbol", "Locus ID", "Cytoband"], ~{all_thresholded_by_genes_merged})
            merge_tables(~{broad_values_by_arm_array}, "Chromosome Arm", ~{broad_values_by_arm_merged})
        CODE
    >>>

    output {
        File all_data_by_genes = all_data_by_genes_merged
        File all_thresholded_by_genes = all_thresholded_by_genes_merged
        File broad_values_by_arm = broad_values_by_arm_merged
    }

    runtime {
        docker: docker
        bootDiskSizeGb: 12
        memory: 1024 + " MB"
        disks: "local-disk " + 1 + " HDD"
        preemptible: 1
        maxRetries: 1
        cpu: 1
    }
}