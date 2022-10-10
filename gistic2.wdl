version development


workflow Gistic2 {
    input {
        File seg_file
        File refgene_file
        File? markers_file
        File? cnv_files

        Float amp_thresh = 0.3
        Float del_thresh = 0.1
        Float qv_thresh = 0.1
        Float cap = 1.5
        Float broad_length_cutoff = 0.7
        Float conf_level = 0.9
        Int join_segment_size = 10
        Int max_sample_segs = 20000
        Int max_marker_spacing = 10000
        Int arm_peel = 1
        Int do_gene_gistic = 1
        Int remove_X = 1
        String gene_collapse_method = "extreme"

        Int memoryMB = 2048
        Int disk_size = 1
        Int preemptible = 1
    }

    call tool_gistic2 {
        input:
            seg_file = seg_file,
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

    output {
        File all_data_by_genes = tool_gistic2.all_data_by_genes
        File all_lesions = tool_gistic2.all_lesions
        File all_thresholded_by_genes = tool_gistic2.all_thresholded_by_genes
        File amp_genes = tool_gistic2.amp_genes
        File amp_qplot_png = tool_gistic2.amp_qplot_png
        File arraylistfile = tool_gistic2.arraylistfile
        File broad_values_by_arm = tool_gistic2.broad_values_by_arm
        File broad_significance_results = tool_gistic2.broad_significance_results
        File del_genes = tool_gistic2.del_genes
        File del_qplot_png = tool_gistic2.del_qplot_png
        File gistic_inputs = tool_gistic2.gistic_inputs
        File gistic_version = tool_gistic2.gistic_version
        File segmented_copy_number_png = tool_gistic2.segmented_copy_number_png
        File gistic_scores = tool_gistic2.gistic_scores
    }
}

task tool_gistic2 {
    input {
        File seg_file
        File refgene_file
        File? markers_file
        File? cnv_files

        Float amp_thresh = 0.1
        Float del_thresh = 0.1
        Float qv_thresh = 0.25
        Float cap = 1.5
        Float broad_length_cutoff = 0.98
        Float conf_level = 0.75
        Int join_segment_size = 4
        Int max_sample_segs = 2500  # This is a bit low by default
        Int max_marker_spacing = 10000
        Int arm_peel = 0
        Int do_gene_gistic = 0
        Int remove_X = 1
        String gene_collapse_method = "mean"

        Int memoryMB = 2096
        Int disk_size = 10
        Int preemptible = 1
    }

    parameter_meta {
        seg_file: "A six-column tab-delimited file containing the segmented data for all tumor/normal pairs in the pair set."
        markers_file: "A three-column tab-delimited file identifying the names and positions all markers."
        refgene_file: "Contains information about the location of genes and cytobands on a given build of the genome.  These files are created in MATLAB."
        cnv_files: "A file specifying germ line CNVs to be excluded from the analysis."
        amp_thresh: "Threshold for copy number amplifications.  Regions with a log2 ratio above this value are considered amplified. (Recommended: 0.1)"
        del_thresh: "Threshold for copy number deletions.  Regions with a log2 ratio below the negative of this value are considered deletions.(Recommended: 0.1)"
        qv_thresh: "Significance threshold for Q-values.  Regions with Q-values below this number are considered signficant. (Recommended: 0.1)"
        cap: "Minimum and maximum cap values on analyzed data. Regions with a log2 ratio greater than the cap are set to the cap value; regions with a log2 ratio less than -cap value are set to -cap. (DEFAULT=1.5)"
        broad_length_cutoff: "Threshold used to distinguish broad from focal events, given in units of fraction of chromsome arm.  (Recommended: 0.7)"
        remove_X: "0/1 flag indicating whether to remove data from the X chromosome before analysis.  (Recommended: 0)"
        conf_level: "Confidence level used to calculate region containing the driver.  (Recommended: 0.99)"
        join_segment_size: "Smallest number of markers to allow in segements from the segemented data.  Segements that contain a number of markers less than or equal to this number are joined to the adjacent segement, closest in copy number. (Recommended: 4)"
        arm_peel: "0/1 flag indicating whether to perform arm-level peel off, wich helps separate peaks and clean up noise.  (Recommended: 1)"
        max_sample_segs: "Maximum number of segements allowed for a smaple in the input data.    Samples with more segments than this are excluded from the analysis. (Recommended: 2000)"
        do_gene_gistic: "0/1 flag indicating tht the gene GISTIC aglrithm should be used to calculate signficance of deletions at the gene level instead of a marker level. (Recommended: 1)"
        gene_collapse_method: "Method for reducing marker-level copy number data to the gene-level copy number data in the gene tables. Markers contained in the gene are used when available, otherwise the flanking marker or markers are used. Allowed values are mean, median, min, max or extreme. The extreme method chooses whichever of min or max is furthest from diploid. (Recommended: extreme)"
        memoryMB: "Integer value specifying the minimum memory requirements (in MB) for the virtual machine running the GISTIC2 task."
        preemptible: "Integer value specifying the maximum number of times Cromwell should request a preemptible machine for this task before defaulting back to a non-preemptible one."
    }

    command {
        set -euo pipefail

        if ~{!defined(cnv_files)} ; then
            echo -e "foo\t1\t1\t100\t1\t200" > dummy_cnv_file.txt
        fi

        # The link_conf_wrapper creates generic symlinks to files that specify
        # the confidence level in their names.
        /src/link_conf_wrapper.sh /src/call_gistic2 \
            . \
            4 \
            ~{seg_file} \
            ~{default="./this_file_does_not_exist.txt" markers_file} \
            ~{refgene_file} \
            ~{default="dummy_cnv_file.txt" cnv_files} \
            ~{amp_thresh} \
            ~{del_thresh} \
            ~{qv_thresh} \
            ~{cap} \
            ~{broad_length_cutoff} \
            ~{remove_X} \
            ~{conf_level} \
            ~{join_segment_size} \
            ~{arm_peel} \
            ~{max_sample_segs} \
            ~{do_gene_gistic} \
            ~{gene_collapse_method} \
            ~{max_marker_spacing} \
            ./version.txt
    }

    output {
        File all_data_by_genes = "all_data_by_genes.txt"
        File all_lesions = "all_lesions.txt"
        File all_thresholded_by_genes = "all_thresholded.by_genes.txt"
        File amp_genes = "amp_genes.txt"
        File amp_qplot_png = "amp_qplot.png"
        File arraylistfile = "arraylistfile.txt"
        File broad_significance_results = "broad_significance_results.txt"
        File broad_values_by_arm = "broad_values_by_arm.txt"
        File del_genes = "del_genes.txt"
        File del_qplot_png = "del_qplot.png"
        File gistic_inputs = "gisticInputs.txt"
        File gistic_version = "gisticVersion.txt"
        File segmented_copy_number_png = "raw_copy_number.png"
        File gistic_scores = "scores.gistic"
    }

    runtime {
        docker: "broadgdac/tool_gistic2:4"
        memory: memoryMB + " MB"
        disks: "local-disk " + disk_size + " HDD"
        preemptible: preemptible
    }
}
