process bamstats {
    cpus 4
    input:
        tuple path(xam), path(xam_idx), val(xam_meta)
        tuple path(ref), path(ref_idx), path(ref_cache), env(REF_PATH)

    output:
        tuple val(xam_meta), path("*.readstats.tsv.gz"), emit: read_stats
        tuple val(xam_meta), path("*.flagstat.tsv"), emit: flagstat
    script:
    def cores = task.cpus > 1 ? task.cpus - 1 : 1
    """
    bamstats ${xam} -s ${xam_meta.sample} --threads ${cores} -u -f ${xam_meta.sample}_${xam_meta.type}.flagstat.tsv | gzip > ${xam_meta.sample}_${xam_meta.type}.readstats.tsv.gz
    """
}

process mosdepth {
    cpus 2
    input:
        tuple path(xam), path(xam_idx), val(xam_meta)
        file target_bed
        tuple path(ref), path(ref_idx), path(ref_cache), env(REF_PATH)
    output:
        tuple val(xam_meta), \
            path("${xam_meta.sample}_${xam_meta.type}.regions.bed.gz"),
            path("${xam_meta.sample}_${xam_meta.type}.mosdepth.global.dist.txt"),
            path("${xam_meta.sample}_${xam_meta.type}.thresholds.bed.gz"), emit: mosdepth_tuple
        tuple val(xam_meta), path("${xam_meta.sample}_${xam_meta.type}.mosdepth.summary.txt"), emit: summary
        tuple val(xam_meta), path("${xam_meta.sample}_${xam_meta.type}.per-base.bed.gz"), emit: perbase, optional: true
    script:
        def perbase_args = params.depth_intervals ? "" : "--no-per-base"
        """
        export REF_PATH=${ref}
        export MOSDEPTH_PRECISION=3
        # Convert bed into windows of given size [CW-2015]
        bedtools makewindows -b ${target_bed} -w ${params.depth_window_size} > cut.bed
        # Run mosdepth
        mosdepth \\
            -x \\
            -t ${task.cpus} \\
            -b cut.bed \\
            --thresholds 1,10,20,30 \\
            ${perbase_args} \\
            ${xam_meta.sample}_${xam_meta.type} \\
            ${xam}
        """
}


// Get coverage to a channel
process get_coverage {
    cpus 1
    input:
        tuple val(meta), path(mosdepth_summary)

    output:
        tuple val(meta.sample), val(meta), env(passes), env(value), emit: pass

    shell:
        '''
        # Check if the first column is "total_region", skipping the header (NR>1) and if the value is above the 
        # threshold for the right type. This is defined with the elvis ops checking the metadata for the bam, getting 
        # the file type (tumor or normal) and, therefore, the right threshold. If it is >= than the threshold, return "true"
        # otherwise return "false".
        passes=$( awk 'BEGIN{v="false"}; NR>1 && $1=="total_region" && $4>=!{meta.type == "tumor" ? params.tumor_min_coverage : params.normal_min_coverage} && v=="false" {v="true"}; END {print v}' !{mosdepth_summary} )

        # Same as above, but simply return the coverage value for the bam file in the "total_region".
        value=$( awk 'BEGIN{v=0}; NR>1 && $1=="total_region" && $4>v {v=$4}; END {print v}' !{mosdepth_summary} )
        '''
}

// Get coverage to a channel
process discarded_sample {
    cpus 1
    input:
        tuple val(sample), val(meta), val(passing), val(coverage)

    output:
        tuple val(sample), val(meta), val(coverage), env(threshold), emit: failed

    shell:
        '''
        threshold=!{meta.type =='tumor' ? params.tumor_min_coverage : params.normal_min_coverage}
        '''
}

// Make report
process makeQCreport {
    input: 
        tuple val(meta), 
            path("readstats_normal.tsv.gz"),
            path("flagstat_normal.tsv"),
            path("summary_depth_normal.tsv"),
            path("depth_normal.tsv.gz"),
            path("readstats_tumor.tsv.gz"),
            path("flagstat_tumor.tsv"),
            path("summary_depth_tumor.tsv"),
            path("depth_tumor.tsv.gz"),
            path("ref.fa.fai")
        path "versions.txt"
        path "params.json"

    output:
        tuple val(meta), path("${meta.sample}.wf-somatic-variation-readQC*.html")

    script:
        // If no *_min_coverage provided, or set to null by mistake, set it to 0.
        def tumor_cvg = params.tumor_min_coverage ?: 0
        def normal_cvg = params.normal_min_coverage ?: 0
        """
        workflow-glue report_qc \\
            --window_size ${params.depth_window_size} \\
            --tumor_cov_threshold ${tumor_cvg} \\
            --normal_cov_threshold ${normal_cvg} \\
            --sample_id ${meta.sample} \\
            --name ${meta.sample}.wf-somatic-variation-readQC \\
            --read_stats_normal readstats_normal.tsv.gz \\
            --read_stats_tumor readstats_tumor.tsv.gz \\
            --flagstat_tumor flagstat_tumor.tsv \\
            --flagstat_normal flagstat_normal.tsv \\
            --mosdepth_summary_tumor summary_depth_tumor.tsv \\
            --mosdepth_summary_normal summary_depth_normal.tsv \\
            --depth_tumor depth_tumor.tsv.gz \\
            --depth_normal depth_normal.tsv.gz \\
            --reference_fai ref.fa.fai \\
            --versions versions.txt \\
            --params params.json
        """
}


process output_qc {
    // publish inputs to output directory
    publishDir (
        params.out_dir,
        mode: "copy",
        saveAs: { dirname ? "$dirname/$fname" : fname }
    )
    input:
        tuple path(fname), val(dirname)
    output:
        path fname
    """
    """
}


workflow alignment_stats {
    take:
        bamfiles
        ref
        bed
        versions
        parameters
    
    main:
        // Compute bam statistics and depth
        stats = bamstats(bamfiles, ref.collect())
        depths = mosdepth(bamfiles, bed.collect(), ref.collect())

        // Combine the outputs for the different statistics.
        // Fo the reporting we will need:
        // 1. Read stats
        // 2. flagstats
        // 3. Per-base depth
        // 4. Depth summary
        stats.read_stats.map{meta, rstats->[meta, rstats]}
            .combine(stats.flagstat.map{meta, fstats->[meta, fstats]}, by:0)
            .combine(depths.summary.map{meta, depth_sum->[meta, depth_sum]}, by:0)
            .combine(depths.mosdepth_tuple.map{meta, reg, dist, thresh ->[meta, reg]}, by:0)
            .set{ for_report }

        // Cross the results for T/N pairs
        for_report
            .branch{
                tumor: it[0].type == 'tumor'
                normal: it[0].type == 'normal'
            }
            .set{forked_channel}
        forked_channel.normal
            .map{ it -> [ it[0].sample ] + it } 
            .cross(
                forked_channel.tumor.map{ it -> [ it[0].sample ] + it } 
            )
            .map { normal, tumor ->
                    [tumor[1]] + normal[2..-1] + tumor[2..-1]
                } 
            .combine(ref.map{it[1]})
            .set{paired_samples}
        
        makeQCreport(paired_samples, versions, parameters)

        // Prepare output channel
        // Send the output to the specified sub-directory of params.out_dir.
        // If null is passed, send it to out_dir/ directly.
        if (params.depth_intervals){
            makeQCreport.out.map{it -> [it[1], null]}
                .concat(stats.flagstat.map{meta, fstats -> [fstats, "${meta.sample}/qc/readstats"]})
                .concat(stats.read_stats.map{meta, rstats -> [rstats, "${meta.sample}/qc/readstats"]})
                .concat(depths.summary.map{meta, depth_sum -> [depth_sum, "${meta.sample}/qc/coverage"]})
                .concat(depths.mosdepth_tuple
                            .map {meta, reg, dist, thresh  -> [meta, [reg, dist, thresh]] }
                            .transpose()
                            .map{meta, fname -> [fname, "${meta.sample}/qc/coverage"]})
                .concat(depths.perbase.map{meta, pbase ->[pbase, "${meta.sample}/qc/coverage"]})
                .set{outputs}
        } else {
            makeQCreport.out.map{meta, report -> [report, null]}
                .concat(stats.flagstat.map{meta, fstats -> [fstats, "${meta.sample}/qc/readstats"]})
                .concat(stats.read_stats.map{meta, rstats -> [rstats, "${meta.sample}/qc/readstats"]})
                .concat(depths.summary.map{meta, depth_sum -> [depth_sum, "${meta.sample}/qc/coverage"]})
                .concat(depths.mosdepth_tuple
                            .map {meta, reg, dist, thresh  -> [meta, [reg, dist, thresh]] }
                            .transpose()
                            .map{meta, fname -> [fname, "${meta.sample}/qc/coverage"]})
                .set{outputs}
        }

        emit:
            outputs = outputs
            coverages = depths.summary
            paired_qc = paired_samples
}