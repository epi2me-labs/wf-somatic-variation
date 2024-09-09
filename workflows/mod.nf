import groovy.json.JsonBuilder

process getVersions {
    label "wf_somatic_mod"
    cpus 1
    memory 4.GB
    output:
        path "versions_tmp.txt"
    script:
    """
    modkit --version | tr -s ' ' ',' >> versions_tmp.txt
    bgzip --version | awk 'NR==1 {print \$1","\$3}' >> versions_tmp.txt
    """
}


process rVersions {
    label "dss"
    cpus 2
    memory 4.GB
    input:
        path "versions_tmp.txt"
    output:
        path "versions.txt"
    script:
    """
    cat versions_tmp.txt > versions.txt
    R --version | awk 'NR==1 {print "R,"\$3}' >> versions.txt
    R -e "packageVersion('DSS')" | awk '\$1=="[1]" {print "DSS,"\$2}' >> versions.txt
    """
}


process getParams {
    cpus 1
    memory 4.GB
    output:
        path "params.json"
    script:
        def paramsJSON = new JsonBuilder(params).toPrettyString()
    """
    # Output nextflow params object to JSON
    echo '$paramsJSON' > params.json
    """
}


process validate_modbam {
    label "wf_common"
    input:
        tuple path(alignment), 
            path(alignment_index), 
            val(meta), 
            path(reference), 
            path(reference_index), 
            path(reference_cache), 
            env(REF_PATH)
    output:
        tuple path(alignment), 
            path(alignment_index), 
            val(meta), 
            path(reference), 
            path(reference_index), 
            path(reference_cache), 
            env(REF_PATH),
            env(valid)

    script:
    """
    valid=0
    workflow-glue check_valid_modbam ${alignment} || valid=\$?

    # Allow EX_OK and EX_DATAERR, otherwise explode
    if [ \$valid -ne 0 ] && [ \$valid -ne 65 ]; then
        exit 1
    fi
    """
}


// Pre-compute sample probabilities
process sample_probs {
    label "wf_somatic_mod"
    // Using 4 threads on a 90X takes ~30sec to complete
    memory { 8.GB * task.attempt - 1.GB }
    maxRetries 1
    errorStrategy = {task.exitStatus in [137,140] ? 'retry' : 'finish'}
    input:
        tuple path(alignment), 
            path(alignment_index), 
            val(meta), 
            path(reference), 
            path(reference_index), 
            path(reference_cache), 
            env(REF_PATH)
    output:
        tuple path(alignment), 
            path(alignment_index), 
            val(meta), 
            path(reference), 
            path(reference_index), 
            path(reference_cache), 
            env(REF_PATH),
            env(probs)

    script:
    // Set `--interval-size` to 5Mb to speed up sampling, and `--only-mapped -p 0.1` to be consistent with `pileup`
    """
    probs=\$( modkit sample-probs ${alignment} -p 0.1 --interval-size 5000000 --only-mapped --threads ${task.cpus} 2> /dev/null | awk 'NR>1 {ORS=" "; print "--filter-threshold "\$1":"\$3}' )
    """
}

process modkit {
    label "wf_somatic_mod"
    cpus params.modkit_threads
    memory {(1.GB * params.modkit_threads * task.attempt) + 3.GB}
    maxRetries 1
    errorStrategy {task.exitStatus in [137,140] ? 'retry' : 'finish'}
    input:
        tuple path(alignment), 
            path(alignment_index), 
            val(meta), 
            path(reference), 
            path(reference_index), 
            path(reference_cache), 
            env(REF_PATH),
            val(probs),
            path(bed),
            val(contig)
    output:
        tuple val(meta),
            val(contig),
            path("${meta.alias}.wf-somatic-mod.${meta.type}.${contig}.bedmethyl.gz")

    script:
    def options = params.force_strand ? '':'--combine-strands --cpg'
    if (params.modkit_args){
        options = "${params.modkit_args}"
    }
    """
    modkit pileup \\
        ${alignment} \\
        "${meta.alias}.wf-somatic-mod.${meta.type}.${contig}.bedmethyl" \\
        --ref ${reference} \\
        --region ${contig} \\
        --interval-size 1000000 \\
        --log-filepath modkit.log \\
        ${probs} \\
        --threads ${task.cpus} ${options}
    bgzip "${meta.alias}.wf-somatic-mod.${meta.type}.${contig}.bedmethyl"
    """
}

process concat_bedmethyl {
    label "wf_somatic_mod"
    cpus 4
    memory 8.GB

    input:
        tuple val(meta),
            val(contigs),
            path("bedmethyls/*")
    output:
        tuple val(meta),
            path("${params.sample_name}.wf-somatic-mod.${meta.type}.bedmethyl.gz")

    script:
    // Concatenate the bedMethyl, sort them and compress them
    // def label = group != '*' ? "${group}." : ""
    """
    zcat -f bedmethyls/* | \
        sort -k 1,1 -k 2,2n --parallel ${task.cpus} | \
        bgzip -c -@ ${task.cpus} > "${meta.alias}.wf-somatic-mod.${meta.type}.bedmethyl.gz"
    """
}

process bedmethyl_split {
    label "wf_common"
    cpus 1
    memory 4.GB
    input:
        tuple val(meta), 
            path(bed)
    output:
        tuple val(meta), 
            path("*.${meta.alias}.wf-somatic-mod.${meta.type}.bedmethyl.gz"), 
            emit: mod_outputs

    script:
    """
    workflow-glue mod_split ${bed}
    for i in *.bedmethyl; do
        bgzip \$i
    done
    """
}


process summary {
    label "wf_somatic_mod"
    cpus 4
    memory { 8.GB * task.attempt - 1.GB }
    maxRetries 1
    errorStrategy = {task.exitStatus in [137,140] ? 'retry' : 'finish'}
    input:
        tuple path(alignment), 
            path(alignment_index), 
            val(meta),
            path(reference), 
            path(reference_index), 
            path(reference_cache), 
            env(REF_PATH)
    output:
        tuple val(meta), 
            path("${meta.alias}.${meta.type}.mod_summary.tsv"), 
            emit: mod_summary

    script:
    // Modkit summary prints out a progress bar that cannot be avoided
    """    
    modkit summary -t ${task.cpus} ${alignment} | \
        awk 'BEGIN{OFS="\t"}; {print \$1,\$2,\$3,\$4,\$5,\$6}' > "${meta.alias}.${meta.type}.mod_summary.tsv"
    """
}

// Convert to DSS
process bed2dss {
    label "wf_somatic_mod"
    cpus 2
    memory 4.GB
    input:
        tuple val(meta), 
            val(mod),
            path("modified.bed.gz")

    output:
        tuple val(meta), 
            val(mod),
            path("*.${meta.alias}_${meta.type}.dss.tsv")

    shell:
    '''
    zcat modified.bed.gz | awk -v OFS='\t' 'BEGIN{{print "chr","pos","N","X"}}{{print $1,$2,$10,$12}}' > "!{mod}.!{meta.alias}_!{meta.type}.dss.tsv"
    '''
}

// Run DSS to compute DMR/L
process dss {
    label "dss"
    cpus { params.dss_threads <= 4 ? params.dss_threads : 4 }
    // Set memory to 16G/core + 2GB for buffer.
    // Benchmark shows that DSS ends up using more memory than predicted,
    // with spikes up to >74GB with 4 cores. In these cases we ignore the raised errors.
    memory { (task.cpus * 19.GB) }
    errorStrategy 'ignore'
    input:
        tuple val(meta), 
            val(mod),
            path("normal.bed"),
            path("tumor.bed")

    output:
        tuple val(meta), 
            val(mod),
            path("${meta.alias}.${mod}.dml.tsv"), 
            emit: dml,
            optional: true
        tuple val(meta), 
            val(mod),
            path("${meta.alias}.${mod}.dmr.tsv"), 
            emit: dmr,
            optional: true

    script:
    """
    #!/usr/bin/env Rscript
    library(DSS)
    require(bsseq)
    require(data.table)
    # Disable scientific notation
    options(scipen=999)

    # Import data and create input structure
    BSobj <- makeBSseqData( 
        list(
            fread("tumor.bed", sep = '\\t', header = T),
            fread("normal.bed", sep = '\\t', header = T)
        ), c("Tumor", "Normal")
    )

    # DML testing
    dmlTest = DMLtest(BSobj, 
        group1=c("Tumor"), 
        group2=c("Normal"),
        equal.disp = FALSE,
        smoothing=TRUE,
        smoothing.span=500,
        ncores=${task.cpus})
    # Compute DMLs
    dmls = callDML(dmlTest,
        delta=0.25,
        p.threshold=0.001)
    # Compute DMRs
    dmrs = callDMR(dmlTest,
        delta=0.25,
        p.threshold=0.001,
        minlen=100,
        minCG=5,
        dis.merge=1500,
        pct.sig=0.5)
    # Write output files
    write.table(dmls, '${meta.alias}.${mod}.dml.tsv', sep='\\t', quote=F, col.names=T, row.names=F)
    write.table(dmrs, '${meta.alias}.${mod}.dmr.tsv', sep='\\t', quote=F, col.names=T, row.names=F)
    """
}

// Make report.
process makeModReport {
    label "wf_common"
    cpus 1
    memory { 6.GB * task.attempt }
    maxRetries 1
    errorStrategy {task.exitStatus in [137,140] ? 'retry' : 'finish'}
    input: 
        tuple val(meta), 
            path("normal_summary/*"),
            path("tumor_summary/*"),
            path("DML/*"),
            path("DMR/*"),
            path("ref.fa.fai"),
            path("versions.txt"),
            path("params.json")

    output:
        tuple val(meta), path("${meta.alias}.wf-somatic-mod-report.html")

    script:
        def genome = meta.genome_build ? "--genome ${meta.genome_build} " : ""
        def diff_mod = params.diff_mod && params.bam_normal ? "--diff_mod" : ""
        """
        workflow-glue report_mod \\
            ${meta.alias}.wf-somatic-mod-report.html \\
            --normal_summary normal_summary/ \\
            --tumor_summary tumor_summary/ \\
            --dml DML/ \\
            --dmr DMR/ \\
            ${diff_mod} \\
            --reference_fai ref.fa.fai \\
            --sample_name ${meta.alias} \\
            --versions versions.txt \\
            --params params.json ${genome} \\
            --workflow_version ${workflow.manifest.version}
        """
}

process publish_modbase {
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

workflow mod {
    take:
        alignment
        reference
        chromosome_codes
    main:
        OPTIONAL = file("$projectDir/data/OPTIONAL_FILE")
        // Check for conflicting parameters
        if (params.force_strand && params.modkit_args){
            throw new Exception(colors.red + "You cannot use --force_strand with --modkit_args." + colors.reset)
        }
        if (params.modkit_args){
            log.warn "--modkit_args will override any preset we defined."
        }
        chroms = reference
        | map{ fasta, fai, cache, env -> fai }
        | splitCsv(sep:'\t')
        | map{ chrom, len, os1, os2, os3 -> chrom }

        if (!params.include_all_ctgs){
            chroms = chroms | filter{ it in chromosome_codes }
        }

        // Call modified bases
        alignment.combine( reference ) | validate_modbam

        // Warn of invalid bam files
        validate_modbam.out.branch{
            stdbam: it[-1] == '65'
            modbam: it[-1] == '0'
            }.set{validated_bam}
        validated_bam.stdbam.subscribe{
            it -> log.warn "Input ${it[1]} is not a modified bam file."
        }

        // Define input bed, if provided
        if (params.bed){
            bed = Channel.fromPath(params.bed)
        } else {
            bed = Channel.fromPath("$projectDir/data/OPTIONAL_FILE")
        }

        // First, compute the sample_probs
        modkit_probs = sample_probs(validated_bam.modbam.map{it[0..-2]})
        | combine(bed)
        | combine(chroms)

        // Run modkit on the valid files.
        modkit_outputs = modkit(modkit_probs) | groupTuple(by:0) | concat_bedmethyl

        // Split modkit, and use file to rename the outputs
        bedmethyl_split(modkit_outputs)
        // Each channel has a nested tuple. Transpose linearize it,
        // and then we extract the modification type with simpleName.
        modbed = bedmethyl_split.out.mod_outputs
            .transpose()
            .map{ meta, tsv -> 
                [meta, tsv.simpleName, tsv]
            }

        // Compute summary stats
        alignment.combine( reference ) | summary

        // If required, compute the differentially modified regions/loci with DSS
        if (params.diff_mod){
            // Convert to DSS
            dss_inputs = modbed | bed2dss
            
            // Combine the outputs to perform DMR analyses
            dss_inputs.branch{
                tumor: it[0].type == 'tumor'
                normal: it[0].type == 'normal'
            }.set{forked_mb}

            // Output methylation data
            forked_mb.normal
                .map{it -> [it[0].sample, it[1], it[2], it[0]]}
                .combine(
                    forked_mb.tumor.map{it -> [it[0].sample, it[1], it[2], it[0]]}, 
                    by:[0,1]
                )
                .map{sample, mod, norm_bed, norm_meta, tum_bed, tum_meta -> [tum_meta, mod, norm_bed, tum_bed ] }
                .set{ paired_beds }
            paired_beds | dss
            dml_ch = dss.out.dml
            dmr_ch = dss.out.dmr
        } else {
            dss_inputs = Channel.empty()
        }

        // Get versions and params
        software_versions = getVersions() | rVersions
        workflow_params = getParams()

        // Make report
        forked_sum = summary.out.mod_summary.branch{
            tumor: it[0].type == 'tumor'
            normal: it[0].type == 'normal'
        }

        // If the normal bam is provided, pass the normal output to the reporting process.
        if (params.bam_normal){
            // Combine the summaries first.
            forked_sum.normal
                | map{meta, summary -> [meta.alias, summary, meta]}
                | combine(
                    forked_sum.tumor.map{
                        meta, summary ->
                        [meta.alias, summary, meta]
                    }, by:0
                )
                | map{
                    sample, norm_summary, norm_meta, tum_summary, tum_meta ->
                    [tum_meta, norm_summary, tum_summary]
                }
                | set{combined_summaries}
        // Otherwise, pass OPTIONAL_FILE
        } else {
            // Combine the summaries first.
            forked_sum.tumor
                | map{
                    t_meta, t_summary ->
                    [t_meta, OPTIONAL, t_summary]
                }
                | set{combined_summaries}
        }

        // If DSS is run, create a channel with the results in it for the report.
        if (params.diff_mod && params.bam_normal){
            dss_outputs = dml_ch
                | combine(dmr_ch, by: [0,1])
                | groupTuple(by:0)            
                // Allow remainder to enable merging an empty tuple with the summaries
                // allowing to then generate the report at the end
                | join(combined_summaries, by: 0, remainder: true)
                | combine(reference)
                // if all DSS processes fail, the tuple will start with a metadata, have
                // a null value, and then the rest of the tuple. Here we check if that is
                // the case, and if so add optional files to fill in the gaps.
                | map{
                    it[1] ? it : [it[0], null, OPTIONAL, OPTIONAL] + it[2..-1]
                }
                | map{
                    meta, mod, dms, dmr, sum_n, sum_t, ref, fai, cache, ref_path -> 
                    [meta, sum_n, sum_t, dms, dmr, fai]
                }
                | set{ for_report }
        // Otherwise, ignore DSS outputs and pass the summaries (or OPTIONAL_FILE) to the report.
        } else {
            combined_summaries
                | combine(reference)
                | map{
                    meta, n_summary, t_summary, ref, fai, cache, ref_path ->                     
                    [meta, n_summary, t_summary, OPTIONAL, OPTIONAL, fai]
                }
                | set{ for_report }
        }

        // Create output report
        for_report
        | combine(software_versions)
        | combine(workflow_params)
        | makeModReport

        // Create output directory
        // Save raw bedMethyl and summaries in main folder,
        // and the rest in the {{ sample }}/mod directory.
        artifacts = modkit_outputs
        | map{
                meta, file -> [file, null]
        }
        | mix(
            modbed.map{
                meta, mod, file -> [file, "${meta.alias}/mod/${mod}/bedMethyl/"]
            }
        )
        | mix(
            summary.out.mod_summary.map{
                meta, summary -> [summary, null]
            }
        )
        | mix(
            dss_inputs.map{
                meta, mod, file -> [file, "${meta.alias}/mod/${mod}/DSS/"]
            }
        ) 
        | mix(
            workflow_params.map{
                params -> [params, "info/mod/"]
            }
        )
        | mix(
            software_versions.map{
                versions -> [versions, "info/mod/"]
            }
        )
        | mix(
            makeModReport.out.map{
                meta, report -> [report, null]
            }
        )
        if (params.diff_mod){
            artifacts = artifacts
            | mix(
                dml_ch.map{
                    meta, mod, file -> [file, "${meta.alias}/mod/${mod}/DML"]
                }
            )
            | mix(
                dmr_ch.map{
                    meta, mod, file -> [file, "${meta.alias}/mod/${mod}/DMR"]
                }
            )
        }
        artifacts | publish_modbase

    emit:
        modbam2bed = bedmethyl_split.out.mod_outputs
        dss = dss_inputs
        mod_summaries = combined_summaries
        report_mod = makeModReport.out
}
