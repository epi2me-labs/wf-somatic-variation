#!/usr/bin/env nextflow
import groovy.json.JsonBuilder
nextflow.enable.dsl = 2

// Load base modules
include {bam_ingress as bam_ingress_normal; bam_ingress as bam_ingress_tumor} from './lib/bamingress.nf'
include {
    somatic_sv
    } from './workflows/wf-somatic-sv.nf'
include {
    index_ref_gzi;
    index_ref_fai;
    cram_cache;
    decompress_ref;
    getAllChromosomesBed;
    getVersions;
    getParams
    } from './modules/local/common'
include {lookup_clair3_model; output_snv} from './modules/local/wf-somatic-snv'
include {alignment_stats; get_coverage; output_qc; discarded_sample} from './workflows/bamstats.nf'
include {snv} from './workflows/wf-somatic-snv'


// This is the only way to publish files from a workflow whilst
// decoupling the publish from the process steps.
process output {
    // publish inputs to output directory
    publishDir "${params.out_dir}", mode: 'copy', pattern: "*"
    input:
        path fname
    output:
        path fname
    """
    echo "Writing output files."
    """
}


// entrypoint workflow
WorkflowMain.initialise(workflow, params, log)
workflow {

    Map colors = NfcoreTemplate.logColours(params.monochrome_logs)

    if (workflow.profile.contains("conda")) {
        throw new Exception(colors.red + "Sorry, wf-human-variation is not compatible with --profile conda, please use --profile standard (Docker) or --profile singularity." + colors.reset)
    }

    can_start = true
    if (!params.snv && !params.sv && !params.methyl) {
        log.error (colors.red + "No work to be done! Choose one or more workflows to run from [--snv, --sv, --methyl]" + colors.reset)
        can_start = false
    }
    if (!params.bam_normal || !params.bam_tumor) {
        log.error (colors.red + "The workflow cannot run without passing both --bam_normal and --bam_tumor" + colors.reset)
        can_start = false
    }
    if (!file(params.bam_normal).exists()) {
        log.error (colors.red + "The workflow cannot run without passing a valid bam normal file" + colors.reset)
        can_start = false
    }
    if (!file(params.bam_tumor).exists()) {
        log.error (colors.red + "The workflow cannot run without passing a valid bam tumor file" + colors.reset)
        can_start = false
    }

    // Check ref and decompress if needed
    ref = null
    ref_index = null
    if (params.ref.toLowerCase().endsWith("gz")) {
        // gzipped ref not supported by some downstream tools (pyfaidx, cram_cache)
        // easier to just decompress and pass it around rather than confusing the user
        decompress_ref(file(params.ref))
        ref = decompress_ref.out.decompressed_ref
        ref_index = file("${ref}.fai").exists() ? Channel.of(file("${ref}.fai")) : index_ref_fai(ref).reference_index // Create fai index channel
    }
    else {
        ref = Channel.fromPath(params.ref, checkIfExists: true)
        ref_index = file(params.ref + ".fai").exists() ? Channel.of(file(params.ref + ".fai")) : index_ref_fai(ref).reference_index // Create fai index channel
    }

    // ************************************************************************
    // Bail from the workflow for a reason we should have already specified
    if (!can_start){
        throw new Exception("The workflow could not be started.")
    }
    // ************************************************************************

    // Dummy optional file
    OPTIONAL = file("$projectDir/data/OPTIONAL_FILE")


    if (!params.disable_ping) {
        Pinguscript.ping_post(workflow, "start", "none", params.out_dir, params)
    }

    // Build ref cache for CRAM steps that do not take a reference
    ref_cache = cram_cache(ref)

    // canonical ref and BAM channels to pass around to all processes
    ref_channel = ref.concat(ref_index).concat(ref_cache).buffer(size: 3)

    // Get software versions
    versions = getVersions()
    parameters = getParams()

    /*
    * Start processing the bam files
    * It accepts two bam files: 
    * 1. Control bam
    * 2. Cancer bam
    */
    bam_normal = bam_ingress_normal(
            ref,
            ref_index,
            params.bam_normal,
        )
    bam_tumor = bam_ingress_tumor(
            ref,
            ref_index,
            params.bam_tumor,
        )   
    all_bams = bam_normal
                .map{ bam, bai, meta -> 
                    meta.sample = params.sample_name
                    meta.type = 'normal'
                    [bam, bai, meta]
                }
                .mix(bam_tumor.map{bam, bai, meta -> 
                    meta.sample = params.sample_name
                    meta.type = 'tumor'
                    return [bam, bai, meta]
                })

    // Check input region bed file.
    // If it doesn't exists, then extract the regions from
    // the reference faidx file.
    bed = null
    default_bed_set = false
    if(params.bed){
        bed = Channel.fromPath(params.bed, checkIfExists: true)
    }
    else {
        default_bed_set = true
        bed = getAllChromosomesBed(ref_channel).all_chromosomes_bed
    }

    //Compute QC metrics and output QC statistics
    qcdata = alignment_stats(all_bams, ref_channel, bed, versions, parameters)
    output_qc( qcdata.outputs )

    // Apply bam coverage hard threshold to the pair
    // The dataset will fail if at least one of the bam has
    // coverage below the specified values. To account for different 
    // sequencing design, the two coverages are specified 
    // independently. 
    if (params.tumor_min_coverage > 0 || params.normal_min_coverage > 0){
        // Define if a dataset passes or not the filtering
        get_coverage(qcdata.coverages)
        
        // Branch filters on T/N class
        get_coverage.out.branch{
            tumor: it[1].type == 'tumor'
            normal: it[1].type == 'normal'
        } .set { branched_checks }

        // Cross the values and apply filter 
        branched_checks.normal
            .cross(branched_checks.tumor)
            .branch{
                pass: it[0][2] == "true" && it[1][2] == "true"
                not_pass: it[0][2] != "true" || it[1][2] != "true"
            }
            .set{ depth_filtered }

        // Combine the bam and branch them by whether they
        // pass the depth filter.
        depth_filtered.pass
            .map{normal, tumor -> normal[0]}
            .combine(all_bams.map{it -> [it[2].sample] + it}, by:0)
            .map{it ->
                it.size > 0 ? [it[1], it[2], it[3]] : it
            }
            .set{ pass_bam_channel }

        // If it doesn't pass the minimum depth required, 
        // emit a bam channel of discarded bam files.
        // Log out an error of failed coverage.
        // The method is much more convoluted to print only the type that is 
        // failing, whether it is normal or tumor, and their respective thresholds.
        depth_filtered.not_pass.map{normal, tumor -> normal}
            .mix(depth_filtered.not_pass.map{normal, tumor -> tumor}) |
            discarded_sample |
            subscribe {
                log.error "ERROR: Sample ${it[0]} (${it[1]}) will not be processed by the workflow as the detected coverage of ${it[2]}x is below the minimum coverage threshold of ${it[3]}x required for analysis."
            }
    } else {
        // If the bam_min_depth is 0, then run regardless.
        all_bams.set{ pass_bam_channel }
    }

    // Run snv workflow if requested
    if (params.snv) {
        // TODO: consider implementing custom modules in future releases.
        lookup_table = Channel
                            .fromPath("${projectDir}/data/clairs_models.tsv", checkIfExists: true)
                            .splitCsv(sep: '\t', header: true)
                            .map{ it -> [it.basecall_model_name, it.clair3_model_name] }
        clairs_model = Channel
                            .value(params.basecaller_cfg)
                            .cross(lookup_table) 
                            .filter{ caller, info -> info[1] != '-' }
                            .map{caller, info -> info[1] }
        lookup_table = Channel
                            .fromPath("${projectDir}/data/clair3_models.tsv", checkIfExists: true)
                            .splitCsv(sep: '\t', header: true)
                            .map{ it -> [it.basecall_model_name, it.clair3_model_name] }
        clair3_model = Channel
                            .value(params.basecaller_cfg)
                            .cross(lookup_table) 
                            .filter{ caller, info -> info[1] != '-' }
                            .map{caller, info -> info[1] }


        clair_vcf = snv(
            pass_bam_channel,
            bed,
            ref_channel,
            // mosdepth_stats,
            // bam_stats,
            clairs_model,
            clair3_model,
        )
        
        // Publish outputs in the appropriate folder
        clair_vcf | output_snv
   
    }
    
    // Start SV calling workflow
    if (params.sv){
        somatic_sv(pass_bam_channel, ref_channel, OPTIONAL)
    }

    // Emit reference and its index
    output(ref_channel.concat(versions).concat(parameters))

}

if (params.disable_ping == false) {
    workflow.onComplete {
        Pinguscript.ping_post(workflow, "end", "none", params.out_dir, params)
    }

    workflow.onError {
        Pinguscript.ping_post(workflow, "error", "$workflow.errorMessage", params.out_dir, params)
    }
}
