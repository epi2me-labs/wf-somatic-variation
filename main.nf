#!/usr/bin/env nextflow
import groovy.json.JsonBuilder
nextflow.enable.dsl = 2

// Load base modules
include {bam_ingress as bam_ingress_control; bam_ingress as bam_ingress_cancer} from './lib/bamingress.nf'
include {
    index_ref_gzi;
    index_ref_fai;
    cram_cache;
    decompress_ref;
    getAllChromosomesBed
    } from './modules/local/common'
include {lookup_clair3_model; output_snp} from './modules/local/wf-somatic-snp'
include {snp} from './workflows/wf-somatic-snp'


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
    if (!params.snp && !params.sv && !params.methyl) {
        log.error (colors.red + "No work to be done! Choose one or more workflows to run from [--snp, --sv, --methyl]" + colors.reset)
        can_start = false
    }
    if (!params.bam_control || !params.bam_cancer) {
        log.error (colors.red + "The workflow cannot run without passing both --bam_control and --bam_cancer" + colors.reset)
        can_start = false
    }
    if (!file(params.bam_control).exists()) {
        log.error (colors.red + "The workflow cannot run without passing a valid bam control file" + colors.reset)
        can_start = false
    }
    if (!file(params.bam_cancer).exists()) {
        log.error (colors.red + "The workflow cannot run without passing a valid bam cancer file" + colors.reset)
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

    /*
    * Start processing the bam files
    * It accepts two bam files: 
    * 1. Control bam
    * 2. Cancer bam
    */
    bam_control = bam_ingress_control(
            ref,
            ref_index,
            params.bam_control,
        )
    bam_cancer = bam_ingress_cancer(
            ref,
            ref_index,
            params.bam_cancer,
        )   
    all_bams = bam_control
                .map{ bam, bai, meta -> 
                    meta.sample = params.sample_name
                    meta.type = 'control'
                    [bam, bai, meta]
                }
                .mix(bam_cancer.map{bam, bai, meta -> 
                    meta.sample = params.sample_name
                    meta.type = 'cancer'
                    return [bam, bai, meta]
                })


    // wf-somatic-snp
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
    // Run snp workflow if requested
    if (params.snp) {
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


        clair_vcf = snp(
            all_bams,
            bed,
            ref_channel,
            // mosdepth_stats,
            // bam_stats,
            clairs_model,
            clair3_model,
        )
        
        // Publish outputs in the appropriate folder
        clair_vcf | output_snp
    }

    // Emit reference and its index
    output(ref_channel)

}

if (params.disable_ping == false) {
    workflow.onComplete {
        Pinguscript.ping_post(workflow, "end", "none", params.out_dir, params)
    }

    workflow.onError {
        Pinguscript.ping_post(workflow, "error", "$workflow.errorMessage", params.out_dir, params)
    }
}
