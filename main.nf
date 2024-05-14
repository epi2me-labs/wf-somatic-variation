#!/usr/bin/env nextflow
import groovy.json.JsonBuilder
nextflow.enable.dsl = 2

// Load base modules
include {
    ingress as ingress_normal;
    ingress as ingress_tumor
    } from './lib/_ingress.nf'
include {
    somatic_sv as sv
    } from './workflows/wf-somatic-sv.nf'
include {
    index_ref_gzi;
    index_ref_fai;
    cram_cache;
    decompress_ref;
    getAllChromosomesBed;
    getVersions;
    getVersions_somvar;
    getParams;
    getGenome
    } from './modules/local/common'
include {lookup_clair3_model; output_snv} from './modules/local/wf-somatic-snv'
include {
    alignment_stats; get_coverage; get_region_coverage; 
    output_qc; get_shared_region
} from './workflows/bamstats.nf'
include {snv} from './workflows/wf-somatic-snv'
include {mod} from './workflows/mod.nf'


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
    if (!params.snv && !params.sv && !params.mod) {
        log.error (colors.red + "No work to be done! Choose one or more workflows to run from [--snv, --sv, --mod]" + colors.reset)
        can_start = false
    }
    if (!file(params.bam_tumor).exists()) {
        log.error (colors.red + "The workflow cannot run without passing a valid bam tumor file" + colors.reset)
        can_start = false
    }
    if (!params.bam_normal && (params.snv || params.sv)) {
        log.error (colors.red + "The tumor-only mode is only available with --mod" + colors.reset)
        can_start = false
    }
    if (params.bam_normal && !file(params.bam_normal).exists()){
        log.error (colors.red + "The workflow cannot run without passing a valid bam normal file" + colors.reset)
        can_start = false
    }
    if (!params.germline) {
        log.warn ("The workflow is running in somatic-only mode, germline calling will be skipped")
    }
    if (params.normal_vcf) {
        if (!file("${params.normal_vcf}", checkifExists: true)){
            throw new Exception("--normal_vcf is specified, but the file doesn't exist: ${params.normal_vcf}")
        }
        if (params.normal_vcf.endsWith('.gz') && !file("${params.normal_vcf}.tbi", checkifExists: true)){
            throw new Exception("No TBI index for VCF file: ${params.normal_vcf}")
        }
        log.info ("Pre-computed VCF for the normal sample provided; running germline calling only for tumor sample")
    }
    if (params.snv && params.hybrid_mode_vcf && params.genotyping_mode_vcf){
        throw new Exception("Can run --hybrid_mode_vcf or --genotyping_mode_vcf, not both. Choose one and try again.")
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


    Pinguscript.ping_start(nextflow, workflow, params)

    // Build ref cache for CRAM steps that do not take a reference
    ref_cache = cram_cache(ref)
    ref_cache = cram_cache.out.ref_cache
    ref_path = cram_cache.out.ref_path
    // canonical ref and BAM channels to pass around to all processes
    ref_channel = ref.concat(ref_index).concat(ref_cache).concat(ref_path).buffer(size: 4)

    // Get software versions
    versions = getVersions() | getVersions_somvar
    parameters = getParams()

    /*
    * Start processing the bam files
    * It accepts two bam files: 
    * 1. Tumor bam
    * 2. Control bam (optional for mod)
    */
    // If running in tumor-only mode, create an empty channel.
    if (params.bam_normal){
        bam_normal = ingress_normal(
                ref,
                ref_index,
                params.bam_normal,
            )
    } else {
        bam_normal = Channel.empty()
    }

    // Import the tumor, which is always required.
    bam_tumor = ingress_tumor(
            ref,
            ref_index,
            params.bam_tumor,
        )

    // Combine everything
    all_bams = bam_normal
                .map{
                    meta, bam, bai -> 
                    [bam, bai, meta + [sample: params.sample_name, type: 'normal']]
                }
                .mix(
                    bam_tumor.map{
                        meta, bam, bai -> 
                        [bam, bai, meta + [sample: params.sample_name, type: 'tumor']]
                        }
                    )

    // Add genome build information
    // CW-2491: make this optional, allowing any genome to be processed
    // CW-3830: perform this before the QC as changing the metadata causes the
    //          `-resume` to break.
    if (params.annotation){
        getGenome(all_bams)
        getGenome.out.genome_build.map{
                bam, bai, meta, g_build -> 
                    [bam, bai, meta + [genome_build: g_build]]
            }.set{all_bams}
    } else {
        all_bams
            .map{
                bam, bai, meta ->
                [bam, bai, meta + [genome_build: null]]
            }
            .set{all_bams}
    }

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
    qc_outputs = qcdata.outputs

    // Apply bam coverage hard threshold to the pair
    // The dataset will fail if at least one of the bam has
    // coverage below the specified values. To account for different 
    // sequencing design, the two coverages are specified 
    // independently. 
    if (params.tumor_min_coverage > 0 || params.normal_min_coverage > 0){
        // Define if a dataset passes or not the filtering
        if (params.bed){
            // Filter out the data based on the individual region's coverage
            coverage_check = qcdata.mosdepth_tuple.combine(bed) | get_region_coverage
            // Unlike humvar, the paired nature of the wf requires to branch and 
            // match the input bed files to detect the shared regions between T/N
            // first we separate tumor and normal
            coverage_check.filt_bed.branch{
                tumor: it[0].type == 'tumor'
                normal: it[0].type == 'normal'
            }.set{branched_filtered_beds}

            // Then, we cross T/N based on the sample name
            filt_bed = branched_filtered_beds.tumor
                .map{
                    meta, bed -> [meta.sample, bed]
                }.combine(
                    branched_filtered_beds.normal.map{
                        meta, bed -> [meta.sample, bed]
                    }, by: 0
                ) | get_shared_region // and then process them to intersect the retained regions

            qc_outputs
                .mix(
                    filt_bed.bed_tuple.map{
                        sample, fname -> [fname, "${sample}/qc/coverage"]
                    }
                )
                .mix(
                    coverage_check.mosdepth_tuple.map{
                        meta, filt_cvg, dist, threshold -> [filt_cvg, "${meta.sample}/qc/coverage"]
                    }
                )
                .set{ qc_outputs }
            
            // Replace the bed with the filtered bed
            bed = filt_bed.bed_file
        } else {
            // Define if a dataset passes or not the filtering
            coverage_check = get_coverage(qcdata.coverages)
        }


        // Branch filters on T/N class
        coverage_check.pass.branch{
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
        depth_filtered.not_pass
            .flatMap()
            .map{
                sample, meta, passing, coverage ->
                def threshold = meta.type =='tumor' ? params.tumor_min_coverage as float : params.normal_min_coverage as float
                def logged = coverage as float < threshold ? 
                "will not be processed by the workflow as the detected coverage of ${coverage}x is below the minimum coverage threshold of ${threshold}x required for analysis" : 
                "will not be processed as the matching bam is below the minimum coverage threshold required for analysis"
                [sample, meta, coverage, threshold, logged]
            }
            .subscribe {
                log.error "ERROR: Sample ${it[0]} (${it[1].type}) ${it[4]}."
            }
    } else {
        // If the bam_min_depth is 0, then run regardless.
        all_bams.set{ pass_bam_channel }
    }
    // Output QC data
    output_qc( qc_outputs )

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
            clairs_model,
            clair3_model
        )
        
        // Publish outputs in the appropriate folder
        clair_vcf | output_snv
   
    }
    
    // Start SV calling workflow
    if (params.sv){
        sv(pass_bam_channel, ref_channel, OPTIONAL)
    }

    // Extract modified bases
    if (params.mod){
        mod(
            pass_bam_channel,
            ref_channel,
        )
    }

    // Emit version and parameters
    output(versions.concat(parameters))


}

workflow.onComplete {
    Pinguscript.ping_complete(nextflow, workflow, params)
}
workflow.onError {
    Pinguscript.ping_error(nextflow, workflow, params)
}
