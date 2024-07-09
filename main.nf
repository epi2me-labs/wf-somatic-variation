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
    getGenome;
    report
    } from './modules/local/common'
include {lookup_clair3_model; publish_snv} from './modules/local/wf-somatic-snv'
include {
    alignment_stats; get_coverage; get_region_coverage; 
    publish_qc; get_shared_region
} from './workflows/bamstats.nf'
include {snv} from './workflows/wf-somatic-snv'
include {snv as snv_to} from './workflows/wf-somatic-snv-to'
include {mod} from './workflows/mod.nf'


// This is the only way to publish files from a workflow whilst
// decoupling the publish from the process steps.
process publish {
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

    def run_tumor_only = !params.bam_normal 
    can_start = true
    if (!params.snv && !params.sv && !params.mod) {
        log.error (colors.red + "No work to be done! Choose one or more workflows to run from [--snv, --sv, --mod]" + colors.reset)
        can_start = false
    }
    if (!file(params.bam_tumor).exists()) {
        log.error (colors.red + "The workflow cannot run without passing a valid bam tumor file" + colors.reset)
        can_start = false
    }
    if (run_tumor_only && params.sv) {
        log.error (colors.red + "The tumor-only mode is not available with --sv" + colors.reset)
        can_start = false
    }
    if (run_tumor_only && params.snv && params.liquid_tumor) {
        log.warn "The SNV tumor-only mode currently has no specific presets for liquid tumors."
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

    // Programmatically define chromosome codes.
    // note that we avoid interpolation (eg. "${chr}N") to ensure that values
    // are Strings and not GStringImpl, ensuring that .contains works.
    ArrayList chromosome_codes = []
    ArrayList chromosomes = [1..22] + ["X", "Y", "M", "MT"]
    for (N in chromosomes.flatten()){
        chromosome_codes += ["chr" + N, "" + N]
    }

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

    // populate output json with ingressed runids
    ArrayList ingressed_run_ids = []
    qcdata.runids.splitText().subscribe(
        onNext: {
            ingressed_run_ids += it.strip()
        },
        onComplete: {
            params.wf["ingress.run_ids"] = ingressed_run_ids
        }
    )

    // Apply bam coverage hard threshold to the pair
    // The dataset will fail if at least one of the bam has
    // coverage below the specified values. To account for different 
    // sequencing design, the two coverages are specified 
    // independently.
    if (params.tumor_min_coverage > 0 || (params.bam_normal && params.normal_min_coverage > 0)){
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
            if (run_tumor_only){
                // If no normal BAM is provided, use the tumor outputs only.
                filt_bed = branched_filtered_beds.tumor
                    .map{
                        meta, bed -> [meta.sample, bed]
                    }
                bed = branched_filtered_beds.tumor
                    .map{
                        meta, bed -> bed
                    }
            } else {
                branched_filtered_beds.tumor
                    .map{
                        meta, bed -> [meta.sample, bed]
                    }.combine(
                        branched_filtered_beds.normal.map{
                            meta, bed -> [meta.sample, bed]
                        }, by: 0
                    ) | get_shared_region // and then process them to intersect the retained regions
                filt_bed = get_shared_region.out.bed_tuple
                // Prepare the filtered bed.
                bed = get_shared_region.out.bed_file
            }

            // Add more outputs
            qc_outputs
                .mix(
                    filt_bed.map{
                        sample, fname -> [fname, "${sample}/qc/coverage"]
                    }
                )
                .mix(
                    coverage_check.mosdepth_tuple.map{
                        meta, filt_cvg, dist, threshold -> [filt_cvg, "${meta.sample}/qc/coverage"]
                    }
                )
                .set{ qc_outputs }
        } else {
            // Define if a dataset passes or not the filtering
            coverage_check = get_coverage(qcdata.coverages)
        }


        // Branch filters on T/N class
        coverage_check.pass.branch{
            tumor: it[1].type == 'tumor'
            normal: it[1].type == 'normal'
        } .set { branched_checks }

        // Apply filters depending on the presence of normal or not
        if (run_tumor_only){
            // Apply filter
            branched_checks.tumor
                .branch{
                    sample, meta, pass, value ->
                    pass: pass == "true"
                    fail: true
                }
                .set{ depth_filtered }
            // Create temporary pass channel keeping only the sample name.
            tmp_pass_ch = depth_filtered.pass.map{ sample, meta, pass, value -> sample }
            // Define non-passing channel
            // Tumor-only is not a nested tuple, so avoid flatmapping
            fail_bam_channel = depth_filtered.fail
        } else {
            // Cross the values and apply filter.
            // Crossing creates a nested tuple of [ normal, tumor ], where
            // where normal is a tuple with structure:
            //   [
            //     sample ID,
            //     meta,
            //     boolean for whether the bam passes coverage check,
            //     coverage value
            //   ]
            // When branching, we therefore check the normal and the tumor
            // both pass the filtering threshold.
            branched_checks.normal
                .cross(branched_checks.tumor)
                .branch{
                    normal, tumor ->
                    pass: normal[2] == "true" && tumor[2] == "true"
                    fail: true
                }
                .set{ depth_filtered }
            // Create temporary pass channel keeping only the sample name.
            tmp_pass_ch = depth_filtered.pass.map{normal, tumor -> normal[0]}
            // Define non-passing channel
            // Being a nested tuple, we need to flatMap it first
            fail_bam_channel = depth_filtered.fail.flatMap()
        }
        // Add the bam and branch them based on passing/failing
        // the depth filter.
        tmp_pass_ch
            .combine(all_bams.map{it -> [it[2].sample] + it}, by:0)
            .map{ sample, xam, xai, meta -> [xam, xai, meta] }
            .set{ pass_bam_channel }

        // If it doesn't pass the minimum depth required, 
        // emit a bam channel of discarded bam files.
        // Log out an error of failed coverage.
        // The method is much more convoluted to print only the type that is 
        // failing, whether it is normal or tumor, and their respective thresholds.
        fail_bam_channel
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
    publish_qc( qc_outputs )

    // Create minimal channel for joint report
    for_joint_report = qcdata.report_qc

    // Run snv workflow if requested
    if (params.snv) {
        // Add back basecaller models, if available.
        // Combine each BAM channel with the appropriate basecaller file
        // Fetch the unique basecaller models and, if these are more than the
        // ones in the metadata, add them in there.
        // We do it in the snv scope as it is the only workflow relying on the
        // model, and given it has to wait for the readStats process, we try
        // minimizing the waits
        pass_bam_channel = pass_bam_channel 
        | map{ xam, xai, meta -> [meta, xam, xai] }
        | combine( qcdata.basecallers, by:0 )
        | map{
            meta, xam, xai, bc ->
            def models = bc.splitText().collect { it.strip() }
            [xam, xai, meta + [basecall_models: models]]
        }
    
        // attempt to pull out basecaller_cfg from metadata, keeping unique values
        metamap_basecaller_cfg = pass_bam_channel
            | map { xam, bai, meta ->
                meta["basecall_models"]
            }
            | flatten  // squash lists
            | unique

        // Ensure that the two BAM have the same basecaller configuration.
        // Check that there is exactly one metamap.
        // Specific for somvar: fail also with >1 metadata, as the two BAM need the
        // same model.
        metamap_basecaller_cfg
                | count
                | map { int n_models ->
                    if (n_models == 0){
                        if (params.basecaller_cfg) {
                            log.warn "Found zero basecall_model in the input alignment header, falling back to the model provided with --basecaller_cfg: ${params.basecaller_cfg}"
                        }
                        else {
                            String input_data_err_msg = '''\
                            ################################################################################
                            # INPUT DATA PROBLEM
                            Your input alignments does not indicate the basecall model in the header and
                            you did not provide an alternative with --basecaller_cfg.

                            wf-somatic-variation requires the basecall model in order to automatically
                            select an appropriate SNP calling model.

                            ## Next steps
                            You must re-run the workflow specifying the basecaller model with the
                            --basecaller_cfg option.
                            ################################################################################
                            '''.stripIndent()
                            error input_data_err_msg
                        }
                    } else if (n_models > 1){
                        String input_data_err_msg = '''\
                        ################################################################################
                        # INPUT DATA PROBLEM
                        Your input tumor and/or normal BAM files indicate two different basecall models
                        in the header.

                        wf-somatic-variation requires the basecall model used to call the tumor and 
                        normal samples to be identical.

                        ## Next steps
                        You must re-run the workflow specifying the same basecaller model with the
                        --basecaller_cfg option.
                        ################################################################################
                        '''.stripIndent()
                        error input_data_err_msg
                    }
                }
        
        // Define basecaller config
        basecaller_cfg = metamap_basecaller_cfg
                | ifEmpty(params.basecaller_cfg)
                | first  // unpack from list

        if (run_tumor_only){
            // Import the table of ClairS-TO models, retaining:
            // 1. basecaller model name
            // 2. ClairS-TO model name
            lookup_table_cls = Channel
                                .fromPath("${projectDir}/data/clairs_models_to.tsv", checkIfExists: true)
                                | splitCsv(sep: '\t', header: true)
                                | map{ it -> [it.basecall_model_name, it.clairs_to_model_name, "-", it.clairs_to_nomodel_reason] }
        } else {
            lookup_table_cls = Channel
                .fromPath("${projectDir}/data/clairs_models.tsv", checkIfExists: true)
                | splitCsv(sep: '\t', header: true)
                | map{ it -> [it.basecall_model_name, it.clairs_model_name, it.liquid_model_name_override, it.clairs_nomodel_reason] }
        }
        // Check that the provided basecaller_cfg is in the table.
        // If the user asks liquid_tumor, and the column has a valid model, then use that.
        // Otherwise, use regular ClairS model.
        clairs_model_ch = basecaller_cfg
            | cross(lookup_table_cls)
        clairs_model = clairs_model_ch
            | filter{ caller, info -> info[1] != '-' }
            | map{
                caller, info -> 
                model = params.liquid_tumor && info[2] != '-' && params.bam_normal ? info[2] : info[1]
            }
        // Log which models have been chosen.
        // If the model is not supported, throw error with informative message.
        clairs_model_ch
        | map {
            caller, info -> 
            clairs_model = params.liquid_tumor && info[2] != '-' ? info[2] : info[1]
            clairs_nomodel_msg = info[3]
            [clairs_model, clairs_nomodel_msg]
        }
        | subscribe{
            model, clairs_nomodel_reason -> 
            if (clairs_nomodel_reason != '-'){
                throw new Exception(colors.red + "${clairs_nomodel_reason}" + colors.reset)
            }
        }


        // Import the table of Clair3 models, retaining:
        // 1. basecaller model name
        // 2. Clair3 model name
        lookup_table_cl3 = Channel
            .fromPath("${projectDir}/data/clair3_models.tsv", checkIfExists: true)
            | splitCsv(sep: '\t', header: true)
            | map{ it -> [it.basecall_model_name, it.clair3_model_name, it.clair3_nomodel_reason] }
        // Check that the provided basecaller_cfg is in the table.
        clair3_model = basecaller_cfg
            | cross(lookup_table_cl3)
            | filter{ caller, info -> info[1] != '-' }
            | map{caller, info -> info[1] }

        if (run_tumor_only){
            clair_vcf = snv_to(
                pass_bam_channel,
                bed,
                ref_channel,
                clairs_model,
                clair3_model
            )
        } else {
            clair_vcf = snv(
                pass_bam_channel,
                bed,
                ref_channel,
                clairs_model,
                clair3_model
            )
        }
        
        // Publish outputs in the appropriate folder
        clair_vcf.outputs | publish_snv
        snv_joint_report = clair_vcf.report_snv
    } else {
        snv_joint_report = Channel.empty()
    }
    
    // Start SV calling workflow
    if (params.sv){
        sv_result = sv(pass_bam_channel, ref_channel, OPTIONAL)
        sv_joint_report = sv_result.report_sv
    } else {
        sv_joint_report = Channel.empty()
    }

    // Extract modified bases
    if (params.mod){
        modkit_output = mod(
            pass_bam_channel,
            ref_channel,
            chromosome_codes
        )
        mod_joint_report = modkit_output.report_mod
    } else {
        mod_joint_report = Channel.empty()
    }

    // Collect all the reports
    for_joint_report
        .mix(snv_joint_report)
        .mix(sv_joint_report)
        .mix(mod_joint_report)
        // Remove the null values for the missing reports
        // Create a nested tuple with all the reports in it
        .groupTuple(by:0)
        // Add additional data types
        .combine(versions)
        .combine(parameters) | report

    // Emit version and parameters
    publish(versions.concat(parameters).concat(report.out))
}

workflow.onComplete {
    Pinguscript.ping_complete(nextflow, workflow, params)
}
workflow.onError {
    Pinguscript.ping_error(nextflow, workflow, params)
}
