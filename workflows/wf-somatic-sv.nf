import groovy.json.JsonBuilder

include {
    severus;
    finaliseVCF;
    getVersions;
    getParams;
    report;
    publish_sv;
} from "../modules/local/wf-somatic-sv.nf"
include {
    decompress;
    bgzipper;
    tabixer;
    annotate_vcf as annotate_sv;
} from "../modules/local/common.nf"

workflow somatic_sv {
    take:
        input_xam
        reference
        optional_file
    main:
        // Filter SV calls, removing everything in tandem repeat
        // if an appropriate BED file is provided.
        if (params.tr_bed) {
            if (params.tr_bed.endsWith('.gz')){
                tr_bed = Channel.fromPath(params.tr_bed, checkIfExists: true) | decompress
            } else {
                tr_bed = Channel.fromPath(params.tr_bed, checkIfExists: true)
            }
        } else {
            tr_bed = Channel.fromPath("$projectDir/data/OPTIONAL_FILE")
        }

        // User provided PoN file
        if (params.pon_file) {
            pon_file = Channel.fromPath(params.pon_file, checkIfExists: true)
        } else {
            pon_file = Channel.fromPath("$projectDir/data/OPTIONAL_FILE")
        }
        boolean tumor_only = !params.bam_normal 
        // Create input channel of paired BAMs
        input_xam
            | map{xam, xai, meta -> [meta.sample, meta, xam, xai]}
            | branch{
                tumor: it[1].type == 'tumor'
                normal: it[1].type == 'normal'
            }
            | set{forked_xam}

        // split out the two bam streams for readability
        def normals = forked_xam.normal
        def tumors  = forked_xam.tumor

        if (tumor_only) {
            paired_xam = tumors.map{
                _sample, meta_t, xam_t, xai_t-> 
                [meta_t, optional_file, optional_file, xam_t, xai_t]
            }
        } else {
            paired_xam = normals
            | combine(tumors, by: 0)
            | map{
                _sample, _meta_n, xam_n, xai_n, meta_t, xam_t, xai_t ->
                [meta_t, xam_n, xai_n, xam_t, xai_t]
            }
        }
        // Run severus
        // Collect the reference to avoid queueing issues, and allowing future
        // implementation of multiple samples
        severus(
            paired_xam,
            reference.collect(),
            tr_bed,
            tumor_only,
            pon_file
        )

        // annotation bed files, always staged but only used if tumor only is true
        def seg_dup_annotation_files = channel.of([ 
                file("${workflow.projectDir}/data/hg38.segdups.bed.gz"),
                file("${workflow.projectDir}/data/hg38.segdups.bed.gz.tbi") 
        ])                

        // Sort and filter output
        // If tumor only, annotate with SegDup regions
        finaliseVCF(severus.out.vcf, tumor_only, seg_dup_annotation_files)

        // Add snpEff annotation if requested
        if (params.annotation){
            vcf_to_annotate = finaliseVCF.out.vcf_gz
                .combine(finaliseVCF.out.vcf_tbi, by:0)
                .map{ it << '*' }
            annotate_sv(vcf_to_annotate, 'somatic-sv')
            ch_vcf = annotate_sv.out.annot_vcf.map{meta, vcf, tbi -> [meta, vcf]}
            ch_tbi = annotate_sv.out.annot_vcf.map{meta, vcf, tbi -> [meta, tbi]}
        // Otherwise, create optional file from the vcf channel to preserve the structure
        } else {
            ch_vcf = finaliseVCF.out.vcf_gz
            ch_tbi = finaliseVCF.out.vcf_tbi
        }

        // Prepare reports and outputs
        software_versions = getVersions()
        workflow_params = getParams()
        report(
            ch_vcf,
            ch_tbi,
            optional_file,
            software_versions, 
            workflow_params)

        // Output everything
        // The report gets saved in the main outut directory (null out subfolder)
        // whereas the rest goes in the sv subdir.
        // This saves the final VCF+TBI after post-processing (single-sample, sorted and annotated)
        ch_vcf
            .map{ meta, vcf -> [ vcf, null ] }
            .concat(
                ch_tbi.map{
                    meta, tbi -> [tbi, null]
                    })
            // This saves the outputs from severus
            .concat(
                ch_vcf.map{
                    meta, txt -> [txt, "${meta.sample}/sv/vcf"]
                    })
            .concat(
                severus.out.all_outputs.map{
                    meta, outdir -> [outdir, "${meta.sample}/sv/"]
                    })
            // Add information
            .concat(
                workflow_params.map{
                    params -> [params, "info/sv/"]
                    })
            .concat(
                software_versions.map{
                    versions -> [versions, "info/sv/"]
                    })
            .concat(
                report.out.html.map{
                    meta, it -> [it, null]
                })
            .set{outputs}
        outputs | publish_sv
    emit:
        soma_sv = report.out.html.concat(
                ch_vcf,
                ch_tbi
            )
        sv_vcf = ch_vcf.combine(ch_tbi, by: 0)
        report_sv = report.out.html
}
