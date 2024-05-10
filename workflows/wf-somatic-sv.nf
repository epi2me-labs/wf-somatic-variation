import groovy.json.JsonBuilder

include {
    severus;
    sortVCF;
    getVersions;
    getParams;
    report;
    output_sv;
} from "../modules/local/wf-somatic-sv.nf"
include {
    decompress_ref as decompress;
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

        // Create input channel of paired BAMs
        input_xam
            | map{xam, xai, meta -> [meta.sample, meta, xam, xai]}
            | branch{
                tumor: it[1].type == 'tumor'
                normal: it[1].type == 'normal'
            }
            | set{forked_xam}
        paired_xam = forked_xam.normal
        | combine(forked_xam.tumor, by: 0)
        | map{
            sample, meta_n, xam_n, xai_n, meta_t, xam_t, xai_t ->
            [meta_t, xam_n, xai_n, xam_t, xai_t]
        }
        // Run severus
        // Collect the reference to avoid queueing issues, and allowing future
        // implementation of multiple samples
        severus(
            paired_xam,
            reference.collect(),
            tr_bed
        )

        // Sort output
        sortVCF(severus.out.vcf)

        // Add snpEff annotation if requested
        if (params.annotation){
            vcf_to_annotate = sortVCF.out.vcf_gz
                .combine(sortVCF.out.vcf_tbi, by:0)
                .map{ it << '*' }
            annotate_sv(vcf_to_annotate, 'somatic-sv')
            ch_vcf = annotate_sv.out.annot_vcf.map{meta, vcf, tbi -> [meta, vcf]}
            ch_tbi = annotate_sv.out.annot_vcf.map{meta, vcf, tbi -> [meta, tbi]}
        // Otherwise, create optional file from the vcf channel to preserve the structure
        } else {
            ch_vcf = sortVCF.out.vcf_gz
            ch_tbi = sortVCF.out.vcf_tbi
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
                    it -> [it, null]
                })
            .set{outputs}
        outputs | output_sv
    emit:
        soma_sv = report.out.html.concat(
                ch_vcf,
                ch_tbi
            )
}
