import groovy.json.JsonBuilder

include {
    nanomonsv_parse;
    nanomonsv_filter;
    nanomonsv_classify;
    nanomonsv_get;
    annotate_classify;
    annotate_filter;
    sortVCF;
    getVersions;
    getParams;
    report;
    output_sv;
    bwa_index;
} from "../modules/local/wf-somatic-sv.nf"
include {bgzipper; tabixer} from "../modules/local/common.nf"


workflow somatic_sv {
    take:
        input_bam
        reference
        optional_file
    main:
        // Run nanomonsv parse on all the files
        nanomonsv_parse(input_bam, reference.collect())

        // Separate normal and tumor
        nanomonsv_parse.out.branch{
            tumor: it[3].type == 'tumor'
            normal: it[3].type == 'normal'
        }.set{forked_channel}

        // This section combine tumor/normal pairs for each sample.
        // Use cross to combine based on the sample ID, returning a tuple with this structure:
        // [ 
        //  [sampleID, cram, crai], 
        //  [sampleID, tissue, cram, crai] 
        // ] 
        // The first entry is always the normal, the second is the tumor tissue.
        // Then uses map to create a metadata for each sample x tissue combination and 
        // returns a tuple with [metadata, CRAM normal, CRAI normal, CRAM tumor, CRAI tumor] 
        forked_channel.normal
            .map{ bam, bai, cache, meta -> [ meta.sample, bam, bai, cache, meta ] } 
            .cross(
                forked_channel.tumor.map{ bam, bai, cache, meta -> [ meta.sample, bam, bai, cache, meta ] }
            )
            .map { normal, tumor ->
                    [tumor[4], normal[1], normal[2], normal[3], tumor[1], tumor[2], tumor[3]]
                } 
            .set{paired_samples}

        // Run standard nanomonsv get if 1 fragment is specified
        nanomonsv_get(paired_samples, reference)
        ch_vcf = nanomonsv_get.out.vcf
        ch_txt = nanomonsv_get.out.txt
        ch_sbd = nanomonsv_get.out.single_breakend

        // Filter SV calls, removing everything in tandem repeat
        // if an appropriate BED file is provided.
        if (params.tr_bed != null && file("${params.tr_bed}.tbi", checkIfExists: true)) {
            // Define input channel for both the bed.gz and bed.gz.tbi
            // Check if it is compressed. 
            if (params.tr_bed.endsWith('.gz')){
                tr_bed = Channel.fromPath(params.tr_bed, checkIfExists: true)
                // If it is check if it is indexed. If not, build it.
                if (file("${params.tr_bed}.tbi", checkIfExists: true)){
                    tr_bed_tbi = Channel.fromPath("${params.tr_bed}.tbi")
                } else {
                    tr_bed_tbi = tabixer(tr_bed)
                }
            } else {
                // If it is not compressed, compress and index it.
                Channel.fromPath(params.tr_bed) | bgzipper | tabixer
                tr_bed = bgzipper.out.bgzip
                tr_bed_tbi = bgzipper.out.tbi
            }
            // Filter sites tandem repeat elements
            nanomonsv_filter(ch_txt, ch_vcf, tr_bed, tr_bed_tbi) | annotate_filter
            ch_vcf = annotate_filter.out.vcf
            ch_txt = annotate_filter.out.txt
        }

        // If requested, perform insert classification
        if (params.classify_insert) {
            // Try to fetch the BWA-MEM indexes if they exists.
            // If not, regenerate the indexes.
            if (!file(params.ref + ".amb").exists() || 
                !file(params.ref + ".ann").exists() || 
                !file(params.ref + ".bwt").exists() || 
                !file(params.ref + ".pac").exists() || 
                !file(params.ref + ".sa").exists()){
                bwa_index(reference)
                indexed = bwa_index.out.bwa_ref
            } else {
                reference.map{
                    it -> [it[0], it[1], it[2], file(params.ref + ".amb"), file(params.ref + ".ann"), file(params.ref + ".bwt"), file(params.ref + ".pac"), file(params.ref + ".sa")]
                }.set{indexed}
            }
            // Check if the results have insert size of at least 100bp in size. If not, skip.
            n_valid_inserts = ch_txt.splitCsv(sep: '\t', header: true).filter{ meta, file -> file.Inserted_Seq.length() >= 100 }.count()
            // Run nanomonsv classify
            nanomonsv_classify( ch_txt, ch_vcf, indexed, n_valid_inserts ) | annotate_classify
            ch_txt = annotate_classify.out.txt
            ch_vcf = annotate_classify.out.vcf
        }

        // Sort and index each SV
        sortVCF(ch_vcf)

        // Prepare reports and outputs
        software_versions = getVersions()
        workflow_params = getParams()
        report(
            sortVCF.out.vcf_gz.collect(),
            sortVCF.out.vcf_tbi.collect(),
            optional_file,
            software_versions, 
            workflow_params)

        // Output everything
        // The report gets saved in the main outut directory (null out subfolder)
        // whereas the rest goes in the sv subdir.
        sortVCF.out.vcf_gz
            .map{ meta, vcf -> [ vcf, null ] }
            .concat(
                sortVCF.out.vcf_tbi.map{
                    meta, tbi -> [tbi, null]
                    })
            .concat(
                ch_txt.map{
                    meta, txt -> [txt, "${meta.sample}/sv/txt"]
                    })
            .concat(
                ch_sbd.map{
                    meta, sbd -> [sbd, "${meta.sample}/sv/single_breakend"]
                })
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
                sortVCF.out.vcf_gz,
                sortVCF.out.vcf_tbi,
                ch_txt,
                ch_sbd
            )
}
