import groovy.json.JsonBuilder

include {
    nanomonsv_parse;
    nanomonsv_filter;
    nanomonsv_classify;
    annotate_classify;
    nanomonsv_get;
    annotate_filter;
    postprocess_nanomon_vcf;
    sortVCF;
    getVersions;
    getParams;
    report;
    output_sv;
    bwa_index;
} from "../modules/local/wf-somatic-sv.nf"
include {bgzipper; tabixer; annotate_vcf as annotate_sv} from "../modules/local/common.nf"

workflow somatic_sv {
    take:
        input_xam
        reference
        optional_file
    main:
        // Create mismatched control panel channel
        // nanomonsv requires to provide the input in the format
        // root_dir/root_name
        // Where root_name is the common part of all the files in the
        // root_dir after stripping file extensions (which we do here
        // with simpleName). We can ignore the other files as only the
        // common path is required by nanomonsv.
        if (params.control_panel){
            control_panel = Channel
                .fromPath("${params.control_panel}/*", checkIfExists: true)
                .map{fname -> [file(params.control_panel), fname.simpleName]}
                .first()
        } else {
            control_panel = Channel
                .fromPath("$projectDir/data/OPTIONAL_FILE")
                .map{
                    [it, 'OPTIONAL']
                }
        }

        // Run nanomonsv parse on all the files
        nanomonsv_parse(input_xam, reference.collect())

        // Separate normal and tumor
        nanomonsv_parse.out.branch{
            tumor: it[0].type == 'tumor'
            normal: it[0].type == 'normal'
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
            .map{ 
                meta, xam, xam_idx, cache -> 
                [ meta.sample, meta, xam, xam_idx, cache ]
            }
            .cross(
                forked_channel.tumor.map{ 
                    meta, xam, xam_idx, cache -> 
                    [ meta.sample, meta, xam, xam_idx, cache ]
                }
            )
            .map { normal, tumor ->
                [tumor[1]] + normal[2..4] + tumor[2..4]
            } 
            .set{parsed_samples}

        // Run standard nanomonsv get if 1 fragment is specified
        nanomonsv_get(parsed_samples, reference, control_panel)
        ch_vcf = nanomonsv_get.out.vcf
        ch_txt = nanomonsv_get.out.txt
        ch_sbd = nanomonsv_get.out.single_breakend
        read_lists = nanomonsv_get.out.read_lists

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
                    fasta, fai, cache, ref_path -> 
                    [fasta, fai, cache, ref_path, file(params.ref + ".amb"), file(params.ref + ".ann"), file(params.ref + ".bwt"), file(params.ref + ".pac"), file(params.ref + ".sa")]
                }.set{indexed}
            }
            // Check if the results have insert size of at least 100bp in size. If not, skip.
            n_valid_inserts = ch_txt.splitCsv(sep: '\t', header: true).filter{ meta, file -> file.Inserted_Seq.length() >= 100 }.count()
            // Check that the genome type is of an appropriate type
            ch_txt.combine(ch_vcf, by:0).branch{
                proper: it[0].genome_build=='hg19' || it[0].genome_build=='hg38'
                improper: true
            }.set{branched_svs}

            // Log if a file can't be annotated due to wrong build.
            // This will emit a pre-annotation VCF, without collapsing and predisposing for 
            // future multisample with some passing and potentially some failing.
            branched_svs.improper.subscribe{
                log.warn "Genome build incompatible with insert_classify for ${it[0].sample}. The final vcf will not contain additional information on potential transposable and repetitive elements."
            }

            // Run nanomonsv classify
            nanomonsv_classify(branched_svs.proper , indexed, n_valid_inserts ) | annotate_classify

            // Define the right outputs
            post_annot_vcf = annotate_classify.out.annotated.mix(branched_svs.improper)

            ch_vcf = post_annot_vcf.map{meta, txt, vcf -> [meta, vcf]}
            ch_txt = post_annot_vcf.map{meta, txt, vcf -> [meta, txt]}
        }

        // Post-process VCF file
        proc_vcf = postprocess_nanomon_vcf(ch_vcf)

        // Sort and index each SV
        sortVCF(proc_vcf)

        // Add snpEff annotation if requested
        if (params.annotation){
            annotate_sv(sortVCF.out.vcf_gz.combine(sortVCF.out.vcf_tbi, by:0), 'somatic-sv')
            proc_vcf = annotate_sv.out.annot_vcf.map{meta, vcf, tbi -> [meta, vcf]}
            proc_tbi = annotate_sv.out.annot_vcf.map{meta, vcf, tbi -> [meta, tbi]}
            gene_txt = annotate_sv.out.gene_txt
        // Otherwise, create optional file from the vcf channel to preserve the structure
        } else {
            proc_vcf = sortVCF.out.vcf_gz
            proc_tbi = sortVCF.out.vcf_tbi
        }

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
        // This saves the final VCF+TBI after post-processing (single-sample, sorted and annotated)
        proc_vcf
            .map{ meta, vcf -> [ vcf, null ] }
            .concat(
                proc_tbi.map{
                    meta, tbi -> [tbi, null]
                    })
            // This saves the raw VCF files from nanomonsv
            .concat(
                ch_vcf.map{
                    meta, txt -> [txt, "${meta.sample}/sv/vcf"]
                    })
            // This saves the tables from nanomonsv for dual-breakend SVs, single-breakend SVs
            // and also the reads supporting each SV event
            .concat(
                ch_txt.map{
                    meta, txt -> [txt, "${meta.sample}/sv/txt"]
                    })
            .concat(
                read_lists.map{
                    meta, txt -> [txt, "${meta.sample}/sv/txt"]
                    })
            .concat(
                ch_sbd.map{
                    meta, sbd -> [sbd, "${meta.sample}/sv/single_breakend"]
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
        
        if (params.annotation){
            outputs
                .concat(
                    gene_txt.map{
                        meta, gene -> [gene, "${meta.sample}/sv/annot/"]
                        })
                .set{outputs}
        }
        outputs | output_sv
    emit:
        soma_sv = report.out.html.concat(
                sortVCF.out.vcf_gz,
                sortVCF.out.vcf_tbi,
                ch_txt,
                ch_sbd
            )
}
