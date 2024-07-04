// ClairS-TO specific processes
include {
    getParams;
    getVersions;
    wf_build_regions;
    clairs_to_extract_candidates;
    clairs_to_create_affirmative_model_tensors
    clairs_to_create_negational_model_tensors
    clairs_to_predict_pileup;
    clairs_to_pileup;
    clairs_to_merge_pileup;
    clairs_to_tag_non_somatic_db;
    clairs_to_merge_tagged;
    clairs_to_hap_filter;
    clairs_to_merge_hapfilt;
    clairs_to_postprocess;
    clairs_to_select_het_snps;
    clairs_to_phase;
    clairs_to_haplotag;
} from "../modules/local/wf-somatic-snv-to.nf"

// Processes shared with ClairS
include {
    clairs_merge_snv_and_indels;
    vcfStats;
    output_snv;
    change_count;
    makeReport;
} from "../modules/local/wf-somatic-snv.nf"

// Additional processes
include {
    annotate_vcf as annotate_snv;
    concat_vcfs as concat_snp_vcfs;
    sift_clinvar_vcf as sift_clinvar_snp_vcf
} from '../modules/local/common.nf'

// workflow module
workflow snv {
    take:
        xam_channel
        bed
        ref
        clairs_model
        clair3_model
    main:
        // Log chosen model
        clairs_model.subscribe{ log.info " - ClairS-TO model: ${it}" }

        // Prepare xam tuple to have the meta in front
        xam = xam_channel.map{xam, xam_idx, meta -> [meta, xam, xam_idx]}

        // Prepare a channel for VCF to be used in hybrid or genotyping mode.
        // The channel will ship both the appropriate ClairS-TO option and the VCF
        // as follow:
        // [ known_sites.vcf[.gz], "--[hybrid/genotyping]_mode_vcf_fn" ]
        typing_ch = Channel.fromPath("$projectDir/data/OPTIONAL_FILE").map {
            vcf = params.hybrid_mode_vcf ?: params.genotyping_mode_vcf ?: it
            mode = params.hybrid_mode_vcf ? "--hybrid_mode_vcf_fn" : params.genotyping_mode_vcf ? "--genotyping_mode_vcf_fn": null
            [file(vcf, checkifExists: true), mode]
        }.collect()  // Use collect to create a value channel

        // Run run_clairs_to in dryrun mode to generate minimal files required.
        wf_build_regions(xam, ref.collect(), clairs_model, bed, typing_ch )

        // Prepare channels and files to pass around
        chunks_file = wf_build_regions.out.chunks_file
        cmd_file = wf_build_regions.out.cmd_file
        split_bed_files = wf_build_regions.out.split_beds
        split_indel_bed_files = wf_build_regions.out.split_indel_beds
        contig_file = wf_build_regions.out.contigs_file
        contigs_channel = contig_file.splitText() { it -> [it[0], it[1].trim()] }
        
        // Create candidates for each chunk.
        // First, create the input channel.
        chunks = chunks_file
            | splitText(){
                    cols = (it[1] =~ /(.+)\s(.+)\s(.+)/)[0]
                    region_map = ["contig": cols[1], "chunk_id":cols[2], "total_chunks":cols[3]]
                    [it[0], region_map]
                }
            | combine(xam, by: 0)
            | combine(split_bed_files, by: 0)
            | combine(split_indel_bed_files, by: 0)
            | combine(ref)
            | combine(clairs_model)
        // Then, extract the candidates.
        clairs_to_extract_candidates(chunks, typing_ch)

        // Create affirmative and negational model tensors.
        // The processing for SNV and indel is identical, so mix them up,
        // then transpose them to have individual files as inputs.
        for_tensor_generation = clairs_to_extract_candidates.out.candidates_snvs
            | mix(clairs_to_extract_candidates.out.candidates_indels)
            | transpose
        // Generate affirmative and negational tensors.
        affirmative = for_tensor_generation | clairs_to_create_affirmative_model_tensors
        negational =  for_tensor_generation | clairs_to_create_negational_model_tensors

        // Combine the affirmative and negational tensors.
        tensors = affirmative
            | map{
                meta, var_type, region, candidate, model, affirm -> 
                [meta, var_type, region, candidate.getName(), model, candidate, affirm]
            }
            | combine(
                negational
                | map{
                    meta, var_type, region, candidate, model, negat -> 
                    [meta, var_type, region, candidate.getName(), model, candidate, negat]
                },
                by: [0, 1, 2, 3, 4]
            )
            | map{
                meta, var_type, region, cand_name, model, candidate_aff, affirm, candidate_neg, negat -> 
                [meta, var_type, region, candidate_aff, model, affirm, negat] 
            }

        // Perform prediction and pileup for each chunk.
        pileups = tensors
            | clairs_to_predict_pileup
            | combine(ref)
            | clairs_to_pileup
        
        // Merge pileup files for SNV and Indels separately
        // by grouping the tuple by meta and variant type.
        pileup_raw_vcf = pileups
            | groupTuple(by:[0, 1])
            | combine(contig_file, by: 0)
            | combine(ref)
            | clairs_to_merge_pileup
        
        // Haplotag reads.
        // First, extract the heterozygote sites.
        xam_tagged = pileup_raw_vcf
            | filter{meta, var_type, vcf -> var_type == 'snv'}
            | combine( contigs_channel, by: 0 )
            | combine(xam, by: 0)
            | combine(ref)
            | clairs_to_select_het_snps
            | clairs_to_phase
            | clairs_to_haplotag

        // Tag non-somatic sites using the database VCF, then merge
        // tagged SNVs and Indels separately.
        pileup_raw_vcf
            | combine( contigs_channel, by: 0 )
            | combine( ref )
            | clairs_to_tag_non_somatic_db

        tagged_vcf = clairs_to_tag_non_somatic_db.out 
            | groupTuple(by:[0, 1])
            | combine(contig_file, by: 0)
            | combine( ref )
            | clairs_to_merge_tagged

        // Perform haplotype filtering, post-processing and
        // final bgzip-indexing.
        // First, we run the haplotype filtering.
        hap_filt_chr = xam_tagged
            | groupTuple(by: 0)
            | combine(
                pileup_raw_vcf
                    | combine(tagged_vcf, by: [0, 1]),
                by: 0
            )
            | combine( contigs_channel, by: 0 )
            | combine(ref)
            | clairs_to_hap_filter
        // We then combine the hap. filtered variants.
        hap_filt = hap_filt_chr
            | groupTuple(by: [0, 1])
            | combine( ref )
            | combine( contig_file, by: 0 )
            | clairs_to_merge_hapfilt

        // We then apply the ClairS-TO post processing.
        pileup_vcf = hap_filt
            | combine(ref)
            | combine(cmd_file, by: 0)
            | clairs_to_postprocess

        // Finally, we concatenate SNVs and Indels in a single VCF file.
        final_vcf = pileup_vcf 
            | groupTuple(by: 0) 
            | clairs_merge_snv_and_indels

        // Annotate the mutation type in the format XX[N>N]XX
        // where XX are the flanking regions of a given size 
        // For now, only K = 3 is provided.
        final_vcf.pileup_vcf
            | combine(final_vcf.pileup_tbi, by:0)
            | combine(ref)
            | change_count
        ch_vcf = change_count.out.mutype_vcf
        ch_tbi = change_count.out.mutype_tbi
        
        // Add snpEff annotation if requested
        if (params.annotation){
            annotations = annotate_snv(
                ch_vcf.combine(ch_tbi, by:0).combine(contigs_channel, by:0), 'somatic-snv'
            )
            final_vcf = concat_snp_vcfs(annotations.groupTuple(by: 0), "${params.sample_name}.wf-somatic-snv").final_vcf
            clinvar_vcf = sift_clinvar_snp_vcf(final_vcf, "snv").final_vcf_clinvar
            ch_vcf = final_vcf.map{meta, vcf, tbi -> [meta, vcf]}
            ch_tbi = final_vcf.map{meta, vcf, tbi -> [meta, tbi]}
        // Otherwise, create optional file from the vcf channel to preserve the structure
        } else {
            clinvar_vcf = ch_vcf.map{meta, vcf -> [meta, file("$projectDir/data/OPTIONAL_FILE")]}
        }

        // Generate basic statistics for the VCF file
        ch_vcf | combine(ch_tbi, by: 0) | vcfStats

        // Create the report for the variants called
        software_versions = getVersions()
        workflow_params = getParams()
        ch_vcf
            | combine(ch_tbi, by: 0)
            | combine(vcfStats.out[0], by: 0)
            | combine(change_count.out.changes, by: 0)
            | combine(clinvar_vcf, by: 0)
            | combine(software_versions)
            | combine(workflow_params)
            | combine(typing_ch)
            | makeReport

        // Create a single channel with all the outputs
        // Send the output to the specified sub-directory of params.out_dir.
        // If null is passed, send it to out_dir/ directly.
        ch_vcf
            .map{ 
                meta, vcf -> [ vcf, null ]
            }
            .concat(
                ch_tbi.map{
                    meta, tbi -> [tbi, null]
                }
            )
            .concat(
                vcfStats.out.map{
                    meta, stats -> [stats, "${meta.sample}/snv/varstats"]
                }
            )
            .concat(
                change_count.out.changes.map{
                    meta, spectra -> [spectra, "${meta.sample}/snv/change_counts"]
                }
            )
            .concat(
                workflow_params.map{
                    params -> [params, "info/snv/"]
                }
            )
            .concat(
                software_versions.map{
                    versions -> [versions, "info/snv/"]
                }
            )
            .concat(
                makeReport.out.html.map{
                    meta, html -> [html, null]
                }
            )
            .concat(
                final_vcf.map{meta, vcf, tbi -> [vcf, "${meta.sample}/snv/vcf/"]}
                )
            .concat(
                final_vcf.map{meta, vcf, tbi -> [tbi, "${meta.sample}/snv/vcf/"]}
            )
            .set{outputs}
        if (params.annotation){
            outputs
                .concat(
                    clinvar_vcf.map{
                        meta, vcf -> [vcf, "${meta.sample}/snv/annot/"]
                    }
                )
                .set{outputs}
        }

    emit:
       outputs = outputs
       snv_stats = vcfStats.out[0].combine(change_count.out.changes, by: 0)
       report_snv = makeReport.out.html

}
