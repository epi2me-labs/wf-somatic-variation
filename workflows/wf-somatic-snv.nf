
include {
    getParams;
    getVersions;
    vcfStats;
    makeReport;
    output_snv;
    clairs_select_het_snps;
    clairs_phase;
    clairs_haplotag;
    wf_build_regions;
    clairs_extract_candidates;
    clairs_create_paired_tensors;
    clairs_predict_pileup;
    clairs_merge_pileup;
    clairs_create_fullalignment_paired_tensors;
    clairs_predict_full;
    clairs_merge_full;
    clairs_full_hap_filter;  
    clairs_merge_final;
    clairs_create_paired_tensors_indels;
    clairs_predict_pileup_indel;
    clairs_merge_pileup_indels;
    clairs_create_fullalignment_paired_tensors_indels;
    clairs_predict_full_indels;
    clairs_merge_full_indels;
    clairs_merge_final_indels;
    clairs_merge_snv_and_indels;
    add_missing_vars as add_missing_snvs;
    add_missing_vars as add_missing_indels;
    getVariantType as getSNVs;
    getVariantType as getIndels;
    change_count
} from "../modules/local/wf-somatic-snv.nf"

include {
    make_chunks;
    pileup_variants;
    aggregate_pileup_variants;
    select_het_snps;
    phase_contig;
    get_qual_filter;
    create_candidates;
    evaluate_candidates;
    aggregate_full_align_variants;
    merge_pileup_and_full_vars;
    aggregate_all_variants;
} from "../modules/local/wf-clair3-snp.nf"

include {annotate_vcf as annotate_snv} from '../modules/local/common.nf'

// workflow module
workflow snv {
    take:
        bam_channel
        bed
        ref
        clairs_model
        clair3_model
    main:
        // Define the presets for fast/normal mode once 
        def clair3_mode = [:]
        if (params.fast_mode){
            clair3_mode = [min_snps_af:"0.15", min_indels_af: "0.2", min_cvg: "8"]
        } else {
            clair3_mode = [min_snps_af:"0.08", min_indels_af: "0.15", min_cvg: "4"]
        } 

        // Prepare the a channel for hybrid or genotyping mode.
        typing_ch = Channel.of("$projectDir/data/OPTIONAL_FILE").map {
            vcf = params.hybrid_mode_vcf ?: params.genotyping_mode_vcf ?: it
            mode = params.hybrid_mode_vcf ? "--hybrid_mode_vcf_fn" : params.genotyping_mode_vcf ? "--genotyping_mode_vcf_fn": null
            [file(vcf, checkifExists: true), mode]
        }.collect()  // Use collect to create a value channel

        // Branch tumor and normal for downstream works
        bam_channel.branch{
            tumor: it[2].type == 'tumor'
            normal: it[2].type == 'normal'
        }.set{forked_channel}        

        // Initialize contigs and intervals for each pair of tumor/normal bam files
        forked_channel.normal
            .map{ bam, bai, meta -> [ meta.sample, bam, bai, meta ] } 
            .cross(
                forked_channel.tumor.map{ bam, bai, meta -> [ meta.sample, bam, bai, meta ] }
            )
            .map { normal, tumor ->
                    [normal[1], normal[2], tumor[1], tumor[2], tumor[3]]
                } 
            .map{ it -> it.flatten() }
            .set{paired_samples}
        wf_build_regions( paired_samples, ref.collect(), clairs_model, bed, typing_ch )
        
        // Define default and placeholder channels for downstream processes.
        bam_for_germline = params.germline ? bam_channel : Channel.empty()
        paired_vcfs = Channel.empty()
        aggregated_vcfs = bam_channel.map { bam, bai, meta -> [meta, file("$projectDir/data/OPTIONAL_FILE"), file("$projectDir/data/OPTIONAL_FILE")] }
        forked_vcfs = Channel.empty()
        
        // If a normal vcf is provided, skip clair3 calling for the normal sample.
        if (params.normal_vcf){
            // Use the metadata from the bam channel
            // Use the vcf for the normal channel.
            normal_vcf = bam_channel
                .filter{bam, bai, meta -> meta.type == 'normal'}
                .map{ bam, bai, meta -> [meta, file("${params.normal_vcf}"), file("${params.normal_vcf}.tbi")] }
            // Prepare the channel for Clair3 to contain only normal samples.
            bam_for_germline = bam_channel
                .filter{bam, bai, meta -> meta.type != 'normal'}
        } 

        // Extract contigs to feed into make_chunks to keep consistent parameters
        clair3_input_ctgs = wf_build_regions.out.contigs_file.map() { it -> [it[0].sample, it[1]] }

        /* =============================================== */
        /* Run Clair3 functions on each bam independently. */
        /* =============================================== */

        // Prepare the bam channel for the chunking.
        bam_for_germline.map{bam, bai, meta -> [meta.sample, bam, bai, meta]}
            .combine(clair3_input_ctgs, by: 0)
            .map{
                sample, bam, bai, meta, ctgs -> [meta, bam, bai, ctgs]
            }
            .combine(ref)
            .combine(bed)
            .set{bams}

        // Prepare the chunks for each bam file.
        make_chunks(bams, clair3_mode)
        chunks = make_chunks.out.chunks_file
            .splitText(){ 
                cols = (it[1] =~ /(.+)\s(.+)\s(.+)/)[0]
                [it[0], ["contig": cols[1], "chunk_id":cols[2], "total_chunks":cols[3]]]
                } 
        contigs = make_chunks.out.contigs_file.splitText() { it -> [it[0], it[1].trim()] }

        // Run the "pileup" caller on all chunks and collate results
        // > Step 1 
        bam_for_germline
            .combine(ref)
            .combine(bed)
            .map{bam, bai, meta, ref, fai, ref_cache, bed ->
                [meta, bam, bai, ref, fai, ref_cache, bed]
            }
            .combine(chunks, by:0)
            .combine(clair3_model)
            .set{fragments}
        pileup_variants(fragments, clair3_mode)

        // Aggregate outputs
        // Clairs model is required to define the correct var_pct_phasing 
        // value (0.7 for r9, 0.8 for r10).
        pileup_variants.out.pileup_vcf_chunks
            .groupTuple(by: 0)
            .combine(ref)
            .combine(
                make_chunks.out.contigs_file, by: 0
            )
            .combine(
                clair3_model
            ) .set{pileup_vcfs}
        aggregate_pileup_variants(pileup_vcfs)

        // Filter collated results to produce per-contig SNPs for phasing.
        // > Step 2
        aggregate_pileup_variants.out.pileup_vcf
            .combine(aggregate_pileup_variants.out.phase_qual, by: 0)
            .combine(contigs, by: 0)
            .set{ aggregated_pileup_vcf }
        select_het_snps(aggregated_pileup_vcf)
        // Perform phasing for each contig.
        // Combine the het variants with the input bam channels 
        // using the metadata as joining criteria (by:2), and then add 
        // the reference channel.
        // Then run the phasing
        phase_inputs = select_het_snps.out.het_snps_vcf
            .combine(bam_for_germline, by: 2)
            .combine(ref)
        phase_contig(phase_inputs)
        phase_contig.out.phased_bam_and_vcf.set { phased_bam_and_vcf }

        // Find quality filter to select variants for "full alignment"
        // processing, then generate bed files containing the candidates.
        // > Step 5
        get_qual_filter(aggregate_pileup_variants.out.pileup_vcf)
        aggregate_pileup_variants.out.pileup_vcf
            .combine(ref)
            .combine(get_qual_filter.out.full_qual, by: 0)
            .combine(contigs, by: 0)
            .set{pileup_and_qual}
        create_candidates(pileup_and_qual)

        // Run the "full alignment" network on candidates. Have to go through a
        // bit of a song and dance here to generate our input channels here
        // with various things duplicated (again because of limitations on 
        // `each` and tuples).
        // > Step 6
        create_candidates.out.candidate_bed.flatMap {
            x ->
                // output globs can return a list or single item
                y = x[2]; if(! (y instanceof java.util.ArrayList)){y = [y]}
                // effectively duplicate chr for all beds - [chr, bed]
                y.collect { [[x[0], x[1]], it] } }
                .set{candidate_beds}
        // produce something emitting: [[chr, bam, bai, vcf], [chr20, bed], [ref, fai, cache], model]
        bams_beds_and_stuff = phased_bam_and_vcf
            .map{meta, ctg, bam, bai, vcf, tbi -> [ [meta, ctg], bam, bai, vcf, tbi ]}
            .cross(candidate_beds)
            .combine(ref.map {it->[it]})
            .combine(clair3_model)
        // take the above and destructure it for easy reading
        bams_beds_and_stuff.multiMap {
            it ->
                bams: it[0].flatten()
                candidates: it[1].flatten()
                ref: it[2]
                model: it[3]
            }.set { mangled }
        // phew! Run all-the-things
        evaluate_candidates(mangled.bams, mangled.candidates, mangled.ref, mangled.model, clair3_mode)

        // merge and sort all files for all chunks for all contigs
        // gvcf is optional, stuff an empty file in, so we have at least one
        // item to flatten/collect and tthis stage can run.
        evaluate_candidates.out.full_alignment.groupTuple(by: 0)
            .combine(ref)
            .combine(make_chunks.out.contigs_file, by:0)
            .set{to_aggregate}
        aggregate_full_align_variants(to_aggregate)

        // merge "pileup" and "full alignment" variants, per contig
        // note: we never create per-contig VCFs, so this process
        //       take the whole genome VCFs and the list of contigs
        //       to produce per-contig VCFs which are then finally
        //       merge to yield the whole genome results.

        // First merge whole-genome results from pileup and full_alignment
        //   for each contig...
        // note: the candidate beds aren't actually used by the program for ONT
        // > Step 7
        aggregate_pileup_variants.out.pileup_vcf
            .combine(aggregate_full_align_variants.out.full_aln_vcf, by: 0)
            .combine(ref)
            .combine( candidate_beds
                        .map { it->[it[0][0], it[1]] }
                        .groupTuple(by:0), by:0 )
            .combine(contigs, by:0)
            .set{to_aggregate_pileup_and_full}
        merge_pileup_and_full_vars(to_aggregate_pileup_and_full)

        // Finally, aggregate full variants for each sample
        merge_pileup_and_full_vars.out.merged_vcf
            .groupTuple(by:0)
            .combine(Channel.fromPath("$projectDir/data/OPTIONAL_FILE"))
            .combine(ref)
            .combine(make_chunks.out.contigs_file, by: 0)
            .set{final_vcfs}
        aggregate_all_variants( final_vcfs )

        // Before proceeding to ClairS, we need to prepare the appropriate 
        // matched VCFs for tumor/normal pairs.
        // First, we branch based on whether they are tumor or normal:
        // If skip phasing, set channel for downstream compatibility.
        if (params.germline){
            aggregated_vcfs = aggregate_all_variants.out.final_vcf
            aggregated_vcfs
                .branch{
                    tumor: it[0].type == 'tumor'
                    normal: it[0].type == 'normal'
                }.set{forked_vcfs}
        }

        // Then we can combine tumor and normal for the same sample.        
        // If normal VCF is provided, then use it; if not check if germline calling is on, and eventually skip it
        normal_vcf_to_cross = params.normal_vcf ? normal_vcf : params.germline ? forked_vcfs.normal : Channel.empty()
        // Check if phasing is on, and if not use empty channel
        tumor_vcf_to_cross = params.germline ? forked_vcfs.tumor : Channel.empty()
        // Perform VCF pairing
        normal_vcf_to_cross
            .map{ meta, vcf, tbi -> [ meta.sample, vcf, tbi, meta ] } 
            .cross(
                tumor_vcf_to_cross.map{ meta, vcf, tbi -> [ meta.sample, vcf, tbi, meta ] }
            )
            .map { normal, tumor ->
                    [tumor[3], tumor[1], tumor[2], normal[3], normal[1], normal[2], ]
                } 
            .map{ it -> it.flatten() }
            .set{paired_vcfs}

        // If skip phasing, set channel for downstream compatibility.
        if (!params.germline){
            paired_vcfs = Channel.empty()
            aggregated_vcfs = bam_channel.map { bam, bai, meta -> [meta, file("$projectDir/data/OPTIONAL_FILE"), file("$projectDir/data/OPTIONAL_FILE")] }
        } else{
            aggregated_vcfs = aggregate_all_variants.out.final_vcf
        }
        /* ============================================ */
        /* Run ClairS functions from here on.           */
        /* The workflow will run partly in parallel     */
        /* with Clair3, and finishing off when Clair3   */
        /* completes the germline calling.              */
        /* ============================================ */

        /*
        /  Running the pileup model.
        */
        // Combine the channels of each chunk, with the pair of
        // bams and the all the bed for each sample
        wf_build_regions.out.chunks_file
                .splitText(){
                        cols = (it[1] =~ /(.+)\s(.+)\s(.+)/)[0]
                        region_map = ["contig": cols[1], "chunk_id":cols[2], "total_chunks":cols[3]]
                        [it[0], region_map]
                    }
                .combine(
                    paired_samples
                        .map{
                            ctrbm, ctrbi, canbm, canbi, meta -> [meta, ctrbm, ctrbi, canbm, canbi]
                            }, by: 0
                    )
                .combine(wf_build_regions.out.split_beds, by: 0)
                .combine(ref)
                .combine(clairs_model)
                .set{chunks}
        clairs_contigs = wf_build_regions.out.contigs_file.splitText() { it -> [it[0], it[1].trim()] }

        // Extract candidates for the tensor generation. Now with support of hybrid/genotyping mode
        clairs_extract_candidates(chunks, typing_ch)

        // Create all candidates channel for downstream analyses.
        // CW-1949: consider also the candidate for hybrid genotyping.
        // First extract the candidates for SNVs, Indels and hybrid typing.
        candidate_regions = clairs_extract_candidates.out.candidates_snvs.transpose().map{
            meta, region, list, target -> [meta, target]
            }.mix(
                clairs_extract_candidates.out.candidates_indels.transpose().map{
                    meta, region, list, target -> [meta, target]
                    }
            ).mix(
                clairs_extract_candidates.out.candidates_hybrids.transpose().map{
                    meta, region, target -> [meta, target]
                    }
        // Then, the add the lists of candidate files.
            ).mix(
                clairs_extract_candidates.out.candidates_snvs.transpose().map{
                    meta, region, list, target -> [meta, list]
                    }
            ).mix(
                clairs_extract_candidates.out.candidates_indels.transpose().map{
                    meta, region, list, target -> [meta, list]
                    }
            ).unique().groupTuple(by: 0)
        // Collect the bed separately as they need to be in "candidates/bed/" subdir 
        candidates_beds = clairs_extract_candidates.out.candidates_bed.map{
            meta, region, bed -> [meta, bed]
            }.groupTuple(by: 0)

        // Prepare the paired tensors for each tumor/normal pair.
        clairs_create_paired_tensors(chunks.combine(clairs_extract_candidates.out.candidates_snvs.transpose(), by: [0,1]))
        
        // Predict variants based on the paired pileup model
        clairs_predict_pileup(clairs_create_paired_tensors.out)

        // Combine all predicted pileup vcf into one pileup.vcf file.
        clairs_predict_pileup.out
            .groupTuple(by:0)
            .combine(wf_build_regions.out.contigs_file, by: 0)
            .set{collected_vcfs}
        clairs_merge_pileup(collected_vcfs, ref.collect())

        /*
        /  Processing the full alignments
        */
        // Create compatible channels in case no phasing is required
        forked_channel.tumor.map{it -> [it[2], it[0], it[1]]}
                    .combine(clairs_contigs, by: 0)
                    .map{meta, bam, bai, contigs -> [meta.sample, contigs, bam, bai, meta]}
                    .combine(
                        forked_channel.normal.map{it -> [it[2].sample, it[0], it[1]]}, by: 0
                        )
                    .combine(
                        clairs_extract_candidates.out.candidates_snvs.transpose().map{it -> [it[0].sample, it[1].contig, it[3]]}, by: [0,1] )
                    .combine(ref)
                    .combine(clairs_model)
                    .set{paired_phased_channel}

        // Create tagged reads channel for downstream analyses.
        // Set contigs to all to avoid issues of duplicated file names.
        tagged = bam_channel
            .map{bam, bai, meta -> [meta.sample, 'all', bam, bai, meta]}

        // Otherwise, perform the phasing and tagging
        if (params.germline){
            // Extract the germline heterozygote sites using both normal and tumor
            // VCF files
            clairs_select_het_snps(paired_vcfs.combine(clairs_contigs, by: 0))
            // Collect tumor het sites.
            het_tumor = clairs_select_het_snps.out.tumor_hets
                .combine(forked_channel.tumor
                            .map {bam, bai, meta -> [meta, bam, bai]}, by: 0
                            )
                .combine(ref)
            // If normal data are phased, use them. Otherwise, empty channel.
            het_normal = params.phase_normal ? clairs_select_het_snps.out.normal_hets : Channel.empty()
            
            // Combine heterozygote sites.
            het_normal
                .combine(forked_channel.normal
                            .map {bam, bai, meta -> [meta, bam, bai]}, by: 0
                            )
                .combine(ref)
                .mix(het_tumor)
                .set { het_to_phase }

            // Phase and haplotag the selected vcf and bams (tumor-only or both).
            het_to_phase | clairs_phase | clairs_haplotag

            // Create tagged reads channel for downstream analyses
            tagged = clairs_haplotag.out.phased_data

            // Branch to separate tumor and normal tagged bams
            tagged
                .branch{
                    tumor: it[4].type == 'tumor'
                    normal: it[4].type == 'normal'
                }.set{f_phased_channel}
            tagged_bam_tumor = f_phased_channel.tumor
            // If phase normal is specified then use the haplotagged bams, otherwise use the original ones.
            tagged_bam_normal = params.phase_normal ? f_phased_channel.normal : forked_channel.normal.map{it -> [it[2].sample, it[0], it[1], it[2]]}

            // Prepare the channel for the tensor generation.
            tagged_bam_tumor
                .combine(tagged_bam_normal, by: params.phase_normal ? [0, 1] : 0)
                .map{sample, contig, tbam, tbai, tmeta, nbam, nbai, nmeta -> 
                        [sample, contig, tbam, tbai, tmeta, nbam, nbai]
                }
                .combine(
                    clairs_extract_candidates.out.candidates_snvs.transpose().map{it -> [it[0].sample, it[1].contig, it[3]]}, by: [0,1] )
                .combine(ref)
                .combine(clairs_model)
                .set{ paired_phased_channel }
        }

        // Create the full-alignment tensors
        clairs_create_fullalignment_paired_tensors(paired_phased_channel)

        // Prediction of variant on full-alignments
        clairs_predict_full(clairs_create_fullalignment_paired_tensors.out.full_tensors)

        // Combine the full-alignment somatic SNV VCF
        clairs_predict_full.out.full_vcfs
            .groupTuple(by:0)
            .combine(wf_build_regions.out.contigs_file, by: 0)
            .combine(ref)
            .set{collected_full_vcfs}
        clairs_merge_full(collected_full_vcfs)

        // Almost there!!!
        //
        // Perform the haplotype-based filtering. This step can either be:
        //  1. By contig: very fast and resource efficient, but does not provide 
        //     the same results as ClairS.
        //  2. Full vcf: very slow and less resource efficient, but provide 
        //     results identical to ClairS.
        //
        // If split haplotype filter is specified, run by contig:
        if (params.skip_haplotype_filter){
            clairs_merge_pileup.out.pileup_vcf
                .combine(
                    clairs_merge_full.out.full_vcf, by:0
                )
                .combine( ref )
                .set{ clair_all_variants }
            snv_calling = clair_all_variants | clairs_merge_final
        } else {
            // Create channel with the tumor bam, all the VCFs
            // (germline, pileup and full-alignment) and the 
            // reference genome.
            tagged
                .map{ samp, ctg, bam, bai, meta -> [meta, bam, bai] }
                .groupTuple(by: 0)
                .combine(
                    aggregated_vcfs.map{it -> [it[0], it[1]]}, by: 0
                )
                .combine(
                    clairs_merge_pileup.out.pileup_vcf, by: 0
                )
                .combine(
                    clairs_merge_full.out.full_vcf, by:0
                )
                .combine( ref )
                .set{ clair_all_variants }
            // Apply the filtering and create the final VCF.
            snv_calling = clair_all_variants | clairs_full_hap_filter | clairs_merge_final
        }

        // Add missing genotypes if either hybrid_vcf or genotyping_vcf are provided
        if (params.hybrid_mode_vcf || params.genotyping_mode_vcf){
            snv_calling.pileup_vcf
                .combine(snv_calling.pileup_tbi, by: 0)
                .combine(candidate_regions, by: 0)
                .combine(candidates_beds, by: 0)
                .set{pretyping_snvs}
            snv_typed = add_missing_snvs(pretyping_snvs, typing_ch)
        } 

        // Multiple output channels are created here for convenience
        // in downstream analyses. pileup_vcf refers to a VCF file
        // that can undergo changes downstream, like combination with
        // the indels, whereas snv_vcf is emitted as raw output.
        pileup_vcf = params.hybrid_mode_vcf || params.genotyping_mode_vcf ? snv_typed.pileup_vcf : snv_calling.pileup_vcf
        pileup_tbi = params.hybrid_mode_vcf || params.genotyping_mode_vcf ? snv_typed.pileup_tbi : snv_calling.pileup_tbi
        snv_vcf = params.hybrid_mode_vcf || params.genotyping_mode_vcf ? snv_typed.pileup_vcf : snv_calling.pileup_vcf
        snv_tbi = params.hybrid_mode_vcf || params.genotyping_mode_vcf ? snv_typed.pileup_tbi : snv_calling.pileup_tbi
        
        // Perform indel calling if the model is appropriate
        if (params.basecaller_cfg.startsWith('dna_r10')){
            // Create paired tensors for the indels candidates
            clairs_create_paired_tensors_indels(chunks.combine(clairs_extract_candidates.out.candidates_indels.transpose(), by: [0,1]))

            // Create paired tensors for the indels candidates
            clairs_predict_pileup_indel( clairs_create_paired_tensors_indels.out )

            // Merge and sort all the pileup indels
            clairs_predict_pileup_indel.out
                .groupTuple(by:0)
                .combine(wf_build_regions.out.contigs_file, by: 0)
                .set{collected_indels}
            clairs_merge_pileup_indels(collected_indels, ref.collect())

            // Create the full alignment indel paired tensors
            // If skip phasing requested, then use the original bam files 
            // with the contig info...
            // Create default channel in case germline is deactivated
            forked_channel.tumor
                .map{bam, bai, meta -> [meta, bam, bai]}
                .combine(clairs_contigs, by:0)
                .map{meta, bam, bai, contig -> [meta.sample, contig, bam, bai, meta]}
                .combine(forked_channel.normal.map{it -> [it[2].sample, it[0], it[1]]}, by: 0)
                .combine(
                    clairs_extract_candidates.out.candidates_indels.transpose().map{it -> [it[0].sample, it[1].contig, it[3]]}, by: [0,1] )
                .combine(ref)
                .combine(clairs_model)
                .set{paired_phased_indels_channel}

            // Otherwise, use the phased data
            if (params.germline){
                // Check if the normal is requested, and if so consider that as well.
                tagged
                    .branch{
                        tumor: it[4].type == 'tumor'
                        normal: it[4].type == 'normal'
                    }.set{f_phased_channel}
                tagged_bam_tumor = f_phased_channel.tumor
                // If phase normal is specified then use the haplotagged bams, otherwise use the original ones.
                tagged_bam_normal = params.phase_normal ? f_phased_channel.normal : forked_channel.normal.map{it -> [it[2].sample, it[0], it[1], it[2]]}

                // Prepare the channel for the tensor generation.
                tagged_bam_tumor
                    .combine(tagged_bam_normal, by: params.phase_normal ? [0, 1] : 0)
                    .map{sample, contig, tbam, tbai, tmeta, nbam, nbai, nmeta -> 
                            [sample, contig, tbam, tbai, tmeta, nbam, nbai]
                    }
                    .combine(
                        clairs_extract_candidates.out.candidates_indels.transpose().map{it -> [it[0].sample, it[1].contig, it[3]]}, by: [0,1] )
                    .combine(ref)
                    .combine(clairs_model)
                    .set{ paired_phased_indels_channel }
            }

            // Run full alignment stage
            clairs_create_fullalignment_paired_tensors_indels( paired_phased_indels_channel )

            // Predict full alignment indels
            clairs_predict_full_indels(clairs_create_fullalignment_paired_tensors_indels.out.full_tensors)

            // Merge full alignment indels
            clairs_predict_full_indels.out.full_vcfs
                .groupTuple(by:0)
                .combine(wf_build_regions.out.contigs_file, by: 0)
                .combine(ref)
                .set{collected_full_vcfs}
            clairs_merge_full_indels( collected_full_vcfs )

            // Merge final indel file
            clairs_merge_pileup_indels.out.pileup_vcf
                .combine( clairs_merge_full_indels.out.full_vcf, by: 0 )
                .combine(ref)
                .set { merged_indels_vcf }
            called_indels = clairs_merge_final_indels(merged_indels_vcf)

            // Add missing genotypes if either hybrid_vcf or genotyping_vcf are provided
            if (params.hybrid_mode_vcf || params.genotyping_mode_vcf){
                called_indels.indel_vcf 
                    .combine(called_indels.indel_tbi, by: 0)
                    .combine(candidate_regions, by: 0)
                    .combine(candidates_beds, by: 0)
                    .set{pretyping_indels}
                typed_indels = add_missing_indels(pretyping_indels, typing_ch)

                // First, extract the correct variation types to avoid duplicated sites
                // in the joint VCF.
                getSNVs(pileup_vcf.combine(pileup_tbi, by: 0), 'snps')
                getIndels(typed_indels.pileup_vcf.combine(typed_indels.pileup_tbi, by: 0), 'indels')

                // Then, create the correct channel structure for the next process
                getSNVs.out.vcf
                    .join(getIndels.out.vcf, by: 0, remainder: true )
                    .map { it - null }
                    .map { [it[0], it[1..-1]] }
                    .set{ all_vcfs }
                // Repeat with the TBIs
                getSNVs.out.tbi
                    .join(getIndels.out.tbi, by: 0, remainder: true )
                    .map { it - null }
                    .map { [it[0], it[1..-1]] }
                    .set{ all_tbis }

            // If not required, use the base InDels VCF.
            } else {
                // First create a tuple of snv+indels VCFs, group them by metadata,
                // remove all null values, and set them as VCF file
                pileup_vcf
                    .join(called_indels.indel_vcf, by: 0, remainder: true )
                    .map { it - null }
                    .map { [it[0], it[1..-1]] }
                    .set{ all_vcfs }
                // Repeat with the TBIs
                pileup_tbi
                    .join(called_indels.indel_tbi, by: 0, remainder: true )
                    .map { it - null }
                    .map { [it[0], it[1..-1]] }
                    .set{ all_tbis }
            }
            // Create a combined channel of VCFs and TBIs.
            // If only one VCF is passed, then the concatenated vcf will match
            // the SNV VCF. Then, merge them.
            all_vcfs.combine(all_tbis, by: 0).set{snv_and_indels}

            // Add check for empty snv and indels channel.
            clairs_merge_snv_and_indels( snv_and_indels )

            // Prepare the new output VCF files
            indels_vcf = params.hybrid_mode_vcf || params.genotyping_mode_vcf ? typed_indels.pileup_vcf : called_indels.indel_vcf
            indels_tbi = params.hybrid_mode_vcf || params.genotyping_mode_vcf ? typed_indels.pileup_tbi : called_indels.indel_tbi
            pileup_vcf = clairs_merge_snv_and_indels.out.pileup_vcf
            pileup_tbi = clairs_merge_snv_and_indels.out.pileup_tbi
        }

        // Annotate the mutation type in the format XX[N>N]XX
        // where XX are the flanking regions of a given size 
        // For now, only K = 3 is provided.
        pileup_vcf
            .combine(pileup_tbi, by: 0)
            .combine(ref)
            .set{clairs_vcf}
        change_count(clairs_vcf)
        ch_vcf = change_count.out.mutype_vcf
        ch_tbi = change_count.out.mutype_tbi
        
        // Add snpEff annotation if requested
        if (params.annotation){
            annotate_snv(ch_vcf.combine(ch_tbi, by:0), 'somatic-snv')
            ch_vcf = annotate_snv.out.annot_vcf.map{meta, vcf, tbi -> [meta, vcf]}
            ch_tbi = annotate_snv.out.annot_vcf.map{meta, vcf, tbi -> [meta, tbi]}
            clinvar_vcf = annotate_snv.out.annot_vcf_clinvar
            gene_txt = annotate_snv.out.gene_txt
        // Otherwise, create optional file from the vcf channel to preserve the structure
        } else {
            clinvar_vcf = ch_vcf.map{meta, vcf -> [meta, file("$projectDir/data/OPTIONAL_FILE")]}
        }

        // Generate basic statistics for the VCF file
        vcfStats(ch_vcf.combine(ch_tbi, by: 0))

        // Create the report for the variants called
        software_versions = getVersions()
        workflow_params = getParams()
        ch_vcf
            .combine(ch_tbi, by: 0)
            .combine(vcfStats.out[0], by: 0)
            .combine(change_count.out.changes, by: 0)
            .combine(clinvar_vcf, by: 0)
            .combine(software_versions)
            .combine(workflow_params)
            .combine(typing_ch)
            .set{ reporting }
        makeReport(reporting)

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
                    })
            .concat(
                vcfStats.out.map{
                    meta, stats -> [stats, "${meta.sample}/snv/varstats"]
                    })
            .concat(
                change_count.out.changes.map{
                    meta, spectra -> [spectra, "${meta.sample}/snv/change_counts"]
                    })
            .concat(
                workflow_params.map{
                    params -> [params, "info/snv/"]
                    })
            .concat(
                software_versions.map{
                    versions -> [versions, "info/snv/"]
                    })
            .concat(
                makeReport.out.html.map{
                    it -> [it, null]
                })
            .concat(
                snv_vcf.map{meta, vcf -> [vcf, "${meta.sample}/snv/vcf/"]}
                )
            .concat(
                snv_tbi.map{meta, tbi -> [tbi, "${meta.sample}/snv/vcf/"]}
                )
            .set{outputs}
        if (params.germline && !params.fast_mode){
            outputs
                .concat(
                    forked_vcfs.tumor.map{
                        meta, vcf, tbi -> [vcf, "${meta.sample}/snv/vcf/"]
                        })
                .concat(
                    forked_vcfs.tumor.map{
                        meta, vcf, tbi -> [tbi, "${meta.sample}/snv/vcf/"]
                        })
                .concat(
                    forked_vcfs.normal.map{
                        meta, vcf, tbi -> [vcf, "${meta.sample}/snv/vcf/"]
                        })
                .concat(
                    forked_vcfs.normal.map{
                        meta, vcf, tbi -> [tbi, "${meta.sample}/snv/vcf/"]
                        })
                .set{ outputs }
        }
        if (params.annotation){
            outputs
                .concat(
                    gene_txt.map{
                        meta, gene -> [gene, "${meta.sample}/snv/annot/"]
                        })
                .concat(
                    clinvar_vcf.map{
                        meta, vcf -> [vcf, "${meta.sample}/snv/annot/"]
                        })
                .set{outputs}
        }
        if (params.basecaller_cfg.startsWith('dna_r10')){
            outputs
                .concat(
                    indels_vcf.map{meta, vcf -> [vcf, "${meta.sample}/snv/vcf/"]}
                    )
                .concat(
                    indels_tbi.map{meta, tbi -> [tbi, "${meta.sample}/snv/vcf/"]}
                    )
                .set { outputs }
        }

    emit:
       outputs 
}
