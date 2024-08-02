
include {
    getParams;
    getVersions;
    vcfStats;
    makeReport;
    publish_snv;
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
    concat_hap_filtered_vcf;  
    clairs_merge_final;
    add_missing_vars;
    getVariantType;
    clairs_merge_snv_and_indels;
    change_count;
    concat_bams
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

include {
    annotate_vcf as annotate_snv;
    concat_vcfs as concat_snp_vcfs;
    sift_clinvar_vcf as sift_clinvar_snp_vcf
} from '../modules/local/common.nf'

// workflow module
workflow snv {
    take:
        bam_channel
        bed
        ref
        clairs_model
        clair3_model
    main:
        // Log the chosen models
        clair3_model.subscribe{ log.info " - Clair3 model: ${it}" }
        clairs_model.subscribe{ log.info " - ClairS model: ${it}" }

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
        forked_channel = bam_channel
        | branch{
            tumor: it[2].type == 'tumor'
            normal: it[2].type == 'normal'
        }

        // Initialize contigs and intervals for each pair of tumor/normal bam files
        paired_samples = forked_channel.normal
            | map{ bam, bai, meta -> [ meta.sample, bam, bai, meta ] } 
            | cross(
                forked_channel.tumor.map{ bam, bai, meta -> [ meta.sample, bam, bai, meta ] }
            )
            // Extract a tuple of [normal_bam, normal_bai, tumor_bam, tumor_bai, tumor_meta]
            | map { normal, tumor -> [normal[1], normal[2], tumor[1], tumor[2], tumor[3]] } 
            | map{ it -> it.flatten() }
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
        bams = bam_for_germline.map{bam, bai, meta -> [meta.sample, bam, bai, meta]}
            | combine(clair3_input_ctgs, by: 0)
            | map{
                sample, bam, bai, meta, ctgs -> [meta, bam, bai, ctgs]
            }
            | combine(ref)
            | combine(bed)
            | combine(clair3_model)

        // Prepare the chunks for each bam file.
        make_chunks(bams, clair3_mode)
        chunks = make_chunks.out.chunks_file
            .splitText(){ 
                cols = (it[1] =~ /(.+)\s(.+)\s(.+)/)[0]
                [it[0], ["contig": cols[1], "chunk_id":cols[2], "total_chunks":cols[3]]]
                } 
        contigs = make_chunks.out.contigs_file.splitText() { it -> [it[0], it[1].trim()] }
        cmd_file = make_chunks.out.cmd_file
        split_beds = wf_build_regions.out.split_beds

        // Run the "pileup" caller on all chunks and collate results
        // > Step 1 
        fragments = bam_for_germline
            | combine(ref)
            | combine(bed)
            | combine(split_beds)
            | map{bam, bai, meta, ref, fai, ref_cache, ref_path, bed, t_meta, bed_dir ->
                [meta, bam, bai, ref, fai, ref_cache, ref_path, bed, bed_dir]
            }
            | combine(chunks, by:0)
            | combine(clair3_model)
            | combine(cmd_file, by:0)
        pileup_variants(fragments, clair3_mode)

        // Aggregate outputs
        // Clairs model is required to define the correct var_pct_phasing 
        // value (0.7 for r9, 0.8 for r10).
        pileup_vcfs = pileup_variants.out.pileup_vcf_chunks
            | groupTuple(by: 0)
            | combine(ref)
            | combine(make_chunks.out.contigs_file, by: 0)
            | combine(clair3_model) 
            | combine(cmd_file, by: 0)
        aggregate_pileup_variants(pileup_vcfs)

        // Filter collated results to produce per-contig SNPs for phasing.
        // > Step 2
        aggregated_pileup_vcf = aggregate_pileup_variants.out.pileup_vcf
            | combine(aggregate_pileup_variants.out.phase_qual, by: 0)
            | combine(contigs, by: 0)
        select_het_snps(aggregated_pileup_vcf)
        // Perform phasing for each contig.
        // Combine the het variants with the input bam channels 
        // using the metadata as joining criteria (by:2), and then add 
        // the reference channel.
        // Then run the phasing
        phase_inputs = select_het_snps.out.het_snps_vcf
            .combine(bam_for_germline, by: 2)
            .combine(ref)
        phased_bam_and_vcf = phase_contig(phase_inputs)

        // Find quality filter to select variants for "full alignment"
        // processing, then generate bed files containing the candidates.
        // > Step 5
        get_qual_filter(aggregate_pileup_variants.out.pileup_vcf)
        pileup_and_qual = aggregate_pileup_variants.out.pileup_vcf
            | combine(ref)
            | combine(get_qual_filter.out.full_qual, by: 0)
            | combine(contigs, by: 0)
        create_candidates(pileup_and_qual)

        // Run the "full alignment" network on candidates. Have to go through a
        // bit of a song and dance here to generate our input channels here
        // with various things duplicated (again because of limitations on 
        // `each` and tuples).
        // > Step 6
        // This creates a nested tuple with the structure [ [meta, contig], candidate ]
        candidate_beds = create_candidates.out.candidate_bed.flatMap {
            x ->
                // output globs can return a list or single item
                y = x[2]; if(! (y instanceof java.util.ArrayList)){y = [y]}
                // effectively duplicate chr for all beds - [chr, bed]
                y.collect { [[x[0], x[1]], it] }
        }
        // produce something emitting: [[chr, bam, bai, vcf], [chr20, bed], [ref, fai, cache, ref_path], model]
        bams_beds_and_stuff = phased_bam_and_vcf
            | combine(cmd_file, by: 0)
            | map{meta, ctg, bam, bai, vcf, tbi, cmd -> [ [meta, ctg], bam, bai, vcf, tbi, cmd ]}
            | cross(candidate_beds)
            | combine(ref.map {it->[it]})
            | combine(clair3_model)
        // take the above and destructure it for easy reading
        mangled = bams_beds_and_stuff.multiMap {
            it ->
                bams: it[0].flatten()
                candidates: it[1].flatten()
                ref: it[2]
                model: it[3]
            }
        // phew! Run all-the-things
        evaluate_candidates(mangled.bams, mangled.candidates, mangled.ref, mangled.model, clair3_mode)

        // merge and sort all files for all chunks for all contigs
        // gvcf is optional, stuff an empty file in, so we have at least one
        // item to flatten/collect and tthis stage can run.
        to_aggregate = evaluate_candidates.out.full_alignment
            | groupTuple(by: 0)
            | join(
                pileup_variants.out.pileup_gvcf_chunks.groupTuple(by: 0),
                remainder: true,
                by:0
            )
            // If no GVCF was requested, it will be null and replaced with OPTIONAL_FILE
            | map{
                meta, fa_vcfs, gvcfs ->
                n_gvcfs = gvcfs == null ? [file("${projectDir}/data/OPTIONAL_FILE")] : gvcfs
                [meta, fa_vcfs, n_gvcfs]
            }
            | combine(ref)
            | combine(make_chunks.out.contigs_file, by:0)
            | combine(cmd_file, by:0)
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
        to_aggregate_pileup_and_full = aggregate_pileup_variants.out.pileup_vcf
            | combine(aggregate_full_align_variants.out.full_aln_vcf, by: 0)
            | join(aggregate_full_align_variants.out.non_var_gvcf, by: 0, remainder: true)
            // If no GVCF was requested, it will be null and replaced with OPTIONAL_FILE
            | map{meta, p_vcf, p_tbi, f_vcf, f_tbi, gvcf -> 
                n_gvcf = gvcf == null ? file("${projectDir}/data/OPTIONAL_FILE") : gvcf
                [meta, p_vcf, p_tbi, f_vcf, f_tbi, n_gvcf]
                }
            | combine(ref)
            | combine(
                candidate_beds
                | map { it->[it[0][0], it[1]] }  // Extract meta and the candidate 
                | groupTuple(by:0), by:0 )
            | combine(contigs, by:0)
        merge_pileup_and_full_vars(to_aggregate_pileup_and_full)

        // Prepare a GVCF channel for the final combination
        if (params.GVCF){
            merged_gvcfs = merge_pileup_and_full_vars.out.merged_gvcf
            | map{meta, contig, gvcf -> [meta, gvcf]}
            | groupTuple(by: 0)
        // If no GVCF requested, then use metadata from vcf and OPTIONAL_FILE
        } else {
            merged_gvcfs = merge_pileup_and_full_vars.out.merged_vcf
            | groupTuple(by: 0)
            | map{[it[0], file("$projectDir/data/OPTIONAL_FILE")]}

        }

        // Finally, aggregate full variants for each sample
        final_vcfs = merge_pileup_and_full_vars.out.merged_vcf
            | groupTuple(by:0)
            | combine(merged_gvcfs, by: 0)
            | combine(ref)
            | combine(make_chunks.out.contigs_file, by: 0)
            | combine(cmd_file, by: 0)
        aggregate_all_variants( final_vcfs )

        // Before proceeding to ClairS, we need to prepare the appropriate 
        // matched VCFs for tumor/normal pairs.
        // First, we branch based on whether they are tumor or normal:
        // If skip phasing, set channel for downstream compatibility.
        if (params.germline){
            aggregated_vcfs = aggregate_all_variants.out.final_vcf
            // Branch tumor and normal sample by checking the type in meta (it[0])
            forked_vcfs = aggregated_vcfs
            | branch{
                tumor: it[0].type == 'tumor'
                normal: it[0].type == 'normal'
            }
        }

        // Then we can combine tumor and normal for the same sample.        
        // If normal VCF is provided, then use it; if not check if germline calling is on, and eventually skip it
        normal_vcf_to_cross = params.normal_vcf ? normal_vcf : params.germline ? forked_vcfs.normal : Channel.empty()
        // Check if phasing is on, and if not use empty channel
        tumor_vcf_to_cross = params.germline ? forked_vcfs.tumor : Channel.empty()
        // Perform VCF pairing
        paired_vcfs = normal_vcf_to_cross
        | map{ meta, vcf, tbi -> [ meta.sample, vcf, tbi, meta ] } 
        | cross(
            tumor_vcf_to_cross.map{ meta, vcf, tbi -> [ meta.sample, vcf, tbi, meta ] }
        )
        // Extract a tuple with structure
        // [ tumor_meta, tumor_vcf, tumor_tbi, normal_meta, normal_vcf, normal_tbi ]
        | map { normal, tumor ->
            [tumor[3], tumor[1], tumor[2], normal[3], normal[1], normal[2], ]
        } 
        | map{ it -> it.flatten() }

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
        chunks = wf_build_regions.out.chunks_file
        // Split file in three columns as any character combination (`(.+)`)
        // splitting any whitespace character (`\s`).
        | splitText(){
            cols = (it[1] =~ /(.+)\s(.+)\s(.+)/)[0]
            region_map = ["contig": cols[1], "chunk_id":cols[2], "total_chunks":cols[3]]
            [it[0], region_map]
        }
        | combine(
            paired_samples
            | map{
                ctrbm, ctrbi, canbm, canbi, meta -> [meta, ctrbm, ctrbi, canbm, canbi]
            }, by: 0
        )
        | combine(split_beds, by: 0)
        | combine(ref)
        | combine(clairs_model)
        clairs_contigs = wf_build_regions.out.contigs_file.splitText() { it -> [it[0], it[1].trim()] }

        // Extract candidates for the tensor generation. Now with support of hybrid/genotyping mode
        extr_candidates = clairs_extract_candidates(chunks, typing_ch)
        clairs_cand_snvs = extr_candidates.candidates_snvs
        clairs_cand_indels = extr_candidates.candidates_indels
        clairs_cand_hyb = extr_candidates.candidates_hybrids
        clairs_cand_bed = extr_candidates.candidates_bed

        // Create all candidates channel for downstream analyses.
        // CW-1949: consider also the candidate for hybrid genotyping.
        // First extract the candidates for SNVs, Indels and hybrid typing.
        candidate_regions = clairs_cand_snvs
            | transpose
            | map{
                meta, region, var_type, list, target -> [meta, target]
                }
            | mix(
                clairs_cand_indels
                    | transpose
                    | map{
                        meta, region, var_type, list, target -> [meta, target]
                    }
            )
            | mix(
                clairs_cand_hyb
                | transpose
                | map{
                    meta, region, target -> [meta, target]
                }
            )
            // Then, the add the lists of candidate files.
            | mix(
                clairs_cand_snvs
                    | transpose
                    | map{
                        meta, region, var_type, list, target -> [meta, list]
                    }
            )
            | mix(
                clairs_cand_indels.transpose().map{
                    meta, region, var_type, list, target -> [meta, list]
                    }
            )
            | unique
            | groupTuple(by: 0)
        // Collect the bed separately as they need to be in "candidates/bed/" subdir 
        candidates_beds = clairs_cand_bed
            | map{
                meta, region, bed -> [meta, bed]
            }
            | groupTuple(by: 0)

        // Run the pileup step on each tumor/normal pair of candidates for SNV
        // and Indels, if the model is r10.
        // The workflow will create the paired tensors, run the prediction, and
        // then merge the called VCFs.
        merged_pileup_vcfs = chunks
            | combine(
                clairs_cand_snvs
                    | mix(clairs_cand_indels)
                    | transpose
                    // If it is R9, ignore indel type
                    | filter{
                        meta, region, var_type, candidates, intervals ->
                        meta.basecall_models[0].startsWith('dna_r10') || var_type == 'snv'
                    },
                by: [0,1]
            )
            | clairs_create_paired_tensors
            | clairs_predict_pileup
            | groupTuple(by:[0,1])
            | combine(wf_build_regions.out.contigs_file, by: 0)
            | combine(ref)
            | clairs_merge_pileup

        /*
        /  Processing the full alignments
        */
        // Create compatible channels in case no phasing is required
        // Otherwise, perform the phasing and tagging
        if (params.germline){
            // Extract the germline heterozygote sites using both normal and tumor
            // VCF files
            het_sites = paired_vcfs | combine(clairs_contigs, by: 0) | clairs_select_het_snps
            // Collect tumor het sites.
            het_tumor = het_sites.tumor_hets
                | combine(forked_channel.tumor
                            | map {bam, bai, meta -> [meta, bam, bai]}, by: 0
                            )
                | combine(ref)
            // If normal data are phased, use them. Otherwise, empty channel.
            het_normal = params.phase_normal ? het_sites.normal_hets : Channel.empty()
            
            // Combine heterozygote sites.
            het_to_phase  = het_normal
                | combine(
                    forked_channel.normal
                    | map {bam, bai, meta -> [meta, bam, bai]}, by: 0
                )
                | combine( ref )
                | mix( het_tumor )

            // Phase and haplotag the selected vcf and bams (tumor-only or both).
            tagged = het_to_phase | clairs_phase | clairs_haplotag

            // Branch to separate tumor and normal tagged bams
            f_phased_channel = tagged
                | branch{
                    tumor: it[4].type == 'tumor'
                    normal: it[4].type == 'normal'
                }
            tagged_bam_tumor = f_phased_channel.tumor
            // If phase normal is specified then use the haplotagged bams, otherwise use the original ones.
            tagged_bam_normal = params.phase_normal ? f_phased_channel.normal : forked_channel.normal.map{it -> [it[2].sample, it[0], it[1], it[2]]}

            // Prepare the channel for the tensor generation.
            paired_phased_channel = tagged_bam_tumor
            // If phased normal is required, then combine the data by sample ID and sequence.
            // Otherwise, combine by sample ID only
            | combine(tagged_bam_normal, by: params.phase_normal ? [0, 1] : 0)
            | map{sample, contig, tbam, tbai, tmeta, nbam, nbai, nmeta -> 
                    [sample, contig, tbam, tbai, tmeta, nbam, nbai]
            }
            | combine(
                clairs_cand_snvs
                | mix(clairs_cand_indels)
                | transpose
                | map{
                    meta, region, var_type, candidate_list, candidate_file ->
                    [meta.sample, region.contig, var_type, candidate_file]}, by: [0,1]
            )
            // If it is R9, ignore indel type
            | filter{
                sample, contig, tbam, tbai, tmeta, nbam, nbai, var_type, intervals ->
                tmeta.basecall_models[0].startsWith('dna_r10') || var_type == 'snv'
            }
            | combine(ref)
            | combine(clairs_model)
        } else {
            // Create compatible channels in case no phasing is required
            paired_phased_channel = forked_channel.tumor.map{xam, xai, meta -> [meta, xam, xai]}
            | combine(clairs_contigs, by: 0)
            | map{meta, bam, bai, contigs -> [meta.sample, contigs, bam, bai, meta]}
            | combine(
                forked_channel.normal.map{xam, xai, meta -> [meta.sample, xam, xai]}, by: 0
                )
            | combine(
                clairs_cand_snvs
                | mix(clairs_cand_indels)
                | transpose
                | map{
                    meta, region, var_type, candidate_list, candidate_file ->
                    [meta.sample, region.contig, var_type, candidate_file]}, by: [0,1]
            )
            // If it is R9, ignore indel type
            | filter{
                sample, contig, tbam, tbai, tmeta, nbam, nbai, var_type, intervals ->
                tmeta.basecall_models[0].startsWith('dna_r10') || var_type == 'snv'
            }
            | combine(ref)
            | combine(clairs_model)

            // Create tagged reads channel for downstream analyses.
            // Set contigs to all to avoid issues of duplicated file names.
            tagged = bam_channel
            | map{bam, bai, meta -> [meta.sample, 'all', bam, bai, meta]}
        }

        // Create the full-alignment tensors
        merged_full_vcfs = clairs_create_fullalignment_paired_tensors(paired_phased_channel)
        | clairs_predict_full
        | groupTuple(by: [0,1])  // Group by meta and variant type.
        | combine(wf_build_regions.out.contigs_file, by: 0)
        | combine(ref)
        | clairs_merge_full

        // Perform the haplotype-based filtering for the SNVs, unless skip requested.
        if (params.skip_haplotype_filter){
            combined_vcf = merged_pileup_vcfs
            | combine(
                merged_full_vcfs, by: [0, 1]  // combine by both sample name and variant type
            )
            | combine( ref )
            | clairs_merge_final
        } else {
            // Create channel with the tumor bam, all the VCFs
            // (germline, pileup and full-alignment) and the 
            // reference genome.

            // When no germline is required, the workflow generates a single dummy
            // channel for the tagged bam, with 'all' in place of the single contigs. This
            // is necessary to avoid breaking other parts of the wf where multiple copies
            // of the BAM are passed to a single process, causing duplicated inputs with the
            // same name. Therefore, the issue is fixed here by defining whether to use the
            // tagged BAM, or the raw input to which we add the processed contigs.
            if (params.germline){
                xam_for_hap_filter = tagged
                | map{ samp, ctg, xam, xai, meta -> [meta, ctg, xam, xai] }
            } else {
                xam_for_hap_filter = bam_channel
                    | map{xam, xai, meta -> [meta, xam, xai]}
                    | combine(clairs_contigs, by:0)
                    | map{meta, xam, xai, ctg -> [meta, ctg, xam, xai]}
            }

            clair_all_variants = xam_for_hap_filter
            | combine(
                aggregated_vcfs, by: 0
            )
            | combine(
                merged_pileup_vcfs
                | combine(merged_full_vcfs, by:[0,1]),
                by:0
            )
            | combine( ref )

            // Apply the filtering and create the final SNV VCF.
            // 1. Run haplotype_filter on each contig individually
            // 2. Group the outputs by  metadata and variant type
            // 3. Add the contig list file and the reference channel
            // 4. Concatenate the files for both pileup and full-alignment
            // 5. Merge the filtered pileup and full alignment VCFs. 
            // clair_all_variants
            combined_vcf =  clair_all_variants
            | clairs_full_hap_filter
            | groupTuple(by: [0, 1])   // combine by both sample name and variant type
            | combine(wf_build_regions.out.contigs_file, by: 0)
            | combine(ref)
            | concat_hap_filtered_vcf
            | clairs_merge_final
        }

        // Add missing genotypes if either hybrid_vcf or genotyping_vcf are provided
        if (params.hybrid_mode_vcf || params.genotyping_mode_vcf){
            pretyping_snvs = combined_vcf
            | combine(candidate_regions, by: 0)
            | combine(candidates_beds, by: 0)
            combined_vcf = add_missing_vars(pretyping_snvs, typing_ch) | getVariantType
        }
        
        // Add check for empty snv and indels channel.
        snv_and_indels_vcf = combined_vcf
        | map { it - null }
        | map{ meta, var_type, vcf, tbi -> [meta, vcf, tbi] }
        | groupTuple(by:0)
        | clairs_merge_snv_and_indels

        // Concatenate haplotagged bams, if phasing requested
        if (params.germline){
            extensions = Channel.of(['cram', 'crai'])
            tagged_bams = concat_bams(
                tagged
                | map{samp, ctg, bam, bai, meta -> [meta, bam, bai]}
                | groupTuple(by:0)
                | combine(ref)
                | combine(extensions)
            )
        } else {
            tagged_bams = Channel.empty()
        }

        // Annotate the mutation type in the format XX[N>N]XX
        // where XX are the flanking regions of a given size 
        // For now, only K = 3 is provided.
        snv_and_indels_vcf
        | combine(ref)
        | change_count
        ch_vcf = change_count.out.mutype_vcf
        ch_tbi = change_count.out.mutype_tbi
        
        // Add snpEff annotation if requested
        if (params.annotation){
            // snpeff is slow so we'll just pass the whole VCF but annotate per contig
            annotations = annotate_snv(
                ch_vcf.combine(ch_tbi, by:0).combine(clairs_contigs, by:0), 'somatic-snv'
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
                    meta, html -> [html, null]
                })
            .concat(
                combined_vcf
                    | map{meta, var_type, vcf, tbi -> [meta, [vcf, tbi]] }
                    | transpose
                    | map{meta, file -> [file, "${meta.sample}/snv/vcf/"]}
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
                .concat(
                    tagged_bams.map{meta, bam, bai -> [bam, null]}
                    )
                .concat(
                    tagged_bams.map{meta, bam, bai -> [bai, null]}
                    )
                .set{ outputs }
        }
        if (params.annotation){
            outputs
                .concat(
                    clinvar_vcf.map{
                        meta, vcf -> [vcf, "${meta.sample}/snv/annot/"]
                        })
                .set{outputs}
        }
        // Emit GVCF if requested
        if (params.GVCF){
            outputs
                .concat(
                    aggregate_all_variants.out.final_gvcf.map{meta, gvcf, tbi -> [gvcf, "${meta.sample}/snv/gvcf/"]}
                    )
                .concat(
                    aggregate_all_variants.out.final_gvcf.map{meta, gvcf, tbi -> [tbi, "${meta.sample}/snv/gvcf/"]}
                    )
                .set { outputs }
        }

    emit:
       outputs = outputs
       snv_stats = vcfStats.out[0].combine(change_count.out.changes, by: 0)
       report_snv = makeReport.out.html
}
