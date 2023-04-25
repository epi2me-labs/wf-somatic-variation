
process make_chunks {
    // Do some preliminaries. Ordinarily this would setup a working directory
    // that all other commands would make use off, but all we need here are the
    // list of contigs and chunks.
    label "wf_somatic_snv"
    cpus 1
    input:
        tuple val(meta), path(bam), path(bai), path(contigs), path(ref), path(fai), path(ref_cache), path(bed)
    output:
        tuple val(meta), path("clair_output_${meta.sample}_${meta.type}/tmp/CONTIGS"), emit: contigs_file
        tuple val(meta), path("clair_output_${meta.sample}_${meta.type}/tmp/CHUNK_LIST"), emit: chunks_file
        tuple val(meta), path("clair_output_${meta.sample}_${meta.type}/tmp/split_beds/"), emit: split_beds, optional: true
    script:
        def bedargs = params.bed ? "--bed_fn ${bed}" : ""
        def include_ctgs = params.include_all_ctgs ? "--include_all_ctgs" : ""
        """
        mkdir -p clair_output
        CTG_LIST=\$( tr '\\n' ',' < ${contigs} )
        echo "Running on contigs: \$CTG_LIST"
        python \$(which clair3.py) CheckEnvs ${include_ctgs} ${bedargs} \\
            --bam_fn ${bam} \\
            --output_fn_prefix clair_output_${meta.sample}_${meta.type} \\
            --ref_fn ${ref} \\
            --vcf_fn ${params.vcf_fn} \\
            --ctg_name \$CTG_LIST \\
            --chunk_num 0 \\
            --chunk_size 5000000 \\
            --threads 1  \\
            --qual 2 \\
            --sampleName ${meta.sample} \\
            --var_pct_full ${params.clair3_var_pct_full} \\
            --ref_pct_full ${params.clair3_ref_pct_full} \\
            --snp_min_af ${params.clair3_snp_min_af} \\
            --indel_min_af ${params.clair3_indel_min_af} \\
            --min_contig_size ${params.min_contig_size} 
        """
}


process pileup_variants {
    // Calls variants per region ("chunk") using pileup network.
    label "wf_somatic_snv"
    cpus 1
    errorStrategy 'retry'
    input:
        tuple val(meta), 
            path(bam), 
            path(bai), 
            path(ref), 
            path(fai), 
            path(ref_cache), 
            path(bed),
            val(region),
            val(model)
    output:
        // TODO: make this explicit, why is pileup VCF optional?
        tuple val(meta), 
            path("pileup_${meta.sample}_${meta.type}_${region.contig}_${region.chunk_id}.vcf"), 
            optional: true, 
            emit: pileup_vcf_chunks
        tuple val(meta), path("gvcf_tmp_path/*"), optional: true, emit: pileup_gvcf_chunks
    shell:
        // note: the VCF output here is required to use the contig
        //       name since that's parsed in the SortVcf step
        // note: snp_min_af and indel_min_af have an impact on performance
        // TODO Fix REF_PATH
        '''
        export REF_PATH=!{ref_cache}/%2s/%2s/%s
        python $(which clair3.py) CallVariantsFromCffi \\
            --chkpnt_fn ${CLAIR_MODELS_PATH}/clair3_models/!{model}/pileup \\
            --bam_fn !{bam} \\
            --call_fn pileup_!{meta.sample}_!{meta.type}_!{region.contig}_!{region.chunk_id}.vcf \\
            --ref_fn !{ref} \\
            --ctgName !{region.contig} \\
            --chunk_id !{region.chunk_id} \\
            --chunk_num !{region.total_chunks} \\
            --platform ont \\
            --fast_mode False \\
            --snp_min_af !{params.clair3_snp_min_af} \\
            --indel_min_af !{params.clair3_indel_min_af} \\
            --minMQ !{params.clair3_min_mq} \\
            --minCoverage !{params.clair3_min_coverage} \\
            --call_snp_only False \\
            --sampleName !{meta.sample} \\
            --vcf_fn !{params.vcf_fn} \\
            --enable_long_indel False \\
            --bed_fn \\
            --samtools samtools \\
            --gvcf !{params.GVCF} \\
            --temp_file_dir gvcf_tmp_path \\
            --pileup 
        ''' 
}


process aggregate_pileup_variants {
    // Aggregates and sorts all variants (across all chunks of all contigs)
    // from pileup network. Determines quality filter for selecting variants
    // to use for phasing.
    label "wf_somatic_snv"
    cpus 2
    input:
        tuple val(meta), 
            path(vcfs, stageAs: "input_vcfs/*"),
            path(ref), 
            path(fai), 
            path(ref_cache),
            path(contigs_file),
            val(model)
    output:
        tuple val(meta), path("pileup_${meta.sample}_${meta.type}.vcf.gz"), path("pileup_${meta.sample}_${meta.type}.vcf.gz.tbi"), emit: pileup_vcf
        tuple val(meta), path("phase_qual"), emit: phase_qual
    script:
        def var_pct_phasing = "--var_pct_phasing 0.7"
        if (model == "r941_prom_sup_g5014" || model == "r941_prom_hac_g5014") {
            var_pct_phasing = "--var_pct_phasing 0.8"
        }
        """
        pypy \$(which clair3.py) SortVcf \
            --input_dir input_vcfs/ \
            --vcf_fn_prefix pileup_${meta.sample}_${meta.type} \
            --output_fn pileup_${meta.sample}_${meta.type}.vcf \
            --sampleName ${params.sample_name} \
            --ref_fn ${ref} \
            --contigs_fn ${contigs_file}

        # Replaced bgzip with the faster bcftools index -n
        if [ "\$( bcftools index -n pileup_${meta.sample}_${meta.type}.vcf.gz )" -eq 0 ]; \
        then echo "[INFO] Exit in pileup variant calling"; exit 1; fi

        bgzip -@ ${task.cpus} -fdc pileup_${meta.sample}_${meta.type}.vcf.gz | \
            pypy \$(which clair3.py) SelectQual --phase ${var_pct_phasing} --output_fn .
        """
}


process select_het_snps {
    // Filters a VCF by contig, selecting only het SNPs.
    label "wf_somatic_snv"
    cpus 2
    input:
        tuple val(meta), path(pileup_vcf), path(pileup_tbi), path(split, stageAs: "phase_qual"), val(contig)
        // this is used implicitely by the program
        // https://github.com/HKU-BAL/Clair3/blob/329d09b39c12b6d8d9097aeb1fe9ec740b9334f6/preprocess/SelectHetSnp.py#L29
        
    output:
        tuple path("split_folder_${meta.sample}_${meta.type}/${contig}.vcf.gz"), path("split_folder_${meta.sample}_${meta.type}/${contig}.vcf.gz.tbi"), val(meta), val(contig), emit: het_snps_vcf
    shell:
        '''
        mkdir split_folder_!{meta.sample}_!{meta.type}/
        cp phase_qual split_folder_!{meta.sample}_!{meta.type}/
        pypy $(which clair3.py) SelectHetSnp \
            --vcf_fn !{pileup_vcf} \
            --split_folder split_folder_!{meta.sample}_!{meta.type}/ \
            --ctgName !{contig}

        bgzip -c split_folder_!{meta.sample}_!{meta.type}/!{contig}.vcf > split_folder_!{meta.sample}_!{meta.type}/!{contig}.vcf.gz
        tabix split_folder_!{meta.sample}_!{meta.type}/!{contig}.vcf.gz
        '''
}


process phase_contig {
    // Tags reads in an input BAM from heterozygous SNPs
    // The haplotag step was removed in clair-v0.1.11 so this step re-emits
    //   the original BAM and BAI as phased_bam for compatability,
    //   but adds the VCF as it is now tagged with phasing information
    //   used later in the full-alignment model
    label "wf_somatic_snv"
    cpus { params.use_longphase_intermediate ? 4 : 1 }
    input:
        tuple val(meta), 
            path(het_snps), 
            path(het_snps_tbi), 
            val(contig), 
            path(bam), 
            path(bai), 
            path(ref), 
            path(fai), 
            path(ref_cache)
    output:
        tuple val(meta), val(contig), path(bam), path(bai), path("phased_${meta.sample}_${meta.type}_${contig}.vcf.gz"), path("phased_${meta.sample}_${meta.type}_${contig}.vcf.gz.tbi"), emit: phased_bam_and_vcf
    script:
        if (params.use_longphase_intermediate)
        """
        export REF_PATH=${ref_cache}/%2s/%2s/%s
        echo "Using longphase for phasing"
        # longphase needs decompressed 
        bgzip -@ ${task.cpus} -dc ${het_snps} > snps.vcf
        longphase phase --ont -o phased_${meta.sample}_${meta.type}_${contig} \
            -s snps.vcf -b ${bam} -r ${ref} -t ${task.cpus}
        bgzip phased_${meta.sample}_${meta.type}_${contig}.vcf
        tabix -f -p vcf phased_${meta.sample}_${meta.type}_${contig}.vcf.gz
        """
        else
        """
        export REF_PATH=${ref_cache}/%2s/%2s/%s
        echo "Using whatshap for phasing"
        whatshap phase \
            --output phased_${meta.sample}_${meta.type}_${contig}.vcf.gz \
            --reference ${ref} \
            --chromosome ${contig} \
            --distrust-genotypes \
            --ignore-read-groups \
            ${het_snps} \
            ${bam}
        tabix -f -p vcf phased_${meta.sample}_${meta.type}_${contig}.vcf.gz
        """
}


process get_qual_filter {
    // Determines quality filter for selecting candidate variants for second
    // stage "full alignment" calling.
    label "wf_somatic_snv"
    cpus 2
    input:
        tuple val(meta), path("pileup.vcf.gz"), path("pileup.vcf.gz.tbi")
    output:
        tuple val(meta), path("output_${meta.sample}_${meta.type}/qual"), emit: full_qual
    shell:
        '''
        echo "[INFO] 5/7 Select candidates for full-alignment calling"
        mkdir output_!{meta.sample}_!{meta.type}
        bgzip -fdc pileup.vcf.gz | \
        pypy $(which clair3.py) SelectQual \
                --output_fn output_!{meta.sample}_!{meta.type} \
                --var_pct_full !{params.clair3_var_pct_full} \
                --ref_pct_full !{params.clair3_ref_pct_full} \
                --platform ont 
        '''
}


process create_candidates {
    // Create BED files for candidate variants for "full alignment" network
    // from the previous full "pileup" variants across all chunks of all chroms
    //
    // Performed per chromosome; output a list of bed files one for each chunk.
    label "wf_somatic_snv"
    cpus 2
    input:
        tuple val(meta), 
            path("pileup.vcf.gz"), 
            path("pileup.vcf.gz.tbi"), 
            path(ref), 
            path(fai), 
            path(ref_cache), 
            path("qual"), 
            val(contig)
        // this is used implicitely by the program
        // https://github.com/HKU-BAL/Clair3/blob/329d09b39c12b6d8d9097aeb1fe9ec740b9334f6/preprocess/SelectCandidates.py#L146
    output:
        tuple val(meta), val(contig), path("candidate_bed_${meta.sample}_${meta.type}/${contig}.*"), emit: candidate_bed, optional: true
    shell:
        // This creates BED files as candidate_bed/<ctg>.0_14 with candidates
        // along with a file the FULL_ALN_FILE_<ctg> listing all of the BED
        // files.  All we really want are the BEDs, the file of filenames is
        // used for the purposes of parallel in the original workflow.
        // https://github.com/HKU-BAL/Clair3/blob/329d09b39c12b6d8d9097aeb1fe9ec740b9334f6/scripts/clair3.sh#L218

        // TODO: would be nice to control the number of BEDs produced to enable
        // better parallelism.
        '''
        mkdir candidate_bed_!{meta.sample}_!{meta.type}/
        cp qual candidate_bed_!{meta.sample}_!{meta.type}/
        pypy $(which clair3.py) SelectCandidates \
            --pileup_vcf_fn pileup.vcf.gz \
            --split_folder candidate_bed_!{meta.sample}_!{meta.type} \
            --ref_fn !{ref} \
            --var_pct_full !{params.clair3_var_pct_full} \
            --ref_pct_full !{params.clair3_ref_pct_full} \
            --platform ont \
            --ctgName !{contig}
        '''
}


process evaluate_candidates {
    // Run "full alignment" network for variants in a candidate bed file.
    // phased_bam just references the input BAM as it no longer contains phase information.
    label "wf_somatic_snv"
    cpus 1
    errorStrategy 'retry'
    input:
        tuple val(meta), val(contig), path(phased_bam), path(phased_bam_index), path(phased_vcf), path(phased_tbi)
        tuple val(meta_2), val(contig_2), path(candidate_bed)
        tuple path(ref), path(fai), path(ref_cache)
        val model
    output:
        tuple val(meta), path("output_${meta.sample}_${meta.type}/full_alignment_*.vcf"), emit: full_alignment, optional: true
    script:
        filename = candidate_bed.name
        def ref_path = "${ref_cache}/%2s/%2s/%s:" + System.getenv("REF_PATH")
        """
        export REF_PATH=${ref_path}
        mkdir output_${meta.sample}_${meta.type}
        echo "[INFO] 6/7 Call low-quality variants using full-alignment model"
        python \$(which clair3.py) CallVariantsFromCffi \
            --chkpnt_fn \${CLAIR_MODELS_PATH}/clair3_models/${model}/full_alignment \
            --bam_fn ${phased_bam} \
            --call_fn output_${meta.sample}_${meta.type}/full_alignment_${filename}.vcf \
            --sampleName ${meta.sample} \
            --ref_fn ${ref} \
            --full_aln_regions ${candidate_bed} \
            --ctgName ${contig} \
            --add_indel_length \
            --minMQ ${params.clair3_min_mq} \
            --minCoverage ${params.clair3_min_coverage} \
            --platform ont \
            --use_gpu False \
            --samtools samtools \
            --enable_long_indel False \
            --phased_vcf_fn ${phased_vcf} \
            --no_phasing_for_fa False \
            --vcf_fn ${params.vcf_fn} \
            --gvcf ${params.GVCF}
        """
}


process aggregate_full_align_variants {
    // Sort and merge all "full alignment" variants
    label "wf_somatic_snv"
    cpus 2
    input:
        tuple val(meta), 
            path(vcfs, stageAs: "full_alignment/*"), 
            path(ref), 
            path(fai), 
            path(ref_cache), 
            path(contigs)
    output:
        tuple val(meta), path("full_alignment_${meta.sample}_${meta.type}.vcf.gz"), path("full_alignment_${meta.sample}_${meta.type}.vcf.gz.tbi"), emit: full_aln_vcf
        tuple val(meta), path("non_var.gvcf"), optional: true, emit: non_var_gvcf
    shell:
        '''
        pypy $(which clair3.py) SortVcf \
            --input_dir full_alignment \
            --output_fn full_alignment_!{meta.sample}_!{meta.type}.vcf \
            --sampleName !{meta.sample} \
            --ref_fn !{ref} \
            --contigs_fn !{contigs}

        if [ "$( bcftools index -n full_alignment_!{meta.sample}_!{meta.type}.vcf.gz )" -eq 0 ]; then
            echo "[INFO] Exit in full-alignment variant calling"
            exit 0
        fi

        # TODO: this could be a separate process
        if [ "!{params.GVCF}" == "true" ]; then
            pypy $(which clair3.py) SortVcf \
                --input_dir gvcf_tmp_path \
                --vcf_fn_suffix .tmp.gvcf \
                --output_fn non_var.gvcf \
                --sampleName !{meta.sample} \
                --ref_fn !{ref} \
                --contigs_fn !{contigs}
        fi
        '''
}


process merge_pileup_and_full_vars{
    // Merge VCFs
    label "wf_somatic_snv"
    cpus 2
    input:
        tuple val(meta), 
            path(pile_up_vcf), 
            path(pile_up_vcf_tbi), 
            path(full_aln_vcf), 
            path(full_aln_vcf_tbi), 
            path(ref), 
            path(fai), 
            path(ref_cache),
            path(beds, stageAs: "candidate_beds/*"),
            val(contig)
        
        
    output:
        tuple val(meta), val(contig), path("output/merge_${meta.sample}_${meta.type}_${contig}.vcf.gz"), path("output/merge_${meta.sample}_${meta.type}_${contig}.vcf.gz.tbi"), emit: merged_vcf
    shell:
        '''
        mkdir output
        echo "[INFO] 7/7 Merge pileup VCF and full-alignment VCF"
        pypy $(which clair3.py) MergeVcf \
            --pileup_vcf_fn !{pile_up_vcf} \
            --bed_fn_prefix candidate_beds \
            --full_alignment_vcf_fn !{full_aln_vcf} \
            --output_fn output/merge_!{meta.sample}_!{meta.type}_!{contig}.vcf \
            --platform ont \
            --print_ref_calls False \
            --gvcf !{params.GVCF} \
            --haploid_precise False \
            --haploid_sensitive False \
            --non_var_gvcf_fn non_var.gvcf \
            --ref_fn !{ref} \
            --ctgName !{contig}

        bgzip -c output/merge_!{meta.sample}_!{meta.type}_!{contig}.vcf > output/merge_!{meta.sample}_!{meta.type}_!{contig}.vcf.gz
        tabix output/merge_!{meta.sample}_!{meta.type}_!{contig}.vcf.gz
        '''
}


process aggregate_all_variants{
    label "wf_somatic_snv"
    cpus 4
    input:
        tuple val(meta), 
            val(contigs), 
            path(vcfs, stageAs: "merge_output/*"), 
            path(tbis, stageAs: "merge_output/*"), 
            path(gvcf, stageAs: "merge_outputs_gvcf/*"),
            path(ref), 
            path(fai), 
            path(ref_cache), 
            path(contigs)
    output:
        tuple val(meta), path("${meta.sample}_${meta.type}_germline.vcf.gz"), path("${meta.sample}_${meta.type}_germline.vcf.gz.tbi"), emit: final_vcf
    shell:
        '''
        prefix="merge"
        phase_vcf=!{params.clair3_phase_vcf}
        if [[ $phase_vcf == "true" ]]; then
            prefix="phased"
        fi
        ls merge_output/*.vcf.gz | parallel --jobs !{task.cpus} "bgzip -d {}"

        pypy $(which clair3.py) SortVcf \
            --input_dir merge_output \
            --vcf_fn_prefix $prefix \
            --output_fn !{meta.sample}_!{meta.type}_germline.vcf \
            --sampleName !{meta.sample} \
            --ref_fn !{ref} \
            --contigs_fn !{contigs}

        if [ "$( bcftools index -n !{meta.sample}_!{meta.type}_germline.vcf.gz )" -eq 0 ]; then
            echo "[INFO] Exit in all contigs variant merging"
            exit 0
        fi

        echo "[INFO] Finish calling, output file: !{meta.sample}_!{meta.type}_germline.vcf.gz"
        '''
}
