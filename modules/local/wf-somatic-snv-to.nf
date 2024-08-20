import groovy.json.JsonBuilder

// Define memory requirements for phasing
def req_mem = params.use_longphase ? [8.GB, 32.GB, 64.GB] : [4.GB, 8.GB, 12.GB]

process getVersions {
    label "wf_somatic_snv_to"
    cpus 1
    memory 2.GB
    output:
        path "versions.txt"
    script:
        """
        run_clairs_to --version | sed 's/ /,/' >> versions.txt
        bcftools --version | sed -n 1p | sed 's/ /,/' >> versions.txt
        whatshap --version | awk '{print "whatshap,"\$0}' >> versions.txt
        longphase --version | awk 'NR==1 {print "longphase,"\$2}' >> versions.txt
        """
}


process getParams {
    label "wf_somatic_snv_to"
    cpus 1
    memory 2.GB
    output:
        path "params.json"
    script:
        def paramsJSON = new JsonBuilder(params).toPrettyString()
        """
        # Output nextflow params object to JSON
        echo '$paramsJSON' > params.json
        """
}

//
// Individual ClairS-TO subprocesses to better handle the workflow.
//
// Prepare chunks and contigs list for parallel processing.
process wf_build_regions {
    label "wf_somatic_snv_to"
    cpus 1
    memory 4.GB
    input:
        tuple val(meta),
            path(tumor_bam, stageAs: "tumor/*"),
            path(tumor_bai, stageAs: "tumor/*")             
        tuple path(ref),
            path(fai),
            path(ref_cache),
            env(REF_PATH)
        val model
        path bed
        tuple path(typing_vcf), val(typing_opt)
            
    output:
        tuple val(meta), 
            path("${meta.sample}/tmp/CONTIGS"),
            emit: contigs_file
        tuple val(meta),
            path("${meta.sample}/tmp/CHUNK_LIST"),
            emit: chunks_file
        tuple val(meta),
            path("${meta.sample}/tmp/split_beds"),
            emit: split_beds
        tuple val(meta),
            path("${meta.sample}/tmp/split_indel_beds/"),
            emit: split_indel_beds
        tuple val(meta),
            path("${meta.sample}/tmp/CMD"),
            emit: cmd_file
    script:
        // Define additional inputs, such as target contigs, target region as bed file and more.
        def bedargs = params.bed ? "--bed_fn ${bed}" : ""
        def include_ctgs = params.include_all_ctgs ? "--include_all_ctgs" : ""
        def target_ctg = params.ctg_name == "EMPTY" ? "" : "--ctg_name ${params.ctg_name}"
        // Enable hybrid/genotyping mode if passed
        def typing_mode = typing_opt ? "${typing_opt} ${typing_vcf}" : ""
        """
        \${CLAIRS_PATH}/run_clairs_to ${target_ctg} ${include_ctgs} \\
            --tumor_bam_fn ${tumor_bam.getName()} \\
            --ref_fn ${ref} \\
            --platform ${model} \\
            --output_dir ./${meta.sample} \\
            --threads ${task.cpus} \\
            --snv_min_af ${params.snv_min_af} \\
            --min_coverage ${params.min_cov} \\
            --qual ${params.min_qual} \\
            --chunk_size 5000000 \\
            --dry_run \\
            ${bedargs} ${typing_mode}

        # Save empty file to prevent empty directory errors on AWS
        touch ${meta.sample}/tmp/split_beds/EMPTY
        """
}

// Extract candidate regions to process.
process clairs_to_extract_candidates {
    label "wf_somatic_snv_to"
    cpus 2
    memory { 6.GB * task.attempt }
    maxRetries 1
    errorStrategy {task.exitStatus in [137,140] ? 'retry' : 'finish'}
    input:
        tuple val(meta),
            val(region),
            path(tumor_bam),
            path(tumor_bai),
            path(split_beds),
            path(split_indel_beds),
            path(ref), 
            path(fai), 
            path(ref_cache), 
            env(REF_PATH),
            val(model)
        tuple path(typing_vcf), val(typing_opt)
    output:
        tuple val(meta),
            val("snv"),
            val(region),
            path(tumor_bam),
            path(tumor_bai),
            path(ref), 
            path(fai), 
            path(ref_cache), 
            env(REF_PATH),
            val(model),
            path("candidates/${region.contig}.*_snv"),
            emit: candidates_snvs,
            optional: true
        tuple val(meta),
            val("indel"),
            val(region),
            path(tumor_bam),
            path(tumor_bai),
            path(ref), 
            path(fai), 
            path(ref_cache), 
            env(REF_PATH),
            val(model),
            path("candidates/${region.contig}.*_indel"),
            emit: candidates_indels,
            optional: true

    script:
        // Call indels when r10 model is provided.
        // Enable hybrid/genotyping mode if passed
        def typing_mode = typing_opt ? "${typing_opt} ${typing_vcf}" : ""
        """
        # Create output folder structure
        for dir_name in candidates indels hybrid; do
            if [ -e \$dir_name ]; then
                rm -r dir_name
            fi
            mkdir -p \$dir_name
        done
        # Create candidates
        pypy3 \${CLAIRS_PATH}/clairs_to.py extract_candidates_calling \\
            --tumor_bam_fn ${tumor_bam} \\
            --ref_fn ${ref} \\
            --samtools samtools \\
            --snv_min_af ${params.snv_min_af} \\
            --indel_min_af ${params.indel_min_af} \\
            --chunk_id ${region.chunk_id} \\
            --chunk_num ${region.total_chunks} \\
            --ctg_name ${region.contig} \\
            --platform ont \\
            --min_coverage ${params.min_cov} \\
            --min_bq ${params.clairs_to_min_bq} \\
            --bed_fn ${split_beds}/${region.contig} \\
            --call_indels_only_in_these_regions ${split_indel_beds}/${region.contig} \\
            --candidates_folder candidates/ \\
            --output_depth True  \\
            --select_indel_candidates True \\
            ${typing_mode}
        """
}


// Create Tensors for Affirmative Model.
process clairs_to_create_affirmative_model_tensors {
    label "wf_somatic_snv_to"
    cpus 2
    memory { 4.GB * task.attempt }
    maxRetries 1
    errorStrategy {task.exitStatus in [137,140] ? 'retry' : 'finish'}
    input:
        tuple val(meta),
            val(variant_type),
            val(region),
            path(tumor_bam), 
            path(tumor_bai),
            path(ref), 
            path(fai), 
            path(ref_cache), 
            env(REF_PATH),
            val(model),
            path(candidate)
    output:
        tuple val(meta),
            val(variant_type),
            val(region),
            path(candidate),
            val(model),
            path("pileup_tensor_can_affirmative/"),
            emit: full_tensors
            
    script:
        """
        if [ -e pileup_tensor_can_affirmative ]; then rm -r pileup_tensor_can_affirmative; fi
        mkdir -p pileup_tensor_can_affirmative
        pypy3 \${CLAIRS_PATH}/clairs_to.py create_tensor_pileup_calling \\
            --tumor_bam_fn ${tumor_bam} \\
            --ref_fn ${ref} \\
            --ctg_name ${region.contig} \\
            --samtools samtools \\
            --min_bq 20 \\
            --candidates_bed_regions ${candidate} \\
            --tensor_can_fn pileup_tensor_can_affirmative/${candidate.getName()} \\
            --platform ont
        """
}

// Create Tensors for Negational Model.
process clairs_to_create_negational_model_tensors {
    label "wf_somatic_snv_to"
    cpus 2
    memory { 4.GB * task.attempt }
    maxRetries 1
    errorStrategy {task.exitStatus in [137,140] ? 'retry' : 'finish'}
    input:
        tuple val(meta), 
            val(variant_type),
            val(region),
            path(tumor_bam), 
            path(tumor_bai),
            path(ref), 
            path(fai), 
            path(ref_cache), 
            env(REF_PATH),
            val(model),
            path(candidate)
    output:
        tuple val(meta),
            val(variant_type),
            val(region),
            path(candidate),
            val(model),
            path("pileup_tensor_can_negational/"),
            emit: full_tensors
            
    script:
        """
        if [ -e pileup_tensor_can_negational ]; then rm -r pileup_tensor_can_negational; fi
        mkdir -p pileup_tensor_can_negational
        pypy3 \${CLAIRS_PATH}/clairs_to.py create_tensor_pileup_calling \\
            --tumor_bam_fn ${tumor_bam} \\
            --ref_fn ${ref} \\
            --ctg_name ${region.contig} \\
            --samtools samtools \\
            --min_bq 0 \\
            --candidates_bed_regions ${candidate} \\
            --tensor_can_fn pileup_tensor_can_negational/${candidate.getName()} \\
            --platform ont
        """
}

// Perform pileup variant prediction using the paired affirmative-negational tensors.
process clairs_to_predict_pileup {
    label "wf_somatic_snv_to"
    cpus 1
    memory { 6.GB * task.attempt }
    maxRetries 3
    errorStrategy {task.exitStatus in [134,137,140] ? 'retry' : 'finish'}
    input:
        tuple val(meta),
            val(variant_type),
            val(region),
            path(candidate),
            val(model),
            path("pileup_tensor_can_affirmative"),
            path("pileup_tensor_can_negational")
    output:
        tuple val(meta),
            val(variant_type),
            val(region),
            path(candidate),
            val(model),
            path("predict/${candidate.getName()}"), 
            optional: true
            
    script:
        // If requested, hybrid/genotyping mode, or vcf_fn are provided, then call also reference sites and germline sites.
        def print_ref = ""
        if (params.print_ref_calls || params.hybrid_mode_vcf || params.genotyping_mode_vcf || params.vcf_fn != 'EMPTY'){
            print_ref = "--show_ref"
        } 
        def print_ger = ""
        if (params.print_germline_calls || params.hybrid_mode_vcf || params.genotyping_mode_vcf || params.vcf_fn != 'EMPTY'){
            print_ger = "--show_germline"
        }
        def run_gpu = "--use_gpu False"
        // Define option and models if it is snv rather than indel
        def indel_opt = variant_type == "indel" ? "--disable_indel_calling False" : "--disable_indel_calling True"
        def model_path = variant_type == "indel" ? "${model}/indel" : "${model}"
        """
        if [ -e predict ]; then rm -r predict; fi
        mkdir -p predict/
        python3 \${CLAIRS_PATH}/clairs_to.py predict \\
            --tensor_fn_acgt pileup_tensor_can_affirmative/${candidate.getName()} \\
            --tensor_fn_nacgt pileup_tensor_can_negational/${candidate.getName()} \\
            --predict_fn predict/${candidate.getName()} \\
            --chkpnt_fn_acgt \${CLAIR_MODELS_PATH}/${model_path}/pileup_affirmative.pkl \\
            --chkpnt_fn_nacgt \${CLAIR_MODELS_PATH}/${model_path}/pileup_negational.pkl \\
            ${run_gpu} \\
            --platform ont \\
            --ctg_name ${region.contig} \\
            --pileup \\
            ${indel_opt} \\
            ${print_ref} \\
            ${print_ger}
        """
}

// Create pileup VCFs using the predictions.
process clairs_to_pileup {
    label "wf_somatic_snv_to"
    cpus 1
    memory { 4.GB * task.attempt }
    // Add also 134 as it is the error code associated with libomp.so already initialized.
    // Since this error is spurious and, so far, unpredictable, add it here to simply retry
    errorStrategy {task.exitStatus in [134,137,140] ? 'retry' : 'finish'}
    maxRetries 1
    input:
        tuple val(meta),
            val(variant_type),
            val(region),
            path(candidate),
            val(model),
            path("predicted/*"),
            path(ref), 
            path(fai), 
            path(ref_cache), 
            env(REF_PATH)

    output:
        tuple val(meta),
            val(variant_type),
            path("vcf_output/p_${candidate.getName()}.vcf"),
            emit: pileup_vcf,
            optional: true
            
    script:
        def indel_opt =  variant_type == "indel" ? "--disable_indel_calling False" : "--disable_indel_calling True"
        def likelihood_mat_path = variant_type == "indel" ? "${model}/indel" : "${model}"
        """
        if [ -e vcf_output ]; then rm -r vcf_output; fi
        mkdir -p vcf_output
        python3 \${CLAIRS_PATH}/clairs_to.py call_variants \\
            --predict_fn predicted/${candidate.getName()} \\
            --call_fn vcf_output/p_${candidate.getName()}.vcf \\
            --ref ${ref} \\
            --platform ont \\
            --likelihood_matrix_data \${CLAIR_MODELS_PATH}/${likelihood_mat_path}/likelihood_matrix.txt \\
            ${indel_opt}
        """
}


// Merge single-chunk pileups in a single VCF file.
process clairs_to_merge_pileup {
    label "wf_somatic_snv_to"
    cpus 1
    memory 4.GB
    input:
        tuple val(meta), 
            val(variant_type),
            path(vcfs, stageAs: 'vcf_output/*'), 
            path(contig_file), 
            path(ref), 
            path(fai), 
            path(ref_cache), 
            env(REF_PATH)
    output:
        tuple val(meta),
            val(variant_type),
            path("pileup_${variant_type}.vcf"),
            emit: full_vcf
            
    shell:
        '''
        pypy3 ${CLAIRS_PATH}/clairs_to.py sort_vcf \\
            --ref_fn !{ref} \\
            --contigs_fn !{contig_file} \\
            --input_dir vcf_output/ \\
            --vcf_fn_prefix p_ \\
            --output_fn pileup_!{variant_type}.vcf
        '''
}

// Use DBs to tag non-somatic variants.
process clairs_to_tag_non_somatic_db {
    label "wf_somatic_snv_to"
    cpus 1
    memory { 8.GB * task.attempt }
    maxRetries 1
    errorStrategy {task.exitStatus in [137,140] ? 'retry' : 'finish'}
    input:
        tuple val(meta),
            val(variant_type), 
            path('pileup.vcf'), 
            val(contig),
            path(ref), 
            path(fai), 
            path(ref_cache), 
            env(REF_PATH)
    output:
        tuple val(meta),
            val(variant_type),
            path("${variant_type}_nonsomatic_tagging_${contig}.vcf"),
            emit: full_vcf
            
    script:
        // For reasons beyond my understanding, changing the order of the databases causes the workflow to crash,
        // something that makes replacing a more automated approach to the input list less viable.
        // We now create them as a variable, allowing us in the future to expose them as a variable.
        def normal_dbs = "\${CLAIR_DBS_PATH}/gnomad.r2.1.af-ge-0.001.sites.vcf.gz,\${CLAIR_DBS_PATH}/dbsnp.b138.non-somatic.sites.vcf.gz,\${CLAIR_DBS_PATH}/1000g-pon.sites.vcf.gz"
        def normal_allele_matching = "True,True,False"
        """
        # Tag non-somatic variants
        pypy3 \${CLAIRS_PATH}/clairs_to.py nonsomatic_tagging \\
            --pileup_vcf_fn pileup.vcf \\
            --output_vcf_fn ./${variant_type}_nonsomatic_tagging_${contig}.vcf \\
            --ctg_name ${contig} \\
            --pypy3 pypy3 \\
            --parallel parallel \\
            --panel_of_normals ${normal_dbs} \\
            --panel_of_normals_require_allele_matching ${normal_allele_matching}
        """
}

// Merge tagged non-somatic variants.
process clairs_to_merge_tagged {
    label "wf_somatic_snv_to"
    cpus 1
    memory 4.GB
    input:
        tuple val(meta), 
            val(variant_type),
            path(vcfs, stageAs: 'vcf_tagged/*'),
            path(contig_file), 
            path(ref), 
            path(fai), 
            path(ref_cache), 
            env(REF_PATH)

    output:
        tuple val(meta),
            val(variant_type),
            path("${variant_type}_pileup_nonsomatic_tagging.vcf"),
            emit: full_vcf
            
    shell:
        '''
        pypy3 ${CLAIRS_PATH}/clairs_to.py sort_vcf \\
            --ref_fn !{ref} \\
            --contigs_fn !{contig_file} \\
            --input_dir vcf_tagged/ \\
            --vcf_fn_prefix !{variant_type}_nonsomatic_tagging_ \\
            --output_fn !{variant_type}_pileup_nonsomatic_tagging.vcf
        '''
}


// Select heterozygote sites for phasing.
process clairs_to_select_het_snps {
    label "wf_somatic_snv_to"
    cpus 2
    memory 4.GB
    input:
        tuple val(meta),
            val(variant_type),
            path("pileup.vcf"),
            val(contig),
            path(xam),
            path(xam_idx),
            path(ref), 
            path(fai), 
            path(ref_cache), 
            env(REF_PATH)
            
    output:
        tuple val(meta),
            val(variant_type),
            val(contig),
            path("split_folder/${contig}.vcf.gz"),
            path("split_folder/${contig}.vcf.gz.tbi"),
            path(xam),
            path(xam_idx),
            path(ref), 
            path(fai), 
            path(ref_cache),
            env(REF_PATH),
            emit: het_sites
    shell:
        '''
        pypy3 ${CLAIRS_PATH}/clairs_to.py select_hetero_snp_for_phasing \
            --tumor_vcf_fn pileup.vcf \
            --output_folder split_folder/ \
            --ctg_name !{contig}

        bgzip -c split_folder/!{contig}.vcf > split_folder/!{contig}.vcf.gz
        tabix split_folder/!{contig}.vcf.gz
        '''
}

// Run variant phasing on each contig using either longphase or whatshap.
process clairs_to_phase {
    label "wf_somatic_snv_to"
    cpus 4
    // Define memory from phasing tool and number of attempt
    memory { req_mem[task.attempt - 1] }
    maxRetries 2
    errorStrategy {task.exitStatus in [137,140] ? 'retry' : 'finish'}
    input:
        tuple val(meta),
            val(variant_type),
            val(contig),
            path(vcf),
            path(tbi),
            path(xam),
            path(xam_idx),
            path(ref), 
            path(fai), 
            path(ref_cache),
            env(REF_PATH)
            
    output:
        tuple val(meta),
            val(variant_type),
            val(contig),
            path("phased_${meta.sample}_${contig}.vcf.gz"), 
            path("phased_${meta.sample}_${contig}.vcf.gz.tbi"),
            path(xam),
            path(xam_idx),
            path(ref), 
            path(fai), 
            path(ref_cache),
            env(REF_PATH),
            emit: phased_data
    script:
        if (params.use_longphase)
        """
        echo "Using longphase for phasing"
        longphase phase --ont -o tmp \\
            -s ${vcf} -b ${xam} -r ${ref} -t ${task.cpus}
        bcftools sort -O z tmp.vcf > phased_${meta.sample}_${contig}.vcf.gz
        tabix -p vcf phased_${meta.sample}_${contig}.vcf.gz
        """
        else
        """
        echo "Using whatshap for phasing"
        whatshap phase \\
            --output phased_${meta.sample}_${contig}.vcf.gz \\
            --reference ${ref} \\
            --chromosome ${contig} \\
            --distrust-genotypes \\
            --ignore-read-groups \\
            ${vcf} \\
            ${xam}
        tabix -f -p vcf phased_${meta.sample}_${contig}.vcf.gz
        """
}

// Haplotag BAM file using either whatshap or longphase.
process clairs_to_haplotag {
    label "wf_somatic_snv_to"
    cpus 4
    memory 4.GB
    input:
        tuple val(meta),
            val(variant_type),
            val(contig),
            path(vcf),
            path(tbi),
            path(xam),
            path(xam_idx),
            path(ref), 
            path(fai), 
            path(ref_cache), 
            env(REF_PATH)
            
    output:
        tuple val(meta),
            path("${meta.sample}_${contig}.bam"), 
            path("${meta.sample}_${contig}.bam.bai"), 
            emit: haplotagged
    script:
        if (params.use_longphase)
        """
        longphase haplotag -o ${meta.sample}_${contig} --reference ${ref} --region ${contig} \\
            -t ${task.cpus} -s ${vcf} -b ${xam}
        samtools index -@${task.cpus} ${meta.sample}_${contig}.bam
        """
        else
        """
        whatshap haplotag \\
            --reference ${ref} \\
            --regions ${contig} \\
            --ignore-read-groups \\
            ${vcf} \\
            ${xam} \\
        | samtools view -b -1 -@3 -o ${meta.sample}_${contig}.bam##idx##${meta.sample}_${contig}.bam.bai --write-index --no-PG
        """
}

// Filter variants based on haplotype information.
process clairs_to_hap_filter {
    label "wf_somatic_snv_to"
    cpus params.haplotype_filter_threads
    memory { (2.GB * task.cpus) + 3.GB }
    input:
        tuple val(meta),
            path("bams/*"),
            path("bams/*"),
            val(variant_type),
            path("pileup.vcf"),
            path("pileup_nonsomatic_tagging.vcf"),
            val(contig),
            path(ref), 
            path(fai), 
            path(ref_cache), 
            env(REF_PATH)
    output:
        tuple val(meta),
            val(variant_type), 
            path("vcf_output/*pileup_filtering_${contig}.vcf"),
            emit: filtered_vcfs
            
    script:
        def debug = params.clairs_debug ? "--debug" : ""
        // If reference calls are requested, or it is a hybrid/genotyping mode, then add --show_ref option.
        def show_ref = params.print_ref_calls || params.hybrid_mode_vcf || params.genotyping_mode_vcf ? "--show_ref" : ""
        // Define if it is indels.
        def indel_opt = variant_type == 'indel' ? "--is_indel" : ""
        """
        pypy3 \${CLAIRS_PATH}/clairs_to.py haplotype_filtering \\
            --tumor_bam_fn bams/${meta.sample}_ \\
            --ref_fn ${ref} \\
            --pileup_vcf_fn pileup_nonsomatic_tagging.vcf \\
            --germline_vcf_fn pileup.vcf \\
            --output_dir vcf_output/ \\
            --output_vcf_fn vcf_output/${variant_type}_pileup_filtering_${contig}.vcf \\
            --threads ${task.cpus} \\
            --apply_haplotype_filtering True \\
            --ctg_name ${contig} \\
            ${indel_opt} ${show_ref} ${debug}
        """
}

// Merge single-contigs haplotype-filtered VCFs.
process clairs_to_merge_hapfilt {
    label "wf_somatic_snv_to"
    cpus 1
    memory 4.GB
    input:
        tuple val(meta),
            val(variant_type), 
            path("vcf_hapfilt/*"),
            path(ref), 
            path(fai), 
            path(ref_cache), 
            env(REF_PATH),
            path(contig_file) 

    output:
        tuple val(meta),
            val(variant_type),
            path("${variant_type}_pileup_filtering.vcf"),
            emit: hapfilt_vcf
            
    shell:
        '''
        pypy3 ${CLAIRS_PATH}/clairs_to.py sort_vcf \\
            --ref_fn !{ref} \\
            --contigs_fn !{contig_file} \\
            --input_dir vcf_hapfilt/ \\
            --vcf_fn_prefix !{variant_type}_pileup_filtering_ \\
            --output_fn !{variant_type}_pileup_filtering.vcf
        '''
}

// Run ClairS-TO postprocessing.
process clairs_to_postprocess {
    label "wf_somatic_snv_to"
    memory { 8.GB * task.attempt }
    maxRetries 1
    errorStrategy {task.exitStatus in [137,140] ? 'retry' : 'finish'}
    input:
        tuple val(meta),
            val(variant_type), 
            path("pileup_filtering.vcf"),
            path(ref), 
            path(fai), 
            path(ref_cache), 
            env(REF_PATH),
            path(CMD)
    output:
        tuple val(meta), 
            path("postprocessed_${variant_type}.vcf.gz"),
            path("postprocessed_${variant_type}.vcf.gz.tbi"),
            emit: final_vcf
            
    script:
        def debug = params.clairs_debug ? "--debug" : ""
        // If reference calls are requested, or it is a hybrid/genotyping mode, then add --show_ref option.
        def show_ref = params.print_ref_calls || params.hybrid_mode_vcf || params.genotyping_mode_vcf ? "--show_ref" : ""
        // Define if it is indels.
        def indel_opt = variant_type == 'indel' ? "--disable_indel_calling False" : "--disable_indel_calling True"
        """
        pypy3 \${CLAIRS_PATH}/clairs_to.py postprocess_vcf \\
            --ref_fn ${ref} \\
            --pileup_vcf_fn pileup_filtering.vcf \\
            --output_fn postprocessed_${variant_type}.vcf \\
            --platform ont \\
            --qual ${params.min_qual} \\
            --sample_name ${meta.sample} \\
            --cmdline ${CMD} \\
            --qual $params.clairs_to_qual \\
            --qual_cutoff_phaseable_region $params.qual_cutoff_phaseable_region \\
            --qual_cutoff_unphaseable_region $params.qual_cutoff_unphaseable_region \\
            ${indel_opt}
        """
}

// Sort-index VCFs.
process vcfSortIndex {
    label "wf_somatic_snv"
    cpus 2
    memory 4.GB

    input:
    tuple val(meta), path(vcf)

    output:
    tuple val(meta), 
        path("${meta.sample}.wf-snv.vcf.gz"),
        path("${meta.sample}.wf-snv.vcf.gz.tbi"),
        emit: final_vcf

    script:
    """
    bcftools sort ${vcf} | \\
        bcftools view --threads ${task.cpus} -O z > ${meta.sample}.wf-snv.vcf.gz && \\
        bcftools index --threads ${task.cpus} -t ${meta.sample}.wf-snv.vcf.gz
    """

}
