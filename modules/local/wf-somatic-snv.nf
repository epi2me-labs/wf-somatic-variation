import groovy.json.JsonBuilder

// See https://github.com/nextflow-io/nextflow/issues/1636
// This is the only way to publish files from a workflow whilst
// decoupling the publish from the process steps.
process output_snv {
    // publish inputs to output directory
    label "wf_somatic_snv"
    publishDir (
        params.out_dir,
        mode: "copy",
        saveAs: { dirname ? "$dirname/$fname" : fname }
    )
    input:
        tuple path(fname), val(dirname)
    output:
        path fname
    """
    """
}


process getVersions {
    label "wf_somatic_snv"
    cpus 1
    output:
        path "versions.txt"
    script:
        """
        run_clairs --version | sed 's/ /,/' >> versions.txt
        bcftools --version | sed -n 1p | sed 's/ /,/' >> versions.txt
        """
}


process getParams {
    label "wf_somatic_snv"
    cpus 1
    output:
        path "params.json"
    script:
        def paramsJSON = new JsonBuilder(params).toPrettyString()
        """
        # Output nextflow params object to JSON
        echo '$paramsJSON' > params.json
        """
}


process vcfStats {
    label "wf_somatic_snv"
    cpus 1
    input:
        tuple val(meta), path(vcf), path(index)
    output:
        tuple val(meta), path("${meta.sample}.stats")
    """
    bcftools stats $vcf > ${meta.sample}.stats
    """
}


process makeReport {
    input:
        tuple val(meta), 
            path(vcf), 
            path(tbi), 
            path("vcfstats.txt"), 
            path("spectra.csv"), 
            path(clinvar_vcf), 
            path("version.txt"), 
            path("params.json")
    output:
        path "*report.html", emit: html
    script:
        def report_name = "${params.sample_name}.wf-somatic-snv-report.html"
        def clinvar = clinvar_vcf.name == 'OPTIONAL_FILE' ? "" : "--clinvar_vcf ${clinvar_vcf}"
        wfversion = workflow.manifest.version
        if( workflow.commitId ){
            wfversion = workflow.commitId
        }
        def germline = params.germline ? "" : "--no_germline"
        def normal_vcf = params.normal_vcf ? "--normal_vcf ${file(params.normal_vcf).name}" : ""
        """
        workflow-glue report_snv \\
            $report_name \\
            --versions version.txt \\
            --params params.json \\
            --vcf_stats vcfstats.txt \\
            --vcf $vcf \\
            --mut_spectra spectra.csv \\
            ${clinvar} ${germline} ${normal_vcf}
        """
}


process lookup_clair3_model {
    label "wf_somatic_snv"
    input:
        path("lookup_table")
        val basecall_model
    output:
        stdout
    shell:
    '''
    resolve_clair3_model.py lookup_table '!{basecall_model}' 
    '''
}


//
// Individual ClairS subprocesses to better handle the workflow.
// 
process wf_build_regions {
    // Filters a VCF by contig, selecting only het SNPs.
    label "wf_somatic_snv"
    cpus 1
    input:
        tuple path(normal_bam, stageAs: "normal/*"), 
            path(normal_bai, stageAs: "normal/*"),
            path(tumor_bam, stageAs: "tumor/*"), 
            path(tumor_bai, stageAs: "tumor/*"), 
            val(meta) 
        tuple path(ref), path(fai), path(ref_cache)
        val model
        path bed
            
    output:
        tuple val(meta), path("${meta.sample}/tmp/CONTIGS"), emit: contigs_file
        tuple val(meta), path("${meta.sample}/tmp/CHUNK_LIST"), emit: chunks_file
        tuple val(meta), path("${meta.sample}/tmp/split_beds"), emit: split_beds
    script:
        def bedargs = params.bed ? "--bed_fn ${bed}" : ""
        def include_ctgs = params.include_all_ctgs ? "--include_all_ctgs" : ""
        def target_ctg = params.ctg_name == "EMPTY" ? "" : "--ctg_name ${params.ctg_name}"
        def indels_call = params.basecaller_cfg.startsWith('dna_r10') ? "--enable_indel_calling" : ""
        """
        \$CLAIRS_PATH/run_clairs ${target_ctg} ${include_ctgs} \\
            --tumor_bam_fn ${tumor_bam.getName()} \\
            --normal_bam_fn ${normal_bam.getName()} \\
            --ref_fn ${ref} \\
            --platform ${model} \\
            --output_dir ./${meta.sample} \\
            --threads ${task.cpus} \\
            --output_prefix pileup \\
            --snv_min_af ${params.snv_min_af} \\
            --min_coverage ${params.min_cov} \\
            --qual ${params.min_qual} \\
            --chunk_size 5000000 \\
            --dry_run \\
            ${bedargs} ${indels_call}

        # Save empty file to prevent empty directory errors on AWS
        touch ${meta.sample}/tmp/split_beds/EMPTY
        """
}

process clairs_select_het_snps {
    // Filters a VCF by contig, selecting only het SNPs.
    label "wf_somatic_snv"
    cpus 2
    input:
        tuple val(meta_tumor), 
            path(tumor_vcf, stageAs: "tumor.vcf.gz"), 
            path(tumor_tbi, stageAs: "tumor.vcf.gz.tbi"),
            val(meta_normal),
            path(normal_vcf, stageAs: "normal.vcf.gz"),
            path(normal_tbi, stageAs: "normal.vcf.gz.tbi"),
            
            val(contig)
            
    output:
        tuple val(meta_normal), val(contig), path("split_folder/${contig}.vcf.gz"), path("split_folder/${contig}.vcf.gz.tbi"), emit: normal_hets
        tuple val(meta_tumor), val(contig), path("split_folder/${contig}.vcf.gz"), path("split_folder/${contig}.vcf.gz.tbi"), emit: tumor_hets
    shell:
        '''
        pypy3 $CLAIRS_PATH/clairs.py select_hetero_snp_for_phasing \\
            --tumor_vcf_fn tumor.vcf.gz \\
            --normal_vcf_fn normal.vcf.gz \\
            --output_folder split_folder/ \\
            --ctg_name !{contig}

        bgzip -c split_folder/!{contig}.vcf > split_folder/!{contig}.vcf.gz
        tabix split_folder/!{contig}.vcf.gz
        '''
}

// Run haplotag on the bam files
process clairs_phase {
    // Filters a VCF by contig, selecting only het SNPs.
    label "wf_somatic_snv"
    cpus { params.use_longphase_intermediate ? 4 : 1 }
    input:
        tuple val(meta), 
            val(contig), 
            path(vcf), 
            path(tbi),
            path(bam), 
            path(bai), 
            path(ref), 
            path(fai), 
            path(ref_cache)
            
    output:
        tuple val(meta), 
            val(contig), 
            path("phased_${meta.sample}_${meta.type}_${contig}.vcf.gz"), 
            path("phased_${meta.sample}_${meta.type}_${contig}.vcf.gz.tbi"), 
            path(bam), 
            path(bai), 
            path(ref), 
            path(fai), 
            path(ref_cache),
            emit: phased_data
    script:
        if (params.use_longphase)
        """
        export REF_PATH=${ref_cache}/%2s/%2s/%s
        echo "Using longphase for phasing"
        # longphase needs decompressed 
        gzip -dc ${vcf} > variants.vcf
        longphase phase --ont -o phased_${meta.sample}_${meta.type}_${contig} \\
            -s variants.vcf -b ${bam} -r ${ref} -t ${task.cpus}
        bgzip phased_${meta.sample}_${meta.type}_${contig}.vcf
        tabix -f -p vcf phased_${meta.sample}_${meta.type}_${contig}.vcf.gz
        """
        else
        """
        export REF_PATH=!{ref_cache}/%2s/%2s/%s
        echo "Using whatshap for phasing"
        whatshap phase \\
            --output phased_${meta.sample}_${meta.type}_${contig}.vcf.gz \\
            --reference ${ref} \\
            --chromosome ${contig} \\
            --distrust-genotypes \\
            --ignore-read-groups \\
            ${vcf} \\
            ${bam}
        tabix -f -p vcf phased_${meta.sample}_${meta.type}_${contig}.vcf.gz
        """
}

// Run haplotag on the bam files
process clairs_haplotag {
    // Filters a VCF by contig, selecting only het SNPs.
    label "wf_somatic_snv"
    cpus 1
    input:
        tuple val(meta), 
            val(contig), 
            path(vcf), 
            path(tbi),
            path(bam), 
            path(bai), 
            path(ref), 
            path(fai), 
            path(ref_cache)
            
    output:
        tuple val(meta.sample), 
            val(contig), 
            path("${meta.sample}_${meta.type}_${contig}.bam"), 
            path("${meta.sample}_${meta.type}_${contig}.bam.bai"), 
            val(meta),
            emit: phased_data
    shell:
        '''
        export REF_PATH=!{ref_cache}/%2s/%2s/%s
        whatshap haplotag \\
            --output !{meta.sample}_!{meta.type}_!{contig}.bam \\
            --reference !{ref} \\
            --regions !{contig} \\
            --ignore-read-groups \\
            !{vcf} \\
            !{bam}
        samtools index !{meta.sample}_!{meta.type}_!{contig}.bam
        '''
}


// Extract candidate regions
process clairs_extract_candidates {
    // Filters a VCF by contig, selecting only het SNPs.
    label "wf_somatic_snv"
    cpus 2
    input:
        tuple val(meta),
            val(region),
            path(normal_bam, stageAs: "normal/*"),
            path(normal_bai, stageAs: "normal/*"),
            path(tumor_bam, stageAs: "tumor/*"),
            path(tumor_bai, stageAs: "tumor/*"),
            path(split_beds),
            path(ref), 
            path(fai), 
            path(ref_cache),
            val(model)
    output:
    tuple val(meta),
            val(region),
            path('candidates/CANDIDATES_FILE_*'),
            path("candidates/${region.contig}.*"),
            emit: candidates_snvs,
            optional: true
    tuple val(meta),
            val(region),
            path('indels/INDEL_CANDIDATES_FILE_*'),
            path("indels/${region.contig}.*_indel"),
            emit: candidates_indels,
            optional: true

    script:
        def bedfile = params.bed ? "" : ""
        def indel_min_af = model == 'ont_r10' ? "--indel_min_af ${params.indel_min_af}" : "--indel_min_af 1.00"
        def select_indels = model == 'ont_r10' ? "--select_indel_candidates True" : "--select_indel_candidates False"
        """
        export REF_PATH=${ref_cache}/%2s/%2s/%s        
        # Create output folder structure
        mkdir candidates; mkdir indels
        # Create candidates
        pypy3 \$CLAIRS_PATH/clairs.py extract_pair_candidates \\
            --tumor_bam_fn ${tumor_bam.getName()} \\
            --normal_bam_fn ${normal_bam.getName()} \\
            --ref_fn ${ref} \\
            --samtools samtools \\
            --snv_min_af ${params.snv_min_af} \\
            ${indel_min_af} \\
            --chunk_id ${region.chunk_id} \\
            --chunk_num ${region.total_chunks} \\
            --ctg_name ${region.contig} \\
            --platform ont \\
            --min_coverage ${params.min_cov} \\
            --bed_fn ${split_beds}/${region.contig} \\
            --candidates_folder candidates/ \\
            ${select_indels} \\
            --output_depth True 
        for i in `ls candidates/INDEL_*`; do
            echo "Moved \$i"
            mv \$i indels/
        done
        for i in `ls candidates/*_indel`; do 
            echo "Moved \$i"
            mv \$i indels/
        done
        """
}

// Create pileup Paired Tensors 
process clairs_create_paired_tensors {
    // Filters a VCF by contig, selecting only het SNPs.
    label "wf_somatic_snv"
    cpus 2
    input:
        tuple val(meta),
            val(region),
            path(normal_bam, stageAs: "normal/*"),
            path(normal_bai, stageAs: "normal/*"),
            path(tumor_bam, stageAs: "tumor/*"),
            path(tumor_bai, stageAs: "tumor/*"),
            path(split_beds),
            path(ref), 
            path(fai), 
            path(ref_cache),
            val(model),
            path(candidate),
            path(intervals)
    output:
        tuple val(meta),
            val(region),
            path(normal_bam),
            path(normal_bai),
            path(tumor_bam),
            path(tumor_bai),
            path(split_beds),
            path(ref), 
            path(fai), 
            path(ref_cache),
            val(model),
            path(candidate),
            path(intervals),
            path("tmp/pileup_tensor_can/")
            
    script:
        """
        export REF_PATH=${ref_cache}/%2s/%2s/%s
        mkdir -p tmp/pileup_tensor_can
        pypy3 \$CLAIRS_PATH/clairs.py create_pair_tensor_pileup \\
            --normal_bam_fn ${normal_bam.getName()} \\
            --tumor_bam_fn ${tumor_bam.getName()} \\
            --ref_fn ${ref} \\
            --samtools samtools \\
            --ctg_name ${region.contig} \\
            --candidates_bed_regions ${intervals} \\
            --tensor_can_fn tmp/pileup_tensor_can/${intervals.getName()} \\
            --platform ont
        """
}


// Predict pileup variants 
process clairs_predict_pileup {
    // Filters a VCF by contig, selecting only het SNPs.
    label "wf_somatic_snv"
    label "avx2"
    cpus 1
    input:
        tuple val(meta),
            val(region),
            path(normal_bam, stageAs: "normal/*"),
            path(normal_bai, stageAs: "normal/*"),
            path(tumor_bam, stageAs: "tumor/*"),
            path(tumor_bai, stageAs: "tumor/*"),
            path(split_beds),
            path(ref), 
            path(fai), 
            path(ref_cache),
            val(model),
            path(candidate),
            path(intervals),
            path(tensor)
    output:
        tuple val(meta), path("vcf_output/p_${intervals.getName()}.vcf"), optional: true
            
    script:
        def print_ref = params.print_ref_calls ? "--show_ref" : ""
        def print_ger = params.print_germline_calls ? "--show_germline" : ""
        def run_gpu = "--use_gpu False"
        """
        export REF_PATH=${ref_cache}/%2s/%2s/%s
        mkdir vcf_output/
        python3 \$CLAIRS_PATH/clairs.py predict \\
            --tensor_fn ${tensor}/${intervals.getName()} \\
            --call_fn vcf_output/p_${intervals.getName()}.vcf \\
            --chkpnt_fn \${CLAIR_MODELS_PATH}/${model}/pileup.pkl \\
            --platform ont \\
            ${run_gpu} \\
            --ctg_name ${region.contig} \\
            --pileup \\
            ${print_ref} \\
            ${print_ger}
        """
}

// Merge predicted pileup variants 
process clairs_merge_pileup {
    // Filters a VCF by contig, selecting only het SNPs.
    label "wf_somatic_snv"
    cpus 1
    input:
        tuple val(meta), 
            path(vcfs, stageAs: 'vcf_output/*'), 
            path(contig_file)
        tuple path(ref), 
            path(fai), 
            path(ref_cache)
    output:
        tuple val(meta), path("pileup.vcf"), emit: pileup_vcf
            
    shell:
        '''
        export REF_PATH=!{ref_cache}/%2s/%2s/%s
        pypy3 $CLAIRS_PATH/clairs.py sort_vcf \\
            --ref_fn !{ref} \\
            --contigs_fn !{contig_file} \\
            --input_dir vcf_output/ \\
            --vcf_fn_prefix p_ \\
            --output_fn pileup.vcf
        '''
}


// Create full-alignment Paired Tensors 
process clairs_create_fullalignment_paired_tensors {
    // Filters a VCF by contig, selecting only het SNPs.
    label "wf_somatic_snv"
    cpus 2
    input:
        tuple val(sample), 
            val(contig), 
            path(tumor_bam, stageAs: "tumor/*"), 
            path(tumor_bai, stageAs: "tumor/*"),
            val(meta), 
            path(normal_bam, stageAs: "normal/*"), 
            path(normal_bai, stageAs: "normal/*"), 
            path(intervals),
            path(ref), 
            path(fai), 
            path(ref_cache),
            val(model)
    output:
        tuple val(meta),
            val(contig),
            path(normal_bam),
            path(normal_bai),
            path(tumor_bam),
            path(tumor_bai),
            path(ref), 
            path(fai), 
            path(ref_cache),
            path(intervals),
            path("fa_tensor_can/"),
            val(model),
            emit: full_tensors
            
    script:
        """
        mkdir fa_tensor_can
        export REF_PATH=${ref_cache}/%2s/%2s/%s
        pypy3 \$CLAIRS_PATH/clairs.py create_pair_tensor \\
            --samtools samtools \\
            --normal_bam_fn ${normal_bam.getName()} \\
            --tumor_bam_fn ${tumor_bam.getName()} \\
            --ref_fn ${ref} \\
            --ctg_name ${contig} \\
            --candidates_bed_regions ${intervals} \\
            --tensor_can_fn fa_tensor_can/${intervals.getName()} \\
            --platform ont
        """
}



// Predict full-alignment variants 
process clairs_predict_full {
    // Filters a VCF by contig, selecting only het SNPs.
    label "wf_somatic_snv"
    label "avx2"
    cpus 1
    input:
        tuple val(meta),
            val(contig),
            path(normal_bam, stageAs: "normal/*"),
            path(normal_bai, stageAs: "normal/*"),
            path(tumor_bam, stageAs: "tumor/*"),
            path(tumor_bai, stageAs: "tumor/*"),
            path(ref), 
            path(fai), 
            path(ref_cache),
            path(intervals),
            path(tensor),
            val(model)
    output:
        tuple val(meta), path("vcf_output/fa_${intervals.getName()}.vcf"), emit: full_vcfs, optional: true
            
    script:
        def print_ref = params.print_ref_calls ? "--show_ref" : ""
        def print_ger = params.print_germline_calls ? "--show_germline" : ""
        def run_gpu = "--use_gpu False"
        """
        export REF_PATH=${ref_cache}/%2s/%2s/%s
        mkdir vcf_output/
        python3 \$CLAIRS_PATH/clairs.py predict \\
            --tensor_fn ${tensor}/${intervals.getName()} \\
            --call_fn vcf_output/fa_${intervals.getName()}.vcf \\
            --chkpnt_fn \${CLAIR_MODELS_PATH}/${model}/full_alignment.pkl \\
            --platform ont \\
            ${run_gpu} \\
            --ctg_name ${contig} \\
            ${print_ref} ${print_ger}
        """
}


// Merge full-alignment variants 
process clairs_merge_full {
    // Filters a VCF by contig, selecting only het SNPs.
    label "wf_somatic_snv"
    cpus 1
    input:
        tuple val(meta), 
            path(vcfs, stageAs: 'vcf_output/*'), 
            path(contig_file), 
            path(ref), 
            path(fai), 
            path(ref_cache)
    output:
        tuple val(meta), path("full_alignment.vcf"), emit: full_vcf
            
    shell:
        '''
        export REF_PATH=!{ref_cache}/%2s/%2s/%s
        pypy3 $CLAIRS_PATH/clairs.py sort_vcf \\
            --ref_fn !{ref} \\
            --contigs_fn !{contig_file} \\
            --input_dir vcf_output/ \\
            --vcf_fn_prefix fa_ \\
            --output_fn full_alignment.vcf
        '''
}

// Filter variants
process clairs_full_hap_filter {
    // Filters a VCF by contig, selecting only het SNPs.
    label "wf_somatic_snv"
    label "avx2"
    cpus params.haplotype_filter_threads
    input:
        tuple val(meta),
            path(tumor_bams, stageAs: "bams/*"),
            path(tumor_bai, stageAs: "bams/*"),
            path(germline_vcf),
            path(pileup_vcf),
            path(full_alignment_vcf),
            path(ref), 
            path(fai), 
            path(ref_cache)
    output:
        tuple val(meta), 
            path("vcf_output/pileup_filter.vcf"),
            path("vcf_output/full_alignment_filter.vcf"),
            path(ref), 
            path(fai), 
            path(ref_cache),
            emit: filtered_vcfs
            
    script:
        def debug = params.clairs_debug ? "--debug" : ""
        def germline = germline_vcf.baseName != 'OPTIONAL_FILE' ? "--germline_vcf_fn ${germline_vcf}" : "--germline_vcf_fn None"
        """
        mkdir vcf_output/
        export REF_PATH=${ref_cache}/%2s/%2s/%s
        pypy3 \$CLAIRS_PATH/clairs.py haplotype_filtering \\
            --tumor_bam_fn bams/${meta.sample}_${meta.type}_ \\
            --ref_fn ${ref} \\
            ${germline} \\
            --pileup_vcf_fn ${pileup_vcf} \\
            --full_alignment_vcf_fn ${full_alignment_vcf} \\
            --output_dir vcf_output/ \\
            --samtools samtools \\
            --threads ${task.cpus} \\
            ${debug}
        """
}

// Final merging of the sites
process clairs_merge_final {
    // Filters a VCF by contig, selecting only het SNPs.
    label "wf_somatic_snv"
    cpus 1
    input:
        tuple val(meta),
            path(pileup_filter),
            path(full_alignment_filter),
            path(ref), 
            path(fai), 
            path(ref_cache)
    output:    
        tuple val(meta), path("${meta.sample}_somatic_snv.vcf.gz"), emit: pileup_vcf
        tuple val(meta), path("${meta.sample}_somatic_snv.vcf.gz.tbi"), emit: pileup_tbi
            
    script:
        """
        export REF_PATH=${ref_cache}/%2s/%2s/%s
        pypy3 \$CLAIRS_PATH/clairs.py merge_vcf \\
            --ref_fn ${ref} \\
            --pileup_vcf_fn ${pileup_filter} \\
            --full_alignment_vcf_fn ${full_alignment_filter} \\
            --output_fn ${meta.sample}_somatic_snv.vcf \\
            --platform ont \\
            --qual ${params.min_qual} \\
            --sample_name ${meta.sample}

        # Check if the number of variants is > 0
        if [ "\$( bcftools index -n ${meta.sample}_somatic_snv.vcf.gz )" -eq 0 ]; \\
        then echo "[INFO] Exit in pileup variant calling"; exit 1; fi
        """
}

/*
*  Add modules to call indels
*/ 
// TODO: Several of these modules can probably be replaced with a modified version of one module, then imported with several aliases and appropriate switches
// Create indels pileup Paired Tensors 
process clairs_create_paired_tensors_indels {
    // Filters a VCF by contig, selecting only het SNPs.
    label "wf_somatic_snv"
    cpus 2
    input:
        tuple val(meta),
            val(region),
            path(normal_bam, stageAs: "normal/*"),
            path(normal_bai, stageAs: "normal/*"),
            path(tumor_bam, stageAs: "tumor/*"),
            path(tumor_bai, stageAs: "tumor/*"),
            path(split_beds),
            path(ref), 
            path(fai), 
            path(ref_cache),
            val(model),
            path(candidate),
            path(intervals)
    output:
        tuple val(meta),
            val(region),
            path(normal_bam),
            path(normal_bai),
            path(tumor_bam),
            path(tumor_bai),
            path(split_beds),
            path(ref), 
            path(fai), 
            path(ref_cache),
            val(model),
            path(candidate),
            path(intervals),
            path("tmp/indel_pileup_tensor_can/")
            
    script:
        """
        export REF_PATH=${ref_cache}/%2s/%2s/%s
        mkdir -p tmp/indel_pileup_tensor_can
        pypy3 \$CLAIRS_PATH/clairs.py create_pair_tensor_pileup \\
            --normal_bam_fn ${normal_bam.getName()} \\
            --tumor_bam_fn ${tumor_bam.getName()} \\
            --ref_fn ${ref} \\
            --samtools samtools \\
            --ctg_name ${region.contig} \\
            --candidates_bed_regions ${intervals} \\
            --tensor_can_fn tmp/indel_pileup_tensor_can/${intervals.getName()} \\
            --platform ont
        """
}

// Predict indels pileup 
process clairs_predict_pileup_indel {
    // Filters a VCF by contig, selecting only het SNPs.
    label "wf_somatic_snv"
    label "avx2"
    cpus 1
    input:
        tuple val(meta),
            val(region),
            path(normal_bam, stageAs: "normal/*"),
            path(normal_bai, stageAs: "normal/*"),
            path(tumor_bam, stageAs: "tumor/*"),
            path(tumor_bai, stageAs: "tumor/*"),
            path(split_beds),
            path(ref), 
            path(fai), 
            path(ref_cache),
            val(model),
            path(candidate),
            path(intervals),
            path(tensor)
    output:
        tuple val(meta), path("vcf_output/indel_p_${intervals.getName()}.vcf"), optional: true
            
    script:
        def print_ref = params.print_ref_calls ? "--show_ref" : ""
        def print_ger = params.print_germline_calls ? "--show_germline" : ""
        def run_gpu = "--use_gpu False"
        """
        export REF_PATH=${ref_cache}/%2s/%2s/%s
        mkdir vcf_output/
        python3 \$CLAIRS_PATH/clairs.py predict \\
            --tensor_fn ${tensor}/${intervals.getName()} \\
            --call_fn vcf_output/indel_p_${intervals.getName()}.vcf \\
            --chkpnt_fn \${CLAIR_MODELS_PATH}/${model}/indel/pileup.pkl \\
            --platform ont \\
            ${run_gpu} \\
            --ctg_name ${region.contig} \\
            --pileup \\
            --enable_indel_calling True \\
            ${print_ref} \\
            ${print_ger}
        """
}

// Merge predicted indels pileup
process clairs_merge_pileup_indels {
    // Filters a VCF by contig, selecting only het SNPs.
    label "wf_somatic_snv"
    cpus 1
    input:
        tuple val(meta), 
            path(vcfs, stageAs: 'vcf_output/*'), 
            path(contig_file)
        tuple path(ref), 
            path(fai), 
            path(ref_cache)
    output:
        tuple val(meta), path("indels_pileup.vcf"), emit: pileup_vcf
            
    shell:
        '''
        export REF_PATH=!{ref_cache}/%2s/%2s/%s
        pypy3 $CLAIRS_PATH/clairs.py sort_vcf \\
            --ref_fn !{ref} \\
            --contigs_fn !{contig_file} \\
            --input_dir vcf_output/ \\
            --vcf_fn_prefix indel_p_ \\
            --output_fn indels_pileup.vcf
        '''
}

// Create full-alignment indels Paired Tensors 
process clairs_create_fullalignment_paired_tensors_indels {
    // Filters a VCF by contig, selecting only het SNPs.
    label "wf_somatic_snv"
    cpus 2
    input:
        tuple val(sample), 
            val(contig), 
            path(tumor_bam, stageAs: "tumor/*"), 
            path(tumor_bai, stageAs: "tumor/*"),
            val(meta), 
            path(normal_bam, stageAs: "normal/*"), 
            path(normal_bai, stageAs: "normal/*"), 
            path(intervals),
            path(ref), 
            path(fai), 
            path(ref_cache),
            val(model)
    output:
        tuple val(meta),
            val(contig),
            path(normal_bam),
            path(normal_bai),
            path(tumor_bam),
            path(tumor_bai),
            path(ref), 
            path(fai), 
            path(ref_cache),
            path(intervals),
            path("fa_tensor_can_indels/"),
            val(model),
            emit: full_tensors
            
    script:
        """
        mkdir fa_tensor_can_indels
        export REF_PATH=${ref_cache}/%2s/%2s/%s
        pypy3 \$CLAIRS_PATH/clairs.py create_pair_tensor \\
            --samtools samtools \\
            --normal_bam_fn ${normal_bam.getName()} \\
            --tumor_bam_fn ${tumor_bam.getName()} \\
            --ref_fn ${ref} \\
            --ctg_name ${contig} \\
            --candidates_bed_regions ${intervals} \\
            --tensor_can_fn fa_tensor_can_indels/${intervals.getName()} \\
            --platform ont
        """
}

// Predict full-alignment indels
process clairs_predict_full_indels {
    // Filters a VCF by contig, selecting only het SNPs.
    label "wf_somatic_snv"
    label "avx2"
    cpus 1
    input:
        tuple val(meta),
            val(contig),
            path(normal_bam, stageAs: "normal/*"),
            path(normal_bai, stageAs: "normal/*"),
            path(tumor_bam, stageAs: "tumor/*"),
            path(tumor_bai, stageAs: "tumor/*"),
            path(ref), 
            path(fai), 
            path(ref_cache),
            path(intervals),
            path(tensor),
            val(model)
    output:
        tuple val(meta), path("vcf_output/indels_fa_${intervals.getName()}.vcf"), emit: full_vcfs, optional: true
            
    script:
        def print_ref = params.print_ref_calls ? "--show_ref" : ""
        def print_ger = params.print_germline_calls ? "--show_germline" : ""
        def run_gpu = "--use_gpu False"
        """
        export REF_PATH=${ref_cache}/%2s/%2s/%s
        mkdir vcf_output/
        python3 \$CLAIRS_PATH/clairs.py predict \\
            --tensor_fn ${tensor}/${intervals.getName()} \\
            --call_fn vcf_output/indels_fa_${intervals.getName()}.vcf \\
            --chkpnt_fn \${CLAIR_MODELS_PATH}/${model}/indel/full_alignment.pkl \\
            --platform ont \\
            ${run_gpu} \\
            --ctg_name ${contig} \\
            --enable_indel_calling True \\
            ${print_ref} ${print_ger}
        """
}

// Merge full-alignment variants 
process clairs_merge_full_indels {
    // Filters a VCF by contig, selecting only het SNPs.
    label "wf_somatic_snv"
    cpus 1
    input:
        tuple val(meta), 
            path(vcfs, stageAs: 'vcf_output/*'), 
            path(contig_file), 
            path(ref), 
            path(fai), 
            path(ref_cache)
    output:
        tuple val(meta), path("indels_full_alignment.vcf"), emit: full_vcf
            
    shell:
        '''
        export REF_PATH=!{ref_cache}/%2s/%2s/%s
        pypy3 $CLAIRS_PATH/clairs.py sort_vcf \\
            --ref_fn !{ref} \\
            --contigs_fn !{contig_file} \\
            --input_dir vcf_output/ \\
            --vcf_fn_prefix indels_fa_ \\
            --output_fn indels_full_alignment.vcf
        '''
}

// Final merging of the sites
process clairs_merge_final_indels {
    // Filters a VCF by contig, selecting only het SNPs.
    label "wf_somatic_snv"
    cpus 1
    input:
        tuple val(meta),
            path(pileup_indels),
            path(full_alignment_indels),
            path(ref), 
            path(fai), 
            path(ref_cache)
    output:    
        tuple val(meta), path("${meta.sample}_somatic_indels.vcf.gz"), emit: indel_vcf
        tuple val(meta), path("${meta.sample}_somatic_indels.vcf.gz.tbi"), emit: indel_tbi
            
    script:
        """
        export REF_PATH=${ref_cache}/%2s/%2s/%s
        pypy3 \$CLAIRS_PATH/clairs.py merge_vcf \\
            --ref_fn ${ref} \\
            --pileup_vcf_fn ${pileup_indels} \\
            --full_alignment_vcf_fn ${full_alignment_indels} \\
            --output_fn ${meta.sample}_somatic_indels.vcf \\
            --platform ont \\
            --qual ${params.min_qual} \\
            --sample_name ${meta.sample} \\
            --enable_indel_calling True \\
            --indel_calling

        # Check if the number of variants is > 0
        if [ "\$( bcftools index -n ${meta.sample}_somatic_indels.vcf.gz )" -eq 0 ]; \\
        then echo "[INFO] Exit in pileup variant calling"; exit 1; fi
        """
}
/*
 * Post-processing processes
 */
// Merge full-alignment variants 
process clairs_merge_snv_and_indels {
    // Filters a VCF by contig, selecting only het SNPs.
    label "wf_somatic_snv"
    cpus 1
    input:
        tuple val(meta), path(vcfs, stageAs: 'VCFs/*'), path(tbis, stageAs: 'VCFs/*')
    output:
        tuple val(meta), path("${meta.sample}_somatic.vcf.gz"), emit: pileup_vcf
        tuple val(meta), path("${meta.sample}_somatic.vcf.gz.tbi"), emit: pileup_tbi
            
    shell:
        '''
        bcftools concat -a VCFs/*.vcf.gz | \\
            bcftools sort -m 2G -T ./ -O z > !{meta.sample}_somatic.vcf.gz && \\
            tabix -p vcf !{meta.sample}_somatic.vcf.gz
        '''
}

// Annotate the mutation counts
process change_count {
    // Filters a VCF by contig, selecting only het SNPs.
    cpus 1
    input:
        tuple val(meta),
            path("input.vcf.gz"),
            path("input.vcf.gz.tbi"),
            path(ref), 
            path(fai), 
            path(ref_cache)
    output:    
        tuple val(meta), path("${meta.sample}_somatic.vcf.gz"), emit: mutype_vcf
        tuple val(meta), path("${meta.sample}_somatic.vcf.gz"), emit: mutype_tbi
        tuple val(meta), path("${meta.sample}_changes.csv"), emit: changes
        tuple val(meta), path("${meta.sample}_changes.json"), emit: changes_json
            
    script:
        """
        workflow-glue annotate_mutations input.vcf.gz ${meta.sample}_somatic.vcf.gz --json -k 3 --genome ${ref}
        tabix -p vcf ${meta.sample}_somatic.vcf.gz
        """
}