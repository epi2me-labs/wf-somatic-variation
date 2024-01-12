import groovy.json.JsonBuilder

// Generate BWA-mem index for the insert_classify stage
process bwa_index {
    label "wf_somatic_sv"
    input:
        tuple path(ref), path(fai), path(cram_cache), env(REF_PATH)
    output:
        tuple path(ref), path(fai), path(cram_cache), env(REF_PATH), path("${ref}.amb"), path("${ref}.ann"), path("${ref}.bwt"), path("${ref}.pac"), path("${ref}.sa"), emit: bwa_ref
    script:
    """
    bwa index -a bwtsw ${ref}
    """
}


// NOTE VCF entries for alleles with no support are removed to prevent them from
//      breaking downstream parsers that do not expect them
process nanomonsv_parse {
    label "wf_somatic_sv"
    cpus 1
    input:
        tuple path(xam), path(xam_idx), val(meta)
        tuple path(ref), path(fai), path(ref_cache), env(REF_PATH)
        
    output:
        tuple val(meta), path(xam), path(xam_idx), path("${meta.sample}_${meta.type}/*")
    script:
    """
    mkdir ${meta.sample}_${meta.type}/
    nanomonsv parse ${xam} "${meta.sample}_${meta.type}/${meta.sample}"
    """
}


// Run nanomonsv get on paired tumor/normal samples
process nanomonsv_get {
    label "wf_somatic_sv"
    label "avx2"
    // Use at least 2 cores
    cpus params.nanomonsv_get_threads < 2 ? 2 : params.nanomonsv_get_threads
    input:
        tuple val(meta),
            path(xam_normal, stageAs: "norm/*"),
            path(xam_idx_normal, stageAs: "norm/*"),
            path(parsed_normal, stageAs: "parsed_normal/*"),
            path(xam_tumor, stageAs: "tum/*"),
            path(xam_idx_tumor, stageAs: "tum/*"),
            path(parsed_tumor, stageAs: "parsed_tumor/*")
        tuple path(ref), 
            path(fai), 
            path(ref_cache), 
            env(REF_PATH)
        tuple path(control_panel), val(control_root)
    output:
        tuple val(meta), path("parsed_tumor/${meta.sample}.nanomonsv.result.txt"), emit: txt
        tuple val(meta), path("parsed_tumor/${meta.sample}.nanomonsv.result.vcf"), emit: vcf
        tuple val(meta), path("parsed_tumor/${meta.sample}.nanomonsv.sbnd.result.txt"), emit: single_breakend
        tuple val(meta), path("parsed_tumor/${meta.sample}.nanomonsv.supporting_read.txt"), emit: read_lists

    script:
    def ncores = task.cpus - 1 // Use cpu-1 to ensure racon/minimap run as subprocess
    def qv = params.qv ? "--qv${params.qv}" : ""
    // If bam_normal is provided, then process it.
    def control_bam = params.bam_normal ? "--control_bam ${xam_normal} --control_prefix parsed_normal/${meta.sample}" : ""
    // Use input control panel for downstream analyses
    def use_control_panel = params.control_panel ? "--control_panel_prefix ${control_panel}/${control_root}" : ""
    """
    nanomonsv get \\
        "parsed_tumor/${meta.sample}" \\
        ${xam_tumor} \\
        ${ref} \\
        --min_indel_size ${params.min_sv_length} \\
        ${control_bam} \\
        --processes ${ncores} \\
        --single_bnd \\
        --use_racon \\
        --max_memory_minimap2 2 \\
        ${qv} ${use_control_panel}
    """
}

// Filter SVs in tandem repeat
process nanomonsv_filter {
    label "wf_somatic_sv"
    cpus 1
    input:
        tuple val(meta), path(txt)
        tuple val(meta), path(vcf)
        file trbed
        file trtbi
    output:
        tuple val(meta), path(vcf), path("${txt.baseName}.filter.txt"), emit: filtered
    script:
    """
    add_simple_repeat.py \
        $txt \
        ${txt.baseName}.filter.txt \
        $trbed
    """
}

// Annotate filters
process annotate_filter {
    cpus 1
    input:
        tuple val(meta), 
            path(vcf),
            path(filter)
    output:
        tuple val(meta), path("${vcf.baseName}.filter.vcf"), emit: vcf
        tuple val(meta), path(filter), emit: txt
    script:
    """
    workflow-glue extract_filtered_svs \\
        --in_vcf ${vcf} \\
        --filtered ${filter} \\
        --out_vcf ${vcf.baseName}.filter.vcf
    """
}


// Classify mobile elements
// Still in alpha stage, quite buggy
process nanomonsv_classify {
    label "wf_somatic_sv"
    input:
        tuple val(meta), path(txt), path(vcf)
        tuple path(ref), path(fai), path(cram_cache), env(REF_PATH), path(amb), path(ann), path(bwt), path(pac), path(sa)
        val n_valid_inserts
    output:
        tuple val(meta), path(txt), path(vcf), path("${txt.baseName}.annot.txt"), emit: txt
    script:
    if (n_valid_inserts > 0)
    """
    nanomonsv insert_classify --genome_id ${meta.genome_build} ${txt} ${txt.baseName}.annot.txt ${ref} 
    """
    else
    """
    ln -s ${txt} ${txt.baseName}.annot.txt 
    """
}

// Annotate classify
process annotate_classify {
    input:
        tuple val(meta), path(txt), path(vcf), path(annot_txt)
    output:
        tuple val(meta), path(annot_txt), path("${vcf.baseName}.annot.vcf"), emit: annotated
    script:
    """
    workflow-glue classify_vcf_svs \\
        --in_vcf ${vcf} \\
        --original ${txt} \\
        --annotated ${annot_txt} \\
        --out_vcf ${vcf.baseName}.annot.vcf
    """
}


// NOTE This is the last touch the VCF has as part of the workflow,
//  we'll rename it with its desired output name here
process sortVCF {
    label "wf_somatic_sv"
    cpus 1
    input:
        tuple val(meta), path(vcf)
    output:
        tuple val(meta), path("${meta.sample}.nanomonsv.result.wf-somatic_sv.vcf.gz"), emit: vcf_gz
        tuple val(meta), path("${meta.sample}.nanomonsv.result.wf-somatic_sv.vcf.gz.tbi"), emit: vcf_tbi
    script:
    """
    bcftools sort -m 2G -O z -o ${meta.sample}.nanomonsv.result.wf-somatic_sv.vcf.gz -T ./ $vcf 
    bcftools index -t ${meta.sample}.nanomonsv.result.wf-somatic_sv.vcf.gz
    """
}

// NOTE This is the last touch the VCF has as part of the workflow,
//  we'll rename it with its desired output name here
// CW-2702: remove sites with END < POS, as in nanomonsv GitHub [issue](https://github.com/friend1ws/nanomonsv/issues/31)
process postprocess_nanomon_vcf {
    label "wf_somatic_sv"
    cpus 1
    input:
        tuple val(meta), path(vcf, stageAs: "input.vcf")
    output:
        tuple val(meta), path("${params.sample_name}.wf-somatic-sv.vcf")
    script:
    def tumor_only = params.bam_normal ? "" : "--tumor_only"
    """
    # Filter out sites with END < POS
    bcftools filter -e "INFO/END < POS" input.vcf > filtered.vcf
    vcf_nanomon2clairs.py --vcf filtered.vcf --sample_id ${params.sample_name} --output "${params.sample_name}.wf-somatic-sv.vcf" $tumor_only
    """
}


process getVersions {
    label "wf_somatic_sv"
    cpus 1
    output:
        path "versions.txt"
    script:
    """
    trap '' PIPE # suppress SIGPIPE without interfering with pipefail
    python -c "import pysam; print(f'pysam,{pysam.__version__}')" >> versions.txt
    nanomonsv --version | head -n 1 | sed 's/ /,/' >> versions.txt
    bcftools --version | head -n 1 | sed 's/ /,/' >> versions.txt
    samtools --version | head -n 1 | sed 's/ /,/' >> versions.txt
    minimap2 --version | head -n 1 | sed 's/^/minimap2,/' >> versions.txt
    """
}


process getParams {
    label "wf_somatic_sv"
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


process report {
    cpus 1
    input:
        tuple val(meta), file(vcf)
        tuple val(meta), file(tbi)
        path eval_json, stageAs: "eval_json/*"
        file versions
        path "params.json"
    output:
        path "*report.html", emit: html
    script:
        def report_name = "${meta.sample}.wf-somatic-sv-report.html"
        // can't use `path.name` here
        // (https://github.com/nextflow-io/nextflow/issues/3574); needs to be
        // `path.fileName.name` instead
        def evalResults = eval_json.fileName.name == 'OPTIONAL_FILE' ? "" : "--eval_results ${eval_json}"
        def tumorOnly = params.bam_normal ? "" : "--tumor_only"
    """
    workflow-glue report_sv \
        $report_name \
        --vcf $vcf \
        --params params.json \
        --params-hidden 'help,schema_ignore_params,${params.schema_ignore_params}' \
        --versions $versions \
        --revision ${workflow.revision} \
        --commit ${workflow.commitId} \
        $evalResults $tumorOnly
    """
}


// See https://github.com/nextflow-io/nextflow/issues/1636
// This is the only way to publish files from a workflow whilst
// decoupling the publish from the process steps.
process output_sv {
    // publish inputs to output directory
    label "wf_somatic_sv"
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
