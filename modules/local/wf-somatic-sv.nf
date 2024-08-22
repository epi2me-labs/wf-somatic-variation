import groovy.json.JsonBuilder
// Memory for the SV process
def severus_mem = [47.GB, 63.GB]

// Severus process
process severus {
    label "wf_somatic_sv"
    cpus Math.max(4, params.severus_threads)
    // Allow retries for testing purposes
    memory { severus_mem[task.attempt - 1] }
    maxRetries 1
    errorStrategy {task.exitStatus in [137,140] ? 'retry' : 'finish'}
    input:
        tuple val(meta),
            path("normal.bam"),
            path("normal.bam.bai"),
            path("tumor.bam"),
            path("tumor.bam.bai")
        tuple path(ref), path(fai), path(ref_cache), env(REF_PATH)
        path tr_bed
        
    output:
        tuple val(meta), path("severus-output/"), emit: all_outputs
        tuple val(meta), path("severus-output/somatic_SVs/severus_somatic.vcf"), emit: vcf

    script:
    // Define tandem repeat file.
    def tr_arg = ""
    if (tr_bed.name != 'OPTIONAL_FILE'){
        tr_arg = "--vntr-bed ${tr_bed}"
    }
    else if (meta.genome_build) {
        log.warn "Automatically selecting TR BED: ${meta.genome_build}.trf.bed"
        tr_arg = "--vntr-bed \${WFSV_TRBED_PATH}/${meta.genome_build}.trf.bed"
    }
    // Define options.
    String options
    if (params.severus_args){
        log.warn "Overriding severus presets with --severus_args: ${params.severus_args}."
        options = "${params.severus_args}"
    }
    else {
        options = "--single-bp --resolve-overlaps --between-junction-ins"
    }
    options += params.vaf_threshold ? " --vaf-thr ${params.vaf_threshold}" : ""
    options += params.min_sv_length ? " --min-sv-size ${params.min_sv_length}" : ""
    options += params.min_support ? " --min-support ${params.min_support}" : ""
    // Severus takes the sample ID for the VCF from the file name.
    // To be consistent with ClairS, rename:
    // - tumor BAM file > {meta.alias}.{bam|cram} > meta.alias becomes the name in the VCF
    // - normal BAM file > {meta.alias}_normal.{bam|cram}
    """
    # Run severus
    severus \
        --target-bam tumor.bam \
        --control-bam normal.bam \
        --out-dir ./severus-output \
        --threads ${task.cpus} \
        ${tr_arg} \
        ${options}
    """
}

// NOTE This is the last touch the VCF has as part of the workflow,
//  we'll rename it with its desired output name here
process sortVCF {
    cpus 2
    memory 4.GB
    input:
        tuple val(meta), path(vcf)
    output:
        tuple val(meta), path("${meta.alias}.wf-somatic-sv.vcf.gz"), emit: vcf_gz
        tuple val(meta), path("${meta.alias}.wf-somatic-sv.vcf.gz.tbi"), emit: vcf_tbi
    script:
    // Severus uses the input file name as sample ID. Fix this using meta.alias here.
    """
    # CW-4313: some INV created by Severus have INFO/END smaller than
    # POS. Exclude these from the final VCF file.
    # https://github.com/epi2me-labs/wf-somatic-variation/issues/28
    echo "tumor\t${meta.alias}" > sample_rename.txt
    bcftools sort -m 2G -O v ${vcf} \
    | bcftools reheader -s sample_rename.txt - \
    | bcftools filter --threads ${task.cpus} -e "INFO/END < POS & INFO/SVTYPE == 'INV'" \
    | bgzip -c > ${meta.alias}.wf-somatic-sv.vcf.gz 
    bcftools index --threads ${task.cpus} -t ${meta.alias}.wf-somatic-sv.vcf.gz
    """
}


process getVersions {
    label "wf_somatic_sv"
    cpus 1
    memory 4.GB
    output:
        path "versions.txt"
    script:
    """
    trap '' PIPE # suppress SIGPIPE without interfering with pipefail
    python -c "import pysam; print(f'pysam,{pysam.__version__}')" >> versions.txt
    severus --version | awk '{OFS=","; print "Severus",\$0}' >> versions.txt
    """
}


process getParams {
    label "wf_somatic_sv"
    cpus 1
    memory 4.GB
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
    label "wf_common"
    cpus 1
    memory 4.GB
    input:
        tuple val(meta), file(vcf)
        tuple val(meta), file(tbi)
        path eval_json, stageAs: "eval_json/*"
        file versions
        path "params.json"
    output:
        tuple val(meta), path("${meta.alias}.wf-somatic-sv-report.html"), emit: html
    script:
        def report_name = "${meta.alias}.wf-somatic-sv-report.html"
        // can't use `path.name` here
        // (https://github.com/nextflow-io/nextflow/issues/3574); needs to be
        // `path.fileName.name` instead
        def evalResults = eval_json.fileName.name == 'OPTIONAL_FILE' ? "" : "--eval_results ${eval_json}"
    """
    workflow-glue report_sv \
        $report_name \
        --vcf $vcf \
        --params params.json \
        --params-hidden 'help,schema_ignore_params,${params.schema_ignore_params}' \
        --versions $versions \
        --workflow_version ${workflow.manifest.version} \
        --revision ${workflow.revision} \
        --commit ${workflow.commitId} \
        $evalResults
    """
}


// See https://github.com/nextflow-io/nextflow/issues/1636
// This is the only way to publish files from a workflow whilst
// decoupling the publish from the process steps.
process publish_sv {
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
