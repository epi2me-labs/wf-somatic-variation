import groovy.json.JsonBuilder

// NOTE -f required to compress symlink
process decompress {
    label "wf_common"
    cpus 1
    memory 4.GB
    input:
        file compressed
    output:
        path "${compressed.baseName}", emit: decompressed
    """
    gzip -df ${compressed}
    """
}

// Module to convert fai index to bed
process bgzipper {
    label "wf_common"
    cpus 1
    memory 4.GB
    input:
        path infile
    output:
        path "${infile.baseName}.gz", emit: bgzip
        path "${infile.baseName}.gz.tbi", emit: tbi
    script:
    """
    sort -k 1,1 -k2,2n ${infile} | bgzip -c > ${infile.baseName}.gz
    """
}

// Module to convert fai index to bed
process tabixer {
    label "wf_common"
    cpus 1
    memory 4.GB
    input:
        path infile
    output:
        path "${infile.baseName}.tbi", emit: tbi
    script:
    // Define ideal preset based on the suffix
    def preset = ""
    if (infile.extension == 'bed'){
        preset = "-p bed"
    } else if (infile.extension == 'vcf'){
        preset = "-p vcf"
    } else if (infile.extension == 'gff' || infile.extension == 'gff3' || infile.extension == 'gtf'){
        preset = "-p gff"
    } else {
        preset = "-s 1 -b 2"
    }
    """
    tabix ${preset} ${infile.baseName}.gz
    """
}

// Module to convert fai index to bed
process getAllChromosomesBed {
    label "wf_common"
    cpus 1
    memory 4.GB
    input:
        tuple path(reference), path(ref_idx), path(ref_cache), env(REF_PATH)
    output:
        path "allChromosomes.bed", emit: all_chromosomes_bed
    """
    awk '{OFS="\t"; print \$1, "0", \$2}' ${ref_idx} > allChromosomes.bed
    # faidx --transform bed $reference > allChromosomes.bed
    """
}


process annotate_vcf {
    // use SnpEff to generate basic functional annotations, SnpSift annotate to add 
    // ClinVar annotations, and SnpSift filter to produce a separate VCF of Clinvar-annotated 
    // variants - if any variants are present in this file, it is used to populate a table in 
    // the report.
    label "snpeff_annotation"
    cpus 1
    memory 7.GB
    input:
        tuple val(meta), path("input.vcf.gz"), path("input.vcf.gz.tbi"), val(contig)
        val(output_label)
    output:
        tuple val(meta), path("${meta.sample}.wf-${output_label}*.vcf.gz"), path("${meta.sample}.wf-${output_label}*.vcf.gz.tbi"), emit: annot_vcf
    shell:
    '''
    if [ "!{contig}" == '*' ]; then
        # SV is quick to annotate, dont bother breaking it apart
        INPUT_FILENAME=input.vcf.gz
        OUTPUT_LABEL="!{output_label}"
    else
        # SNP is slow to annotate, we'll break it apart by contig
        # and merge it back later. filter the master VCF to current contig
        bcftools view -r !{contig} input.vcf.gz | bgzip > input.chr.vcf.gz
        INPUT_FILENAME=input.chr.vcf.gz
        OUTPUT_LABEL="!{output_label}.!{contig}"
    fi

    # deal with samples which aren't hg19 or hg38
    if [[ "!{meta.genome_build}" != "hg38" ]] && [[ "!{meta.genome_build}" != "hg19" ]]; then
        # return the original VCF and index as the outputs
        cp ${INPUT_FILENAME} !{meta.sample}.wf-${OUTPUT_LABEL}.vcf.gz
        cp ${INPUT_FILENAME}.tbi !{meta.sample}.wf-${OUTPUT_LABEL}.vcf.gz.tbi
    else
        # do some annotation
        if [[ "!{meta.genome_build}" == "hg38" ]]; then
            snpeff_db="GRCh38.p13"
            clinvar_vcf="${CLINVAR_PATH}/clinvar_GRCh38.vcf.gz"

        elif [[ "!{meta.genome_build}" == "hg19" ]]; then
            snpeff_db="GRCh37.p13"
            clinvar_vcf="${CLINVAR_PATH}/clinvar_GRCh37.vcf.gz"
        fi

        snpEff -Xmx!{task.memory.giga - 1}g ann -noStats -noLog $snpeff_db ${INPUT_FILENAME} > !{meta.sample}.intermediate.snpeff_annotated.vcf
        # Add ClinVar annotations
        SnpSift annotate $clinvar_vcf !{meta.sample}.intermediate.snpeff_annotated.vcf | bgzip > !{meta.sample}.wf-${OUTPUT_LABEL}.vcf.gz
        tabix !{meta.sample}.wf-${OUTPUT_LABEL}.vcf.gz

        # tidy up
        rm !{meta.sample}.intermediate*
    fi
    '''
}

process sift_clinvar_vcf {
    label "snpeff_annotation"
    cpus 1
    memory 3.GB
    input:
        tuple val(meta), path("input.vcf.gz"), path("input.vcf.gz.tbi")
        val(output_label)
    output:
        tuple val(meta), path("${meta.sample}.wf-${output_label}_clinvar.vcf"), emit: final_vcf_clinvar
    shell:
    '''
    # deal with samples which aren't hg19 or hg38
    if [[ "!{meta.genome_build}" != "hg38" ]] && [[ "!{meta.genome_build}" != "hg19" ]]; then
        # create an empty ClinVar VCF
        touch !{meta.sample}.wf-!{output_label}_clinvar.vcf
    else
        bcftools view input.vcf.gz | SnpSift filter "( exists CLNSIG )" > !{meta.sample}.wf-!{output_label}_clinvar.vcf
    fi
    '''
}

process concat_vcfs {
    cpus 2
    memory 3.GB
    input:
        tuple val(meta), path("vcfs/*"), path("vcfs/*")
        val(prefix)
    output:
        tuple val(meta), path("${prefix}.vcf.gz"), path("${prefix}.vcf.gz.tbi"), emit: final_vcf
    script:
        def concat_threads = Math.max(task.cpus - 1, 1)
        """
        bcftools concat --threads ${concat_threads} -O u vcfs/*.vcf.gz | bcftools sort -O z - > ${prefix}.vcf.gz
        tabix -p vcf ${prefix}.vcf.gz
        """
}

process getVersions {
    label "wf_common"
    cpus 1
    memory 4.GB
    output:
        path "versions.txt"
    script:
    """
    python --version | tr -s ' ' ',' | tr '[:upper:]' '[:lower:]' > versions.txt
    samtools --version | (head -n 1 && exit 0) | sed 's/ /,/' >> versions.txt
    fastcat --version | sed 's/^/fastcat,/' >> versions.txt
    ezcharts --version | sed 's/ /,/' >> versions.txt
    python -c "import pysam; print(pysam.__version__)" | sed 's/^/pysam,/'  >> versions.txt
    bgzip --version | awk 'NR==1 {print \$1","\$3}' >> versions.txt
    tabix --version | awk 'NR==1 {print \$1","\$3}' >> versions.txt
    """
}

process getVersions_somvar {
    cpus 1
    memory 4.GB
    input:
        path "versions.tmp.txt"
    output:
        path "versions.txt"
    script:
    """
    cat versions.tmp.txt > versions.txt
    minimap2 --version | sed 's/^/minimap2,/' >> versions.txt
    mosdepth --version | sed 's/ /,/' >> versions.txt
    bedtools --version | sed 's/ /,/' >> versions.txt
    """
}


process getParams {
    label "wf_common"
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


process getGenome {
    label "wf_common"
    cpus 1
    memory 4.GB
    input:
        tuple path(xam), path(xam_idx), val(xam_meta)
    output:
        tuple path(xam), path(xam_idx), val(xam_meta), env(genome_build), emit: genome_build, optional: true
     script:
        """
        # use view -H rather than idxstats, as idxstats will still cause a scan of the whole CRAM (https://github.com/samtools/samtools/issues/303)
        samtools view -H ${xam} --no-PG | grep '^@SQ' | sed -nE 's,.*SN:([^[:space:]]*).*LN:([^[:space:]]*).*,\\1\\t\\2,p' > ${xam}_genome.txt
        get_genome.py --chr_counts ${xam}_genome.txt -o output.txt
        genome_build=`cat output.txt`
        """
}

// Generate report file.
process report {
    label "wf_common"
    input:
        tuple val(meta), 
            path(reports),
            path('versions.txt'),
            path('params.json')
    output:
        path "${params.sample_name}.wf-somatic-variation-report.html", emit: html
    script:
        // Define report name.
        def report_name = "${params.sample_name}.wf-somatic-variation-report.html"
        def output_dir = file(params.out_dir)
        """
        workflow-glue report \\
            $report_name \\
            --outdir_path ${output_dir} \\
            --sample_id ${params.sample_name} \\
            --versions versions.txt \\
            --params params.json \\
            --workflow_version ${workflow.manifest.version}
        """
}