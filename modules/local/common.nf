import groovy.json.JsonBuilder

process cram_cache {
    cpus 1
    memory 4.GB
    input:
        path reference
    output:
        path("ref_cache/"), emit: ref_cache
        env(REF_PATH), emit: ref_path
    shell:
    '''
    # Invoke from binary installed to container PATH
    seq_cache_populate.pl -root ref_cache/ !{reference}
    REF_PATH="ref_cache/%2s/%2s/%s"
    '''
}

process index_ref_fai {
    cpus 1
    memory 4.GB
    input:
        file reference
    output:
        path "${reference}.fai", emit: reference_index
    """
    samtools faidx ${reference}
    """
}

process index_ref_gzi {
    cpus 1
    memory 4.GB
    input:
        file reference
    output:
        path "${reference}.gzi", emit: reference_index
    """
    bgzip -r ${reference}
    """
}

// NOTE -f required to compress symlink
process decompress_ref {
    cpus 1
    memory 4.GB
    input:
        file compressed_ref
    output:
        path "${compressed_ref.baseName}", emit: decompressed_ref
    """
    gzip -df ${compressed_ref}
    """
}

process minimap2_ubam {
    cpus {params.ubam_map_threads + params.ubam_sort_threads + params.ubam_bam2fq_threads}
    memory { 15.GB * task.attempt }
    maxRetries 2
    errorStrategy {task.exitStatus in [137,140] ? 'retry' : 'finish'}

    input:
        path reference
        path old_reference
        tuple path(reads), path(reads_idx)
    output:
        tuple path("${params.sample_name}.cram"), path("${params.sample_name}.cram.crai"), emit: alignment
    script:
    def bam2fq_ref = old_reference.name != "OPTIONAL_FILE" ? "--reference ${old_reference}" : ''
    // samtools sort need extra threads. Specify these as N - 1 to avoid using too many cores.
    def sort_threads = params.ubam_sort_threads - 1
    """
    samtools bam2fq -@ ${params.ubam_bam2fq_threads} -T 1 ${bam2fq_ref} ${reads} | minimap2 -y -t ${params.ubam_map_threads} -ax map-ont ${reference} - \
    | samtools sort -@ ${sort_threads} --write-index -o ${params.sample_name}.cram##idx##${params.sample_name}.cram.crai -O CRAM --reference ${reference} -
    """
}

// Module to convert fai index to bed
process bgzipper {
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
    cpus 2
    memory 7.GB
    input:
        tuple val(meta), path("input.vcf.gz"), path("input.vcf.gz.tbi")
        val(output_label)
    output:
        tuple val(meta), path("${meta.sample}.wf-${output_label}.vcf.gz"), path("${meta.sample}.wf-${output_label}.vcf.gz.tbi"), emit: annot_vcf
        tuple val(meta), path("${meta.sample}.wf-${output_label}-snpEff-genes.txt"), emit: gene_txt
        tuple val(meta), path("${meta.sample}.wf-${output_label}-clinvar.vcf"), emit: annot_vcf_clinvar
    shell:
    '''
    # deal with samples which aren't hg19 or hg38
    if [[ "!{meta.genome_build}" != "hg38" ]] && [[ "!{meta.genome_build}" != "hg19" ]]; then
        # return the original VCF and index as the outputs
        cp input.vcf.gz !{params.sample_name}.wf-!{output_label}.vcf.gz
        cp input.vcf.gz.tbi !{params.sample_name}.wf-!{output_label}.vcf.gz.tbi
        # create an empty snpEff_genes file
        touch !{params.sample_name}.wf-!{output_label}-snpEff-genes.txt
        # create an empty ClinVar VCF
        touch !{params.sample_name}.wf-!{output_label}-clinvar.vcf
    else
        # do some annotation
        if [[ "!{meta.genome_build}" == "hg38" ]]; then
            snpeff_db="GRCh38.p13"
            clinvar_vcf="${CLINVAR_PATH}/clinvar_GRCh38.vcf.gz"
            
        elif [[ "!{meta.genome_build}" == "hg19" ]]; then
            snpeff_db="GRCh37.p13"
            clinvar_vcf="${CLINVAR_PATH}/clinvar_GRCh37.vcf.gz"
        fi

        # Specify N-1 G of memory, otherwise SnpEff will crash with the default 1G
        snpEff -Xmx!{task.memory.giga - 1}g ann $snpeff_db input.vcf.gz > !{params.sample_name}.snpeff_annotated.vcf
        # Add ClinVar annotations
        SnpSift annotate $clinvar_vcf !{params.sample_name}.snpeff_annotated.vcf > !{params.sample_name}.wf-!{output_label}.vcf
        # Get the ClinVar-annotated variants into a separate VCF
        cat !{params.sample_name}.wf-!{output_label}.vcf | SnpSift filter "( exists CLNSIG )" > !{params.sample_name}.wf-!{output_label}-clinvar.vcf
    
        bgzip -c !{params.sample_name}.wf-!{output_label}.vcf > !{params.sample_name}.wf-!{output_label}.vcf.gz
        tabix !{params.sample_name}.wf-!{output_label}.vcf.gz
    
        # tidy up
        rm !{params.sample_name}.snpeff_annotated.vcf
        rm !{params.sample_name}.wf-!{output_label}.vcf
        mv snpEff_genes.txt !{params.sample_name}.wf-!{output_label}-snpEff-genes.txt
    fi
    '''
}


process getVersions {
    cpus 1
    memory 4.GB
    output:
        path "versions.txt"
    script:
    """
    python --version | tr -s ' ' ',' | tr '[:upper:]' '[:lower:]' > versions.txt
    minimap2 --version | sed 's/^/minimap2,/' >> versions.txt
    samtools --version | (head -n 1 && exit 0) | sed 's/ /,/' >> versions.txt
    fastcat --version | sed 's/^/fastcat,/' >> versions.txt
    mosdepth --version | sed 's/ /,/' >> versions.txt
    ezcharts --version | sed 's/ /,/' >> versions.txt
    python -c "import pysam; print(pysam.__version__)" | sed 's/^/pysam,/'  >> versions.txt
    bgzip --version | awk 'NR==1 {print \$1","\$3}' >> versions.txt
    tabix --version | awk 'NR==1 {print \$1","\$3}' >> versions.txt
    """
}


process getParams {
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
    cpus 1
    memory 4.GB
    input:
        tuple path(xam), path(xam_idx), val(xam_meta)
    output:
        tuple path(xam), path(xam_idx), val(xam_meta), env(genome_build), emit: genome_build, optional: true
     script:
        """
        samtools idxstats ${xam} > ${xam}_genome.txt
        get_genome.py --chr_counts ${xam}_genome.txt -o output.txt
        genome_build=`cat output.txt`
        """
}
