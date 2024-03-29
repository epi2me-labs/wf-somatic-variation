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
        # use view -H rather than idxstats, as idxstats will still cause a scan of the whole CRAM (https://github.com/samtools/samtools/issues/303)
        samtools view -H ${xam} --no-PG | grep '^@SQ' | sed -nE 's,.*SN:([^[:space:]]*).*LN:([^[:space:]]*).*,\\1\\t\\2,p' > ${xam}_genome.txt
        get_genome.py --chr_counts ${xam}_genome.txt -o output.txt
        genome_build=`cat output.txt`
        """
}
