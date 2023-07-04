import groovy.json.JsonBuilder

process cram_cache {
    input:
        path reference
    output:
        path("ref_cache/"), emit: cram_cache
    script:
    def is_conda = workflow.profile.toLowerCase().contains("conda")
    """
    if [[ "${is_conda}" == "true" ]]; then
        wget https://raw.githubusercontent.com/samtools/samtools/master/misc/seq_cache_populate.pl;
        # Invoke without messing with PATH and +x
        INVOCATION='perl seq_cache_populate.pl'
    else
        # Invoke from binary installed to container PATH
        INVOCATION='seq_cache_populate.pl'
    fi
    \$INVOCATION -root ref_cache/ ${reference}
    """
}

process index_ref_fai {
    cpus 1
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
    input:
        path reference
        path old_reference
        tuple path(reads), path(reads_idx)
    output:
        tuple path("${params.sample_name}.cram"), path("${params.sample_name}.cram.crai"), emit: alignment
    script:
    def bam2fq_ref = old_reference.name != "OPTIONAL_FILE" ? "--reference ${old_reference}" : ''
    """
    samtools bam2fq -@ ${params.ubam_bam2fq_threads} -T 1 ${bam2fq_ref} ${reads} | minimap2 -y -t ${params.ubam_map_threads} -ax map-ont ${reference} - \
    | samtools sort -@ ${params.ubam_sort_threads} --write-index -o ${params.sample_name}.cram##idx##${params.sample_name}.cram.crai -O CRAM --reference ${reference} -
    """
}

// Module to convert fai index to bed
process bgzipper {
    cpus 1
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
    input:
        tuple path(reference), path(ref_idx), path(ref_cache)
    output:
        path "allChromosomes.bed", emit: all_chromosomes_bed
    """
    awk '{OFS="\t"; print \$1, "0", \$2}' ${ref_idx} > allChromosomes.bed
    # faidx --transform bed $reference > allChromosomes.bed
    """
}

process getVersions {
    cpus 1
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
