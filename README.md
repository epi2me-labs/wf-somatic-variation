# Somatic variation workflow
This repository contains a [nextflow](https://www.nextflow.io/) workflow
to identify somatic variation in a paired normal/tumor sample.
This workflow currently perform:
 - Alignment QC and statistics.
 - Somatic short variant calling (SNV and Indels).
 - Somatic structural variants calling (SV).




## Introduction

This workflow enables analysis of somatic variation using the following tools:
1. [ClairS](https://github.com/HKU-BAL/ClairS)
2. [Nanomonsv](https://github.com/friend1ws/nanomonsv)




## Quickstart
The workflow uses [nextflow](https://www.nextflow.io/) to manage compute and 
software resources, as such nextflow will need to be installed before attempting
to run the workflow.

The workflow can currently be run using either
[Docker](https://www.docker.com/products/docker-desktop) or
[Singularity](https://docs.sylabs.io/guides/latest/user-guide/) to provide isolation of
the required software. Both methods are automated out-of-the-box provided
either Docker or Singularity is installed.

It is not required to clone or download the git repository in order to run the workflow.
For more information on running EPI2ME Labs workflows [visit our website](https://labs.epi2me.io/wfindex).


**Workflow options**

To obtain the workflow, having installed `nextflow`, users can run:

```
nextflow run epi2me-labs/wf-somatic-variation --help
```

to see the options for the workflow.

**Somatic short variant calling**

The workflow currently implements a deconstructed version of [ClairS](https://github.com/HKU-BAL/ClairS) (v0.1.0) to identify somatic variants in a paired tumor/normal sample.
This workflow allows to take advantage of the parallel nature of Nextflow, providing the best performance in high-performance, distributed systems.

Currently, ClairS supports the following basecalling models:
 - dna_r10.4.1_e8.2_400bps_sup@v3.5.2
 - dna_r9.4.1_e8_hac@v3.3
 - dna_r9.4.1_e8_sup@v3.3
 - dna_r9.4.1_450bps_hac_prom
 - dna_r9.4.1_450bps_hac
Any other model provided will prevent the workflow to start. 

**Indel calling**

Currently, indel calling is supported only for `dna_r10` basecalling models. When the user specify an r9 model the workflow will automatically skip the indel processes and perform only the SNV calling. 

**Somatic structural variant (SV) calling with Nanomonsv**

The workflow allows for the call of somatic SVs using long-read sequencing data.
Starting from the paired cancer/control samples, the workflow will:
1. Parse the SV signatures in the reads using `nanomonsv parse`
2. Call the somatic SVs using `nanomonsv get`
3. Filter out the SVs in simple repeats using `add_simple_repeat.py` (*optional*)
4. Annotate transposable and repetitive elements using `nanomonsv insert_classify` (*optional*)

**Output folder**

The output directory has the following structure:
```
output/
├── GRCh38_no_alt_chr17.fa
├── GRCh38_no_alt_chr17.fa.fai
├── ref_cache
├── execution # Execution reports
│   ├── report.html
│   ├── timeline.html
│   └── trace.txt
├── qc
│   └── SAMPLE
│       ├── coverage
│       │   ├── SAMPLE_normal.mosdepth.global.dist.txt
│       │   ├── SAMPLE_normal.mosdepth.summary.txt
│       │   ├── SAMPLE_normal.per-base.bed.gz
│       │   ├── SAMPLE_normal.regions.bed.gz
│       │   ├── SAMPLE_normal.thresholds.bed.gz
│       │   ├── SAMPLE_tumor.mosdepth.global.dist.txt
│       │   ├── SAMPLE_tumor.mosdepth.summary.txt
│       │   ├── SAMPLE_tumor.per-base.bed.gz
│       │   ├── SAMPLE_tumor.regions.bed.gz
│       │   └── SAMPLE_tumor.thresholds.bed.gz
│       └── readstats
│           ├── SAMPLE_normal.flagstat.tsv
│           ├── SAMPLE_normal.readstats.tsv.gz
│           ├── SAMPLE_tumor.flagstat.tsv
│           └── SAMPLE_tumor.readstats.tsv.gz
├── snp # ClairS outputs
│   ├── SAMPLE  # ClairS outputs for SAMPLE
│   │   ├── spectra  # Mutational spectra for the workflow; for now, it only works for the SNVs
│   │   │   └── SAMPLE_spectrum.csv
│   │   ├── varstats  # Bcftools stats output
│   │   │   └── SAMPLE.stats
│   │   └── vcf  # VCF outputs
│   │       ├── SAMPLE_somatic_mutype.vcf.gz
│   │       ├── SAMPLE_somatic_mutype.vcf.gz.tbi
│   │       ├── germline  # Clair3 Germline calling for both tumor and normal bams
│   │       │   ├── tumor
│   │       │   │   ├── SAMPLE_tumor_germline.vcf.gz
│   │       │   │   └── SAMPLE_tumor_germline.vcf.gz.tbi
│   │       │   └── normal
│   │       │       ├── SAMPLE_normal_germline.vcf.gz
│   │       │       └── SAMPLE_normal_germline.vcf.gz.tbi
│   │       ├── indels  # VCF containing the indels from ClairS
│   │       │   ├── SAMPLE_somatic_indels.vcf.gz
│   │       │   └── SAMPLE_somatic_indels.vcf.gz.tbi
│   │       └── snv  # VCF containing the SNVs from ClairS
│   │           ├── SAMPLE_somatic_snv.vcf.gz
│   │           └── SAMPLE_somatic_snv.vcf.gz.tbi
│   └── info  # SNV runtime info
│       ├── params.json
│       └── versions.txt
├── sv
│   ├── SAMPLE
│   │   ├── single_breakend
│   │   │   └── SAMPLE.nanomonsv.sbnd.result.txt
│   │   ├── txt
│   │   │   └── SAMPLE.nanomonsv.result.annot.txt
│   │   └── vcf
│   │       ├── SAMPLE.nanomonsv.result.annot.wf_somatic_sv.vcf.gz
│   │       └── SAMPLE.nanomonsv.result.annot.wf_somatic_sv.vcf.gz.tbi
│   └── info  # SV runtime info
│       ├── params.json
│       └── versions.txt
├── SAMPLE.wf-somatic-snp-report.html
├── SAMPLE.wf-somatic-variation-readQC-report.html
├── params.json
└── versions.txt
```
The primary outputs are:
1. `output/snp/SAMPLE/vcf/SAMPLE_somatic_mutype.vcf.gz`: the final VCF file with SNVs and, if r10, InDels
2. `output/snp/SAMPLE/spectra/SAMPLE_spectrum.csv`: the mutation spectrum for the sample
3. `output/snp/SAMPLE/vcf/germline/[tumor/normal]`: the germline calls for both the tumor and normal bam files
4. `output/*.html`: the reports of the SNV pipeline

**Somatic structural variant (SV) calling with Nanomonsv**

The workflow allows for the call of somatic SVs using long-read sequencing data.
Starting from the paired cancer/control samples, the workflow will:
1. Parse the SV signatures in the reads using `nanomonsv parse`
2. Call the somatic SVs using a custom version of `nanomonsv get`
3. Filter out the SVs in simple repeats using `add_simple_repeat.py` (*optional*)
4. Annotate transposable and repetitive elements using `nanomonsv insert_classify` (*optional*)

**Hardware limitations**: the SV calling workflow requires to run on a system supporting AVX2 instructions. please, ensure that 
your system support it before running it.




## Useful links

* [nextflow](https://www.nextflow.io/)
* [docker](https://www.docker.com/products/docker-desktop)
* [singularity](https://docs.sylabs.io/guides/latest/user-guide/)
* [ClairS](https://github.com/HKU-BAL/ClairS)
* [Clair3](https://github.com/HKU-BAL/Clair3)
* [mosdepth](https://github.com/brentp/mosdepth)
* [fastcat](https://github.com/epi2me-labs/fastcat)
* [minimap2](https://github.com/lh3/minimap2)
* [samtools](https://github.com/samtools/samtools)
* [bcftools](https://samtools.github.io/bcftools/bcftools.html)
* [pysam](https://github.com/pysam-developers/pysam)
* [tabix](https://github.com/samtools/htslib)
* [nanomonsv](https://github.com/friend1ws/nanomonsv)
