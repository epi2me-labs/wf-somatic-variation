# Somatic variation workflow
**NOTE: this workflow is currently under active development and still doesn't call somatic variation**
This repository contains a [nextflow](https://www.nextflow.io/) workflow
to identify somatic variation in a paired normal/tumor sample.
This workflow is intended to perform:
 - Somatic short variant calling (SNP).




## Introduction

This workflow enables analysis of somatic variation using the following tools:
1. [ClairS](https://github.com/HKU-BAL/ClairS)



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
│   │       ├── germline  # Germline calling for both tumor and normal bams
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
│   ├── info  # Runtime info
│   │   ├── params.json
│   │   └── versions.txt
│   └── reports  # Output report for the workflow
├── SAMPLE.wf-somatic-snp-report.html
├── SAMPLE.wf-somatic-variation-readQC-report.html
├── params.json
└── versions.txt
```
The primary outputs are:
1. `output/snp/SAMPLE/vcf/SAMPLE_somatic_mutype.vcf.gz`: the final VCF file with SNVs and, if r10, InDels
2. `output/snp/SAMPLE/spectra/SAMPLE_spectrum.csv`: the mutation spectrum for the sample
3. `output/snp/SAMPLE/vcf/germline/[tumor/normal]`: the germline calls for both the tumor and normal bam files
4. `output/*.html`: the reports of the SNP pipeline



## Useful links

* [nextflow](https://www.nextflow.io/)
* [docker](https://www.docker.com/products/docker-desktop)
* [singularity](https://docs.sylabs.io/guides/latest/user-guide/)
* [ClairS](https://github.com/HKU-BAL/ClairS)