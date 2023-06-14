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

**Input and Data preparation**

The workflow relies on three primary input files:
1. A reference genome in [fasta format](https://www.ncbi.nlm.nih.gov/genbank/fastaformat/)
2. An [aligned BAM file](https://samtools.github.io/hts-specs/SAMv1.pdf) for the tumor sample
3. An [aligned BAM file](https://samtools.github.io/hts-specs/SAMv1.pdf) for the normal sample

The reference genome should be either hg19 (GRCh37) or hg38 (GRCh38).
The aligned bam files can be generated starting from:
1. [POD5](https://github.com/nanoporetech/pod5-file-format)/[FAST5](https://github.com/nanoporetech/ont_fast5_api) files using the [wf-basecalling](https://github.com/epi2me-labs/wf-basecalling) workflow, or
2. [fastq](https://www.ncbi.nlm.nih.gov/sra/docs/submitformats/#fastq) files using [wf-alignment](https://github.com/epi2me-labs/wf-alignment).

Both workflows will generate aligned BAM files that are ready to be used with `wf-somatic-variation`.

**Demo data**

The workflow comes with matched demo data accessible [here](https://ont-exd-int-s3-euwst1-epi2me-labs.s3.amazonaws.com/wf-somatic-variation/wf-somatic-variation-demo.tar.gz):
```
wget -q -O demo_data.tar.gz https://ont-exd-int-s3-euwst1-epi2me-labs.s3.amazonaws.com/wf-somatic-variation/wf-somatic-variation-demo.tar.gz
```
This demo is derived from a Tumor/Normal pair of samples, that we have made publicly accessible. Check out our [blog post](https://labs.epi2me.io/colo-2023.05/) for more details.

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
├── execution # Execution reports
│   ├── report.html
│   ├── timeline.html
│   └── trace.txt
├── SAMPLE
│   ├── qc
│   │   ├── coverage
│   │   │   ├── SAMPLE_normal.mosdepth.global.dist.txt
│   │   │   ├── SAMPLE_normal.mosdepth.summary.txt
│   │   │   ├── SAMPLE_normal.per-base.bed.gz
│   │   │   ├── SAMPLE_normal.regions.bed.gz
│   │   │   ├── SAMPLE_normal.thresholds.bed.gz
│   │   │   ├── SAMPLE_tumor.mosdepth.global.dist.txt
│   │   │   ├── SAMPLE_tumor.mosdepth.summary.txt
│   │   │   ├── SAMPLE_tumor.per-base.bed.gz
│   │   │   ├── SAMPLE_tumor.regions.bed.gz
│   │   │   └── SAMPLE_tumor.thresholds.bed.gz
│   │   └── readstats
│   │       ├── SAMPLE_normal.flagstat.tsv
│   │       ├── SAMPLE_normal.readstats.tsv.gz
│   │       ├── SAMPLE_tumor.flagstat.tsv
│   │       └── SAMPLE_tumor.readstats.tsv.gz
│   │
│   ├── snv  # ClairS outputs
│   │   ├── spectra  # Mutational spectra for the workflow; for now, it only works for the SNVs
│   │   │   └── SAMPLE_spectrum.csv
│   │   ├── varstats  # Bcftools stats output
│   │   │   └── SAMPLE.stats
│   │   ├── vcf  # VCF outputs
│   │   │   ├── SAMPLE_tumor_germline.vcf.gz
│   │   │   ├── SAMPLE_tumor_germline.vcf.gz.tbi
│   │   │   ├── SAMPLE_normal_germline.vcf.gz
│   │   │   ├── SAMPLE_normal_germline.vcf.gz.tbi
│   │   │   ├── SAMPLE_somatic_indels.vcf.gz
│   │   │   ├── SAMPLE_somatic_indels.vcf.gz.tbi
│   │   │   ├── SAMPLE_somatic_snv.vcf.gz
│   │   │   └── SAMPLE_somatic_snv.vcf.gz.tbi
│   │   └── info  # SV runtime info
│   │       ├── params.json
│   │       └── versions.txt
│   │
│   └── sv
│       ├── single_breakend
│       │   └── SAMPLE.nanomonsv.sbnd.result.txt
│       ├── txt
│       │   └── SAMPLE.nanomonsv.result.annot.txt
│       └── info  # SV runtime info
│           ├── params.json
│           └── versions.txt
├── SAMPLE_somatic_mutype.vcf.gz
├── SAMPLE_somatic_mutype.vcf.gz.tbi
├── SAMPLE.nanomonsv.result.wf_somatic_sv.vcf.gz
├── SAMPLE.nanomonsv.result.wf_somatic_sv.vcf.gz.tbi
├── SAMPLE.wf-somatic-snp-report.html
├── SAMPLE.wf-somatic-sv-report.html
├── SAMPLE.wf-somatic-variation-readQC-report.html
├── params.json
└── versions.txt
```
The primary outputs are:
1. `output/SAMPLE_somatic_mutype.vcf.gz`: the final VCF file with SNVs and, if r10, InDels
2. `output/SAMPLE.nanomonsv.result.wf_somatic_sv.vcf.gz`: the final VCF with the somatic SVs from nanomonsv
3. `output/*.html`: the reports of the different stages
4. `output/SAMPLE/snp/spectra/SAMPLE_spectrum.csv`: the mutation spectrum for the sample
5. `output/SAMPLE/snp/vcf/germline/[tumor/normal]`: the germline calls for both the tumor and normal bam files
6. `output/SAMPLE/sv/txt/SAMPLE.nanomonsv.result.annot.txt`: the somatic SVs called with nanomonsv in tabular format
7. `output/SAMPLE/sv/single_breakend/SAMPLE.nanomonsv.sbnd.result.txt`: the single break-end SVs called with nanomonsv

**Somatic structural variant (SV) calling with Nanomonsv**

The workflow allows for the call of somatic SVs using long-read sequencing data.
Starting from the paired cancer/control samples, the workflow will:
1. Parse the SV signatures in the reads using `nanomonsv parse`
2. Call the somatic SVs using a custom version of `nanomonsv get`
3. Filter out the SVs in simple repeats using `add_simple_repeat.py` (*optional*)
4. Annotate transposable and repetitive elements using `nanomonsv insert_classify` (*optional*)

**Hardware limitations**: the SV calling workflow requires to run on a system supporting AVX2 instructions. please, ensure that 
your system support it before running it.
