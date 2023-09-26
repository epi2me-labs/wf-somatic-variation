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


**Hardware limitations**: the SV calling workflow requires to run on a system supporting AVX2 instructions. Please, ensure that 
your system supports these before running it.


**Input and data preparation**

The workflow relies on three primary input files:
1. A reference genome in [fasta format](https://www.ncbi.nlm.nih.gov/genbank/fastaformat/)
2. A [BAM file](https://samtools.github.io/hts-specs/SAMv1.pdf) for the tumor sample (either aligned or unaligned)
3. A [BAM file](https://samtools.github.io/hts-specs/SAMv1.pdf) for the normal sample (either aligned or unaligned)

The BAM files can be generated from:
1. [POD5](https://github.com/nanoporetech/pod5-file-format)/[FAST5](https://github.com/nanoporetech/ont_fast5_api) files using the [wf-basecalling](https://github.com/epi2me-labs/wf-basecalling) workflow, or
2. [fastq](https://www.ncbi.nlm.nih.gov/sra/docs/submitformats/#fastq) files using [wf-alignment](https://github.com/epi2me-labs/wf-alignment).
Both workflows will generate aligned BAM files that are ready to be used with `wf-somatic-variation`.


**Working with non-human samples**
The workflow is designed to work with human samples, and in particular with either the hg19 (GRCh37) or hg38 (GRCh38) reference genomes.
Nevertheless, almost all tasks carried out by the workflow are species agnostic, and will work with non-human and non-model species just as well. 
To run the workflow with a non-human organism, the user will need to disable some options as follows:
1. For EPI2ME users: ensure that the `Classify insertion sequence` and `Annotation` options are de-selected. 
2. For CLI users: add the `--classify_insert false --annotation false` options to your command line.


**Demo data**

The workflow comes with matched demo data accessible [here](https://ont-exd-int-s3-euwst1-epi2me-labs.s3.amazonaws.com/wf-somatic-variation/wf-somatic-variation-demo.tar.gz):
```
wget -q -O demo_data.tar.gz https://ont-exd-int-s3-euwst1-epi2me-labs.s3.amazonaws.com/wf-somatic-variation/wf-somatic-variation-demo.tar.gz
```
This demo is derived from a Tumor/Normal pair of samples, that we have made publicly accessible. Check out our [blog post](https://labs.epi2me.io/colo-2023.05/) for more details.


**Somatic short variant calling**

The workflow currently implements a deconstructed version of [ClairS](https://github.com/HKU-BAL/ClairS) (v0.1.5) to identify somatic variants in a paired tumor/normal sample.
This workflow allows to take advantage of the parallel nature of Nextflow, providing the best performance in high-performance, distributed systems.

Currently, ClairS supports the following basecalling models:

| Workflow basecalling model | ClairS model |
|----------------------------|--------------|
| dna_r10.4.1_e8.2_400bps_sup@v4.2.0 | ont_r10_dorado_5khz |
| dna_r10.4.1_e8.2_400bps_sup@v4.1.0 | ont_r10_dorado_4khz |
| dna_r10.4.1_e8.2_400bps_sup@v3.5.2 | ont_r10_guppy |
| dna_r9.4.1_e8_hac@v3.3 | ont_r9_guppy |
| dna_r9.4.1_e8_sup@v3.3 | ont_r9_guppy |
| dna_r9.4.1_450bps_hac_prom | ont_r9_guppy |
| dna_r9.4.1_450bps_hac | ont_r9_guppy |

Any other model provided will prevent the workflow to start. 

**Indel calling**

Currently, indel calling is supported only for `dna_r10` basecalling models. When the user specifies an r9 model the workflow will automatically skip the indel processes and perform only the SNV calling. 

**Tweaking the variant calling**

By default, the workflow uses `Clair3` to call germline sites on both the normal and tumor sample, that are then used internally to refine the somatic variant calling.
This mode is computationally demanding, and it's behaviour can be changed with a few options:
* Reduce the accuracy of the variant calling with `--fast_mode`; 
* Provide a pre-computed VCF file reporting the germline calls for the normal sample with `--normal_vcf`;
* Disable the germline calling altogether with `--germline false`.


**Somatic structural variant (SV) calling with Nanomonsv**

The workflow allows for the calling of somatic SVs using long-read sequencing data.
Starting from the paired cancer/control samples, the workflow will:
1. Parse the SV signatures in the reads using `nanomonsv parse`
2. Call the somatic SVs using `nanomonsv get`
3. Filter out the SVs in simple repeats using `add_simple_repeat.py` (*optional*)
4. Annotate transposable and repetitive elements using `nanomonsv insert_classify` (*optional*)

As of `nanomonsv` v0.7.1 (and v0.4.0 of the workflow), users can provide the approximate single base quality value (QV) for their dataset. To decide which is the most appropriate value for your dataset, visit the `nanomonsv get` [web page](https://github.com/friend1ws/nanomonsv#get), but it can be summarized as follow:

|     Basecaller     |  Quality value  |
|--------------------|-----------------|
|     guppy (v5)     |       10        |
|  guppy (v5 or v6)  |       15        |
|       dorado       |       20        |

To provide the correct qv value, simply use `--qv 20`.

**08/08/2023**: the vcf produced by nanomonSV is now processed to have one sample with the name specified with `--sample_name` (rather than the two-sample `TUMOR`/`CONTROL`). The original VCFs generated by nanomonSV are now saved in:
```
├── SAMPLE
│   ├── sv
│   │   └── vcf  # VCF outputs
│   │       ├── SAMPLE.results.nanomonsv.vcf

```
See below for the full updated output directory structure.

**Modified base calling**

Modified base calling can be performed by specifying `--mod`. The workflow will call modified bases using [modkit](https://github.com/nanoporetech/modkit). 
The default behaviour of the workflow is to run modkit with the `--cpg --combine-strands` options set. It is possible to report strand-aware modifications 
by providing `--force_strand`, which will trigger modkit to run in default mode.
The modkit run can be fully customized by providing `--modkit_args`. This will override any preset, and allow full control over the run of modkit.


**Output folder**

The primary outputs are:
1. `output/SAMPLE.wf-somatic-snv.vcf.gz`: the final VCF file with SNVs and, if r10, InDels
2. `output/SAMPLE.wf-somatic-sv.vcf.gz`: the final VCF with the somatic SVs from nanomonsv
3. `output/SAMPLE.[normal/tumor].mod_summary.tsv`: the modification summary file for the [tumor/normal] sample
4. `output/*.html`: the reports of the different stages
5. `output/*.ht.cram`: haplotagged CRAM files for the tumor (and normal when `--phase_normal` is provided) samples

Additional outputs include:
1. QC subfiles are saved in `output/SAMPLE/qc/`, and include:
    * `coverage/`: additional statistic files created by `mosdepth`
    * `readstats/`: additional statistic files created by `bamstats`
2. SNV specific subfiles are saved in `output/SAMPLE/snv/`, and include: 
    * `change_counts/`: the change counts for the sample
    * `annot/`: additional annotation files, such as the gene table and ClinVar vcf file
    * `varstats/`: variant statistics computed by `bcftools stats`
    * `gvcf/`: germline GVCF files
    * `vcf/`: more VCF files, such as germline, raw SNV and Indels calls 
3. SV specific subfiles are saved in `output/SAMPLE/sv/`, and include: 
    * `vcf/`: the raw somatic SVs called with nanomonsv in VCF format
    * `annot/`: additional annotation files, such as the gene table and ClinVar vcf file
    * `single_breakend/`: the single break-end SVs called with nanomonsv
    * `txt/`: the somatic SVs called with nanomonsv in tabular format
4. Modified-bases specific subfiles are saved in `output/SAMPLE/mod/`, and include:
    * `raw/`: the raw bedMethyl files from `modkit`
    * `[CHANGE]/bedMethyl`: bedMethyl files for change type `CHANGE`
    * `[CHANGE]/DSS/`: DSS inputs for change type `CHANGE`
    * `[CHANGE]/DML`: DSS differentially modified loci output for change type `CHANGE`
    * `[CHANGE]/DMR`: DSS differentially modified regions output for change type `CHANGE`

`CHANGE` refers to the change type analysed (e.g. modC, 5mC, 5hmC, etc) and is automatically detected by `modkit` from the input
`modbam`.

A full-analysis output directory should present a folder structure like the following one:
```
output/
├── execution # Execution reports
│   ├── report.html
│   ├── timeline.html
│   └── trace.txt
│
├── SAMPLE
│   ├── qc
│   │   ├── coverage
│   │   │   ├── filtered.bed
│   │   │   ├── SAMPLE_normal.mosdepth.global.dist.txt
│   │   │   ├── SAMPLE_normal.mosdepth.summary.txt
│   │   │   ├── SAMPLE_normal.per-base.bed.gz
│   │   │   ├── SAMPLE_normal.regions.bed.gz
│   │   │   ├── SAMPLE_normal.regions.filt.bed.gz
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
│   │   ├── annot  # Annotion files
│   │   │   ├── SAMPLE.wf-snp-snpEff-genes.txt
│   │   │   └── SAMPLE.wf-snp-clinvar.vcf
│   │   ├── change_counts  # Mutational change counts for the sample; for now, it only works for the SNVs
│   │   │   └── SAMPLE_changes.csv
│   │   ├── varstats  # Bcftools stats output
│   │   │   └── SAMPLE.stats
│   │   ├── gvcf  # GVCF outputs
│   │   │   ├── SAMPLE_tumor_germline.gvcf.gz
│   │   │   ├── SAMPLE_tumor_germline.gvcf.gz.tbi
│   │   │   ├── SAMPLE_normal_germline.gvcf.gz
│   │   │   └── SAMPLE_normal_germline.gvcf.gz.tbi
│   │   └── vcf  # VCF outputs
│   │       ├── SAMPLE_tumor_germline.vcf.gz
│   │       ├── SAMPLE_tumor_germline.vcf.gz.tbi
│   │       ├── SAMPLE_normal_germline.vcf.gz
│   │       ├── SAMPLE_normal_germline.vcf.gz.tbi
│   │       ├── SAMPLE_somatic_indels.vcf.gz
│   │       ├── SAMPLE_somatic_indels.vcf.gz.tbi
│   │       ├── SAMPLE_somatic_snv.vcf.gz
│   │       └── SAMPLE_somatic_snv.vcf.gz.tbi
│   │
│   ├── sv
│   │   ├── annot  # Annotion files
│   │   │   ├── SAMPLE.wf-sv-snpEff-genes.txt
│   │   │   └── SAMPLE.wf-sv-clinvar.vcf
│   │   ├── vcf  # Raw nanomonsv VCF
│   │   │   └── SAMPLE.nanomonsv.result.vcf
│   │   ├── single_breakend
│   │   │   └── SAMPLE.nanomonsv.sbnd.result.txt
│   │   └── txt
│   │       └── SAMPLE.nanomonsv.result.annot.txt
│   │
│   └── mod
│       ├── 5mC   # Modified bases code
│       │   ├── bedMethyl   # bedMethyl output files
│       │   │   ├── 5mC.SAMPLE_normal.bed.gz
│       │   │   └── 5mC.SAMPLE_tumor.bed.gz
│       │   ├── DML   # Differentially methylated loci
│       │   │   └── SAMPLE.5mC.dml.tsv
│       │   ├── DMR   # Differentially methylated regions
│       │   │   └── SAMPLE.5mC.dmr.tsv
│       │   └── DSS   # DSS input files
│       │       ├── 5mC.SAMPLE_normal.dss.tsv
│       │       └── 5mC.SAMPLE_tumor.dss.tsv
│       └── raw   # Raw outputs from modkit
│           ├── SAMPLE_normal.bed
│           └── SAMPLE_tumor.bed
│
├── info  # single component runtime info
│   ├── mod
│   │   ├── params.json
│   │   └── versions.txt
│   ├── snv
│   │   ├── params.json
│   │   └── versions.txt
│   └── sv
│       ├── params.json
│       └── versions.txt
│
├── SAMPLE.wf-somatic-snv.vcf.gz
├── SAMPLE.wf-somatic-snv.vcf.gz.tbi
├── SAMPLE.wf-somatic-sv.vcf.gz
├── SAMPLE.wf-somatic-sv.vcf.gz.tbi
├── SAMPLE.normal.mod_summary.tsv
├── SAMPLE.tumor.mod_summary.tsv
├── SAMPLE.wf-somatic-snp-report.html
├── SAMPLE.wf-somatic-sv-report.html
├── SAMPLE.wf-somatic-mod-report.html
├── SAMPLE.wf-somatic-variation-readQC-report.html
├── SAMPLE_normal.ht.cram
├── SAMPLE_normal.ht.cram.crai
├── SAMPLE_tumor.ht.cram
├── SAMPLE_tumor.ht.cram.crai
├── params.json
└── versions.txt
```
