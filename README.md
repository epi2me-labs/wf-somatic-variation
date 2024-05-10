# Somatic variation workflow

Nextflow workflow to identify somatic variation.



## Introduction

This workflow calls variants from the alignment files of a paired tumor/normal sample.

This workflow can be used for the following: 

+ Alignment QC and statistics.
+ Somatic short variant calling (SNV and Indels).
+ Somatic structural variants calling (SV).
+ Modified sites calling (mod).




## Compute requirements

Recommended requirements:

+ CPUs = 64
+ Memory = 256GB

Minimum requirements:

+ CPUs = 16
+ Memory = 48GB

Approximate run time: Variable depending on sequencing modality (targeted or whole genome sequencing), as well as coverage and the individual analyses requested. For instance, a complete analysis of a 60X/30X Tumor/Normal pair with default settings takes approximately 6h 30m using the recommended requirements.

ARM processor support: False




## Install and run

These are instructions to install and run the workflow on command line. You can also access the workflow via the [EPI2ME application](https://labs.epi2me.io/downloads/).  

The workflow uses [nextflow](https://www.nextflow.io/) to manage compute and software resources, therefore nextflow will need to be installed before attempting to run the workflow. 

The workflow can currently be run using either [Docker](https://www.docker.com/products/docker-desktop) or
[singularity](https://docs.sylabs.io/guides/3.0/user-guide/index.html) to provide isolation of 
the required software. Both methods are automated out-of-the-box provided 
either docker or singularity is installed. This is controlled by the [`-profile`](https://www.nextflow.io/docs/latest/config.html#config-profiles) parameter as exemplified in the example below. 

It is not required to clone or download the git repository in order to run the workflow. 
More information on running EPI2ME workflows can be found on our [website](https://labs.epi2me.io/wfindex).

The following command can be used to obtain the workflow. This  will pull the repository in to the assets folder of nextflow and provide a list of all parameters available for the workflow as well as an example command:

```
nextflow run epi2me-labs/wf-somatic-variation –help 
```
A demo dataset is provided for testing of the workflow. It can be downloaded using: 
```
wget https://ont-exd-int-s3-euwst1-epi2me-labs.s3.amazonaws.com/wf-somatic-variation/wf-somatic-variation-demo.tar.gz 
tar -xzvf wf-somatic-variation-demo.tar.gz 
```
The workflow can be run with the demo data using: 
```
nextflow run epi2me-labs/wf-somatic-variation \ 
--bam_normal wf-somatic-variation-demo/demo_normal.bam \
--bam_tumor wf-somatic-variation-demo/demo_tumor.bam \
--ref wf-somatic-variation-demo/GCA_000001405.15_GRCh38_no_alt_analysis_set_chr20.fna \
--sv --snv --mod \
--sample_name SAMPLE \
--normal_min_coverage 0 \
--tumor_min_coverage 0 \
-profile standard 
```
For further information about running a workflow on the cmd line see https://labs.epi2me.io/wfquickstart/  



## Related protocols

This workflow is designed to take input sequences that have been produced from [Oxford Nanopore Technologies](https://nanoporetech.com/) devices.

Find related protocols in the [Nanopore community](https://community.nanoporetech.com/docs/).



## Input example

<!---Example of input directory structure, delete and edit as appropriate per workflow.--->
This workflow accepts BAM files (aligned or unaligned) as input.

The `--bam_tumor` and `--bam_normal` input parameters for this workflow accept the path to a single BAM file or folder containing multiple BAM files for the tumor sample and the normal sample, respectively. The normal sample is optional for some components. A sample name can be supplied with `--sample`.

```
(i)                     (ii)    
input_reads.bam     ─── input_directory
                        ├── reads0.bam
                        └── reads1.bam
```



## Input parameters

### Workflow Options

| Nextflow parameter name  | Type | Description | Help | Default |
|--------------------------|------|-------------|------|---------|
| snv | boolean | Call for somatic small variants. | If this option is selected, small variant calling will be carried out using ClairS. | False |
| sv | boolean | Call for somatic structural variants. | If this option is selected, the workflow will call somatic structural variants using nanomonsv. | False |
| mod | boolean | Enable output of differentially modified sites and differentially modified regions [requires input BAMs with ML and MM tags]. | If this option is selected, modified bases will be aggregated with modkit and differential modifications will be computed with DSS. | False |


### Main options

| Nextflow parameter name  | Type | Description | Help | Default |
|--------------------------|------|-------------|------|---------|
| sample_name | string | Sample name to be displayed in workflow outputs. | The sample name will be used from the workflow to correctly name output files. | SAMPLE |
| bam_normal | string | BAM or unaligned BAM (uBAM) files for the normal sample to use in the analysis. | This accepts one of two cases: (i) the path to a single BAM file; (ii) the path to a top-level directory containing BAM files. A sample name can be supplied with `--sample`. |  |
| bam_tumor | string | BAM or unaligned BAM (uBAM) files for the tumor sample to use in the analysis. | This accepts one of two cases: (i) the path to a single BAM file; (ii) the path to a top-level directory containing BAM files. A sample name can be supplied with `--sample`. |  |
| ref | string | Path to a reference FASTA file. | Reference against which to compare reads for variant calling. |  |
| bed | string | An optional BED file enumerating regions to process for variant calling. |  |  |
| tr_bed | string | An optional BED file enumerating simple repeat regions. | This command provides a bed file specifying the location of the simple repetitive elements in the genome of choice. This file should be a standard bed file, as described in the [UCSC specification](https://genome.ucsc.edu/FAQ/FAQformat.html#format1). |  |
| out_dir | string | Directory for output of all workflow results. |  | output |
| annotation | boolean | Perform SnpEff annotation. | If this option is deselected, VCFs will not be annotated with [SnpEff](https://pcingola.github.io/SnpEff/). | True |


### Quality Control Options

| Nextflow parameter name  | Type | Description | Help | Default |
|--------------------------|------|-------------|------|---------|
| depth_intervals | boolean | Output a bedGraph file with entries for each genomic interval featuring homogeneous depth. | The output [bedGraph](https://genome.ucsc.edu/goldenPath/help/bedgraph.html) file will have an entry for each genomic interval in which all positions have the same alignment depth. By default this workflow outputs summary depth information from your aligned reads. Per-base depth outputs are slower to generate but may be required for some downstream applications. | False |
| tumor_min_coverage | number | Minimum read coverage for the tumor sample required to run analysis. | A tumor BAM below this coverage value will fail the coverage threshold, and will not be processed in downstream analyses. | 20 |
| normal_min_coverage | number | Minimum read coverage for the normal sample required to run analysis. | A normal BAM below this coverage value will fail the coverage threshold, and will not be processed in downstream analyses. | 20 |


### Small variant calling options

| Nextflow parameter name  | Type | Description | Help | Default |
|--------------------------|------|-------------|------|---------|
| basecaller_cfg | string | Name of the model to use for converting signal and selecting a small variant calling model. | Required for basecalling and small variant calling. The basecaller configuration is used to automatically select the appropriate small variant calling model. Refer to the [model table on the Dorado repository for selecting a simplex basecalling model](https://github.com/nanoporetech/dorado#available-basecalling-models). | dna_r10.4.1_e8.2_400bps_sup@v3.5.2 |
| hybrid_mode_vcf | string | Enable hybrid calling mode that combines the *de novo* calling results and genotyping results at the positions in the VCF file given. |  |  |
| genotyping_mode_vcf | string | VCF file input containing candidate sites to be genotyped. Variants will only be called at the sites in the VCF file if provided. |  |  |
| normal_vcf | string | VCF file input containing normal sites for the given sample. | Pointing to a pre-computed normal VCF file will prevent the workflow from calling the germline sites for the normal sample, reducing processing time. |  |
| include_all_ctgs | boolean | Call for variants on all sequences in the reference, otherwise small variants will only be called on chr{1..22,X,Y}. | Enabling this option will call for variants on all contigs of the input reference sequence. Typically this option is not required as standard human reference sequences contain decoy and unplaced contigs that are usually omitted for the purpose of variant calling. This option might be useful for non-standard reference sequence databases. | False |
| skip_haplotype_filter | boolean | Skip haplotype filtering of variants. | Setting this will skip haplotype filtering of variants. | False |
| fast_mode | boolean | Fast germline variants calling in Clair3 (does not emit germline calls). | Setting this will speed up the germline calling from Clair3 by relaxing the variant calling parameters; this matches ClairS default behaviour, and therefore will not emit germline VCFs. | False |
| germline | boolean | The workflow will perform germline calling and tumor phasing as default; set to false to disable this (greatly speeds up execution). |  | True |
| phase_normal | boolean | Phase and tag the normal, in addition to the tumor dataset. |  | False |


### Somatic structural variant calling options

| Nextflow parameter name  | Type | Description | Help | Default |
|--------------------------|------|-------------|------|---------|
| min_sv_length | integer | Minimum SV size to call. | Provide the minimum size of the structural variants to call with nanomonsv. | 50 |


### Methylation calling options

| Nextflow parameter name  | Type | Description | Help | Default |
|--------------------------|------|-------------|------|---------|
| force_strand | boolean | Require modkit to call strand-aware modifications. |  | False |


### Multiprocessing Options

| Nextflow parameter name  | Type | Description | Help | Default |
|--------------------------|------|-------------|------|---------|
| ubam_map_threads | integer | Set max number of threads to use for aligning reads from uBAM (limited by config executor cpus). |  | 8 |
| ubam_sort_threads | integer | Set max number of threads to use for sorting and indexing aligned reads from uBAM (limited by config executor cpus). |  | 3 |
| ubam_bam2fq_threads | integer | Set max number of threads to use for uncompressing uBAM and generating FASTQ for alignment (limited by config executor cpus). |  | 1 |
| severus_threads | integer | Total number of threads to use for `Severus` (minimum of 4 and limited by config executor cpus). |  | 8 |
| dss_threads | integer | Total number of threads to use in the DSS differential modification analysis (limited by config executor cpus). |  | 1 |
| modkit_threads | integer | Total number of threads to use in modkit modified base calling (limited by config executor cpus). |  | 4 |
| haplotype_filter_threads | integer | Set max number of threads to use for the haplotype filtering stage in SNV workflow (limited by config executor cpus). |  | 4 |






## Outputs

Output files may be aggregated including information for all samples or provided per sample. Per-sample files will be prefixed with respective aliases and represented below as {{ alias }}.

| Title | File path | Description | Per sample or aggregated |
|-------|-----------|-------------|--------------------------|
| Workflow alignment statistics report | {{ alias }}.wf-somatic-variation-readQC-report.html | Report of the alignment statistics for each tumor/normal paired sample. | per-sample |
| Workflow SNV report | {{ alias }}.wf-somatic-snv-report.html | Report of the SNV for each tumor/normal paired sample. | per-sample |
| Workflow SV report | {{ alias }}.wf-somatic-sv-report.html | Report of the SV for each tumor/normal paired sample. | per-sample |
| Workflow MOD report | {{ alias }}.wf-somatic-mod-report.html | Report of the modified bases for each tumor/normal paired sample. | per-sample |
| Somatic short variant VCF | {{ alias }}.wf-somatic-snv.vcf.gz | VCF file with the somatic SNVs for the sample. | per-sample |
| Somatic short variant VCF index | {{ alias }}.wf-somatic-snv.vcf.gz.tbi | The index of the resulting VCF file with the somatic SNVs for the sample. | per-sample |
| Somatic structural variant VCF | {{ alias }}.wf-somatic-sv.vcf.gz | VCF file with the somatic SVs for the sample. | per-sample |
| Somatic structural variant VCF index | {{ alias }}.wf-somatic-sv.vcf.gz.tbi | The index of the resulting VCF file with the somatic SVs for the sample. | per-sample |
| Modified bases BEDMethyl (normal) | {{ alias }}.wf-somatic-mods.normal.bedmethyl.gz | BED file with the aggregated modification counts for the normal sample. | per-sample |
| Modified bases BEDMethyl (tumor) | {{ alias }}.wf-somatic-mods.tumor.bedmethyl.gz | BED file with the aggregated modification counts for the tumor sample. | per-sample |
| Haplotagged alignment file (normal) | {{ alias }}_normal.ht.bam | BAM or CRAM file with the haplotagged reads for the normal sample. | per-sample |
| Haplotagged alignment file index (normal) | {{ alias }}_normal.ht.bam.bai | The index of the resulting BAM or CRAM file with the haplotagged reads for the normal sample. | per-sample |
| Haplotagged alignment file (tumor) | {{ alias }}_tumor.ht.bam | BAM or CRAM file with the haplotagged reads for the tumor sample. | per-sample |
| Haplotagged alignment file index (tumor) | {{ alias }}_tumor.ht.bam.bai | The index of the resulting BAM or CRAM file with the haplotagged reads for the tumor sample. | per-sample |




## Pipeline overview

This workflow is designed to perform variant calling of small variants, structural
variants and modified bases aggregation from paired tumor/normal BAM files
for a single sample.

Per-sample files will be prefixed with respective aliases and represented
below as {{ alias }}.

### 1. Input and data preparation.

The workflow relies on three primary input files:
1. A reference genome in [fasta format](https://www.ncbi.nlm.nih.gov/genbank/fastaformat/)
2. A single tumor sample in the format of one [BAM file](https://samtools.github.io/hts-specs/SAMv1.pdf), or a folder of BAM files (either aligned or unaligned)
3. An optional single normal sample in the format of one [BAM file](https://samtools.github.io/hts-specs/SAMv1.pdf), or a folder of BAM files (either aligned or unaligned)

The BAM files can be generated from:
1. [POD5](https://github.com/nanoporetech/pod5-file-format)/[FAST5](https://github.com/nanoporetech/ont_fast5_api) files using the [wf-basecalling](https://github.com/epi2me-labs/wf-basecalling) workflow, or
2. [fastq](https://www.ncbi.nlm.nih.gov/sra/docs/submitformats/#fastq) files using [wf-alignment](https://github.com/epi2me-labs/wf-alignment).

Both workflows will generate aligned BAM files that are ready to be used with `wf-somatic-variation`.
It is possible to run the workflow without the BAM file of the "normal" sample. See [tumor-only mode](#6-tumor-only-mode) for more details.

### 2. Data QC and pre-processing.
The workflow starts by performing multiple checks of the input BAM files, as well as computing:
1. The depth of sequencing of each BAM file with [mosdepth](https://github.com/brentp/mosdepth).
2. The read alignment statistics for each BAM file with [fastcat](https://github.com/epi2me-labs/fastcat).

After computing the coverage, the workflow will check that the input BAM files have a depth greater than
`--tumor_min_coverage` and `--normal_min_coverage` for the tumor and normal BAM files, respectively.
It is necessary that **both** BAM files have passed the respective thresholds. In cases where the user
sets the minimum coverage to `0`, the check will be skipped and the workflow will proceed directly to the
downstream analyses.

### 3. Somatic short variants calling with ClairS.

The workflow currently implements a deconstructed version of [ClairS](https://github.com/HKU-BAL/ClairS)
 to identify somatic variants in a paired tumor/normal sample. 
This workflow takes advantage of the parallel nature of Nextflow, providing optimal efficiency in
high-performance, distributed systems.

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

Currently, indel calling is supported only for `dna_r10` basecalling models.
When the user specifies an r9 model the workflow will automatically skip
the indel processes and perform only the SNV calling. 

The workflow uses `Clair3` to call germline sites on both the normal and tumor
sample, which are then used internally to refine the somatic variant calling.
This mode is computationally demanding, and it's behaviour can be changed with
a few options:
* Reduce the accuracy of the variant calling with `--fast_mode`.
* Provide a pre-computed VCF file reporting the germline calls for the normal sample with `--normal_vcf`.
* Disable the germline calling altogether with `--germline false`.

SNVs called using the GRCh37 or GRCh38 genomes can be annotated using [SnpEff](https://pcingola.github.io/SnpEff/)
by setting `--annotation true`. Furthermore, the workflow will add annotations from
the [ClinVar](https://www.ncbi.nlm.nih.gov/clinvar/) database.


### 4. Somatic structural variant (SV) calling with Severus.

The workflow allows for the calling of somatic SVs using long-read sequencing data.
Starting from the paired cancer/control samples, the workflow will call the somatic SVs using [Severus](https://github.com/KolmogorovLab/Severus).

The behaviour of Severus can be tweaked with the options:
1. `--min_sv_length`: minimum size of SVs to call
1. `--min_support`: minimum number of reads to support an SV call
1. `--vaf_threshold`: sites with variant allele frequency (VAF) below this value will be filtered out 
1. `--severus_args`: additional arguments for Severus

SVs called using the GRCh37 or GRCh38 genomes can be annotated using [SnpEff](https://pcingola.github.io/SnpEff/) 
by setting `--annotation true`.


### 5. Modified base calling with modkit

Modified base calling can be performed by specifying `--mod`. The workflow
will aggregate the modified bases using [modkit](https://github.com/nanoporetech/modkit) and
perform differential modification analyses using [DSS](https://bioconductor.org/packages/DSS/). 
The default behaviour of the workflow is to run modkit with the 
`--cpg --combine-strands` options set.
It is possible to report strand-aware modifications by providing `--force_strand`.
Users can further change the behaviour of `modkit` by passing options directly
to modkit via the `--modkit_args` option. This will override any preset,
and allow full control over the run of modkit. For more details on the usage
of `modkit pileup`, checkout the software [documentation](https://nanoporetech.github.io/modkit/).

### 6. Tumor-only mode

It is possible to run a reduced version of the workflow using only the tumor BAM files.
Currently, the following components can run in tumor-only mode:
- base workflow: BAM coverage and QC statistics
- `--mod`: the workflow will run modkit on the tumor BAM file, but will skip the differentially modified region and loci detection

### 7. Run the workflow on a region
When sequencing specific regions or genes, the runtime can vary substantially.
The table below provides a guideline on the average compute time for a given number of genes or regions.
All analyses are run using up to 128GB of RAM and 16 cores, computing `--sv` and `--snv` on a 75X/55X tumor/normal coverage in the regions of interest.

| Number of genes |    Region size    | CPU/h |  Runtime   |
|-----------------|-------------------|-------|------------|
|      1-10       |      200Kb-1Mb    |  ~2.0 |  8m-15m    |
|     11-100      |       1Mb-6Mb     |  ~7.0 | 25m-45m    |
|    100-1000     |      20Mb-60Mb    | ~14.0 | 45m-1h:30m |





## Troubleshooting

+ If the workflow fails please run it with the demo data set to ensure the workflow itself is working. This will help us determine if the issue is related to the environment, input parameters or a bug.
+ See how to interpret some common nextflow exit codes [here](https://labs.epi2me.io/trouble-shooting/).
+ To run the workflow with a non-human organism, proceed as follows:
    * EPI2ME Desktop Application: disable the `Annotation` option. 
    * Command line: set `--annotation false`.
+ Short somatic Indel calling is supported only for `dna_r10` basecalling models.



## FAQ's

+ *Does the workflow calls 5hmC, on top of 5mC?* - Yes, the workflow does call 5hmC, but only if you performed the basecalling with the appropriate module; for more details, check out the [dorado github page](https://github.com/nanoporetech/dorado#dna-models).

+ *Can I run the workflow in tumor-only mode?* - Yes, but currently this mode is only available for `--mod`.




## Related blog posts

+ [Importing third-party workflows into EPI2ME Labs](https://labs.epi2me.io/nexflow-for-epi2melabs/)

See the [EPI2ME website](https://labs.epi2me.io/) for lots of other resources and blog posts.



