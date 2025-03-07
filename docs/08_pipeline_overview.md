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

The workflow can run the `snv` component in tumor-only mode.
See [tumor-only mode](#6-tumor-only-mode) for more details.


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
will aggregate the modified bases using [modkit](https://github.com/nanoporetech/modkit). 
The default behaviour of the workflow is to run modkit with the 
`--cpg --combine-strands` options set.
It is possible to report strand-aware modifications by providing `--force_strand`.
Users can further change the behaviour of `modkit` by passing options directly
to modkit via the `--modkit_args` option. This will override any preset,
and allow full control over the run of modkit. For more details on the usage
of `modkit pileup`, checkout the software [documentation](https://nanoporetech.github.io/modkit/).
The workflow will perform differential modification analyses using [DSS](https://bioconductor.org/packages/DSS/)
when the user provides both tumor and normal samples.
DSS is very resource intensive, and might easily run out of memory. Therefore, it is possible to skip this step by setting `--diff_mod false`, saving compute time and allowing the workflow to run to completion.


### 6. Tumor-only mode

It is possible to run a reduced version of the workflow using only the tumor BAM files.
Currently, only the following components can run in tumor-only mode:
- base workflow: BAM coverage and QC statistics
- `--mod`: the workflow will run modkit on the tumor BAM file, but will skip the differentially modified region and loci detection
- `--snv`: the workflow will use [ClairS-TO](https://github.com/HKU-BAL/ClairS-TO), instead of [ClairS](https://github.com/HKU-BAL/ClairS), to call SNVs.


### 7. Run the workflow on a region
When sequencing specific regions or genes, the runtime can vary substantially.
The table below provides a guideline on the average compute time for a given number of genes or regions.
All analyses are run using up to 128GB of RAM and 16 cores, computing `--sv` and `--snv` on a 75X/55X tumor/normal coverage in the regions of interest.

| Number of genes |    Region size    | CPU/h |  Runtime   |
|-----------------|-------------------|-------|------------|
|      1-10       |      200Kb-1Mb    |  ~2.0 |  8m-15m    |
|     11-100      |       1Mb-6Mb     |  ~7.0 | 25m-45m    |
|    100-1000     |      20Mb-60Mb    | ~14.0 | 45m-1h:30m |
