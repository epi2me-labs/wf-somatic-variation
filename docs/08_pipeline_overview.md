This workflow is designed to perform variant calling of small variants, structural
variants and modified bases aggregation from paired tumor/normal BAM files
for a single sample.

Per-sample files will be prefixed with respective aliases and represented
below as {{ alias }}.

### 1. Input and data preparation.

The workflow relies on three primary input files:
1. A reference genome in [fasta format](https://www.ncbi.nlm.nih.gov/genbank/fastaformat/)
2. A [BAM file](https://samtools.github.io/hts-specs/SAMv1.pdf) for the tumor sample (either aligned or unaligned)
3. An optional [BAM file](https://samtools.github.io/hts-specs/SAMv1.pdf) for the normal sample (either aligned or unaligned)

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
(v0.1.6) to identify somatic variants in a paired tumor/normal sample. 
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

SNVs can be annotated using [SnpEff](https://pcingola.github.io/SnpEff/) by
setting `--annotation true`. Furthermore, the workflow will add annotations from
the [ClinVar](https://www.ncbi.nlm.nih.gov/clinvar/) database.


### 4. Somatic structural variant (SV) calling with Nanomonsv.

The workflow allows for the calling of somatic SVs using long-read sequencing data.
Starting from the paired cancer/control samples, the workflow will:
1. Parse the SV signatures in the reads using `nanomonsv parse`
2. Call the somatic SVs using `nanomonsv get`
3. Filter out the SVs in simple repeats using `add_simple_repeat.py` (*optional*)
4. Annotate transposable and repetitive elements using `nanomonsv insert_classify` (*optional*)

As of `nanomonsv` v0.7.1 (and v0.4.0 of this workflow), users can provide
the approximate single base quality value (QV) for their dataset.
To decide which is the most appropriate value for your dataset, visit the
`get` section of the `nanomonsv` [web page](https://github.com/friend1ws/nanomonsv#get),
but it can be summarized as follow:

|     Basecaller     |  Quality value  |
|--------------------|-----------------|
|     guppy (v5)     |       10        |
|  guppy (v5 or v6)  |       15        |
|       dorado       |       20        |

To provide the correct qv value, simply use `--qv 20`.

The VCF produced by nanomonsv is now processed to have one sample with the
name specified with `--sample_name` (rather than the two-sample
`TUMOR`/`CONTROL`). By default, this file does not report a genotype in the resulting
VCF, unless the user specifies the `--genotype_sv` option. In this case, a genotype
is defined based on the number of reference-supporting reads in the tumor sample:
if the value is >= `min_ref_support`, then it is called as heterozygote, otherwise
it is called as homozygote.
The original VCFs generated by nanomonsv can still be accessed by the user, and is
saved as `{{ alias }}/sv/vcf/{{ alias }}.results.nanomonsv.vcf`.

SVs can be annotated using [SnpEff](https://pcingola.github.io/SnpEff/) by
setting `--annotation true`.

It is possible to provide a non-matching control panel with the `--control_panel` option.
`nanomonsv get` can use a panel of non-matched control data to remove reads carrying common
alleles, as described in the [documentation](https://github.com/friend1ws/nanomonsv#control-panel).
These datasets can be generated by using `nanomonsv merge_control` as described [here](https://github.com/friend1ws/nanomonsv#merge_control).


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

