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
| bam_normal | string | Path to a BAM (or CRAM) containing aligned or unaligned reads for the normal sample. | You may choose to provide a BAM/CRAM, but not both. |  |
| bam_tumor | string | Path to a BAM (or CRAM) containing aligned or unaligned reads for the tumor sample. | You may choose to provide a BAM/CRAM, but not both. |  |
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
| classify_insert | boolean | Perform SV insert classification. | Run nanomonsv insert_classify to annotate transposable and repetitive elements for the inserted SV sequences. | False |
| qv | integer | Approximate single base quality value (QV), one of 10, 15, 20 or 25. | Expected single base quality as described in the [nanomonsv web page](https://github.com/friend1ws/nanomonsv#get). |  |
| control_panel | string | Path to the directory containing the non-matched control panel data generated with `nanomonsv merge_control`. | `nanomonsv get` can use a panel of non-matched control data to remove reads carrying common alleles; see [here](https://github.com/friend1ws/nanomonsv#control-panel) for more details. |  |


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
| nanomonsv_get_threads | integer | Total number of threads to use in `nanomonsv get` (minimum of 2 and limited by config executor cpus). |  | 4 |
| dss_threads | integer | Total number of threads to use in the DSS differential modification analysis (limited by config executor cpus). |  | 1 |
| modkit_threads | integer | Total number of threads to use in modkit modified base calling (limited by config executor cpus). |  | 4 |
| haplotype_filter_threads | integer | Set max number of threads to use for the haplotype filtering stage in SNV workflow (limited by config executor cpus). |  | 4 |


