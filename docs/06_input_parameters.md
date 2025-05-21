### Workflow Options

| Nextflow parameter name  | Type | Description | Help | Default |
|--------------------------|------|-------------|------|---------|
| snv | boolean | Call for somatic small variants. | If this option is selected, small variant calling will be carried out using ClairS. | False |
| sv | boolean | Call for somatic structural variants. | If this option is selected, the workflow will call somatic structural variants using severus. | False |
| mod | boolean | Enable output of differentially modified sites and differentially modified regions [requires input BAMs with ML and MM tags]. | If this option is selected, modified bases will be aggregated with modkit and differential modifications will be computed with DSS. | False |


### Main options

| Nextflow parameter name  | Type | Description | Help | Default |
|--------------------------|------|-------------|------|---------|
| sample_name | string | Sample name to be displayed in workflow outputs. | The sample name will be used from the workflow to correctly name output files. | SAMPLE |
| bam_normal | string | BAM or unaligned BAM (uBAM) files for the normal sample to use in the analysis. | This accepts one of two cases: (i) the path to a single BAM file; (ii) the path to a top-level directory containing BAM files. A sample name can be supplied with `--sample_name`. |  |
| bam_tumor | string | BAM or unaligned BAM (uBAM) files for the tumor sample to use in the analysis. | This accepts one of two cases: (i) the path to a single BAM file; (ii) the path to a top-level directory containing BAM files. A sample name can be supplied with `--sample_name`. |  |
| ref | string | Path to a reference FASTA file. | Reference against which to compare reads for variant calling. |  |
| bed | string | An optional BED file enumerating regions to process for variant calling. |  |  |
| tr_bed | string | An optional BED file enumerating simple repeat regions. | This parameter provides a BED file specifying the location of the simple repetitive elements in the genome of choice. This file should be a standard BED file, as described in the [UCSC specification](https://genome.ucsc.edu/FAQ/FAQformat.html#format1). |  |
| pon_file | string | An optional Panel of Normals file enumerating normal expected structural variants. | This parameter provides a TSV file specifying the location of expected structural variations in the genome of choice. [Severus](https://github.com/KolmogorovLab/Severus/blob/main/scripts/pon_from_vcf.py) includes a script to generate these from a VCF file. A default is used for HG38, generated from the 1000 genomes project. |  |
| out_dir | string | Directory for output of all workflow results. |  | output |
| annotation | boolean | Perform SnpEff annotation. | If this option is deselected, VCFs will not be annotated with [SnpEff](https://pcingola.github.io/SnpEff/). | True |
| include_all_ctgs | boolean | Call for variants on all sequences in the reference, otherwise small variants will only be called on chr{1..22,X,Y,MT}. | Enabling this option will call for variants on all contigs of the input reference sequence. Typically this option is not required as standard human reference sequences contain decoy and unplaced contigs that are usually omitted for the purpose of variant calling. This option might be useful for non-standard reference sequence databases. | False |
| igv | boolean | Visualize outputs in the EPI2ME IGV visualizer. | Enabling this option will visualize the output VCF files in the EPI2ME desktop app IGV visualizer. | False |


### Quality Control Options

| Nextflow parameter name  | Type | Description | Help | Default |
|--------------------------|------|-------------|------|---------|
| depth_intervals | boolean | Output a bedGraph file with entries for each genomic interval featuring homogeneous depth. | The output [bedGraph](https://genome.ucsc.edu/goldenPath/help/bedgraph.html) file will have an entry for each genomic interval in which all positions have the same alignment depth. By default this workflow outputs summary depth information from your aligned reads. Per-base depth outputs are slower to generate but may be required for some downstream applications. | False |
| tumor_min_coverage | number | Minimum read coverage for the tumor sample required to run analysis. | A tumor BAM below this coverage value will fail the coverage threshold, and will not be processed in downstream analyses. | 20 |
| normal_min_coverage | number | Minimum read coverage for the normal sample required to run analysis. | A normal BAM below this coverage value will fail the coverage threshold, and will not be processed in downstream analyses. | 20 |


### Small variant calling options

| Nextflow parameter name  | Type | Description | Help | Default |
|--------------------------|------|-------------|------|---------|
| hybrid_mode_vcf | string | Enable hybrid calling mode that combines the *de novo* calling results and genotyping results at the positions in the VCF file given. |  |  |
| genotyping_mode_vcf | string | VCF file input containing candidate sites to be genotyped. Variants will only be called at the sites in the VCF file if provided. |  |  |
| normal_vcf | string | VCF file input containing normal sites for the given sample. | Pointing to a pre-computed normal VCF file will prevent the workflow from calling the germline sites for the normal sample, reducing processing time. |  |
| skip_haplotype_filter | boolean | Skip haplotype filtering of variants. | Setting this will skip haplotype filtering of variants. | False |
| fast_mode | boolean | Fast germline variants calling in Clair3 (does not emit germline calls). | Setting this will speed up the germline calling from Clair3 by relaxing the variant calling parameters; this matches ClairS default behaviour, and therefore will not emit germline VCFs. | False |
| germline | boolean | The workflow will perform germline calling and tumor phasing as default; set to false to disable this (greatly speeds up execution). |  | True |
| phase_normal | boolean | Phase and tag the normal, in addition to the tumor dataset. |  | False |
| liquid_tumor | boolean | The sample is a liquid tumor | Setting this to true will have ClairS to use specific presets and model (where available) for liquid tumors, increasing accuracy. | False |


### Somatic structural variant calling options

| Nextflow parameter name  | Type | Description | Help | Default |
|--------------------------|------|-------------|------|---------|
| min_sv_length | integer | Minimum SV size to call. | Provide the minimum size of the structural variants to call with severus. | 50 |


### Methylation calling options

| Nextflow parameter name  | Type | Description | Help | Default |
|--------------------------|------|-------------|------|---------|
| force_strand | boolean | Require modkit to call strand-aware modifications. |  | False |
| diff_mod | boolean | Detect differentially modified loci and regions with DSS. |  | True |


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


### Advanced Options

| Nextflow parameter name  | Type | Description | Help | Default |
|--------------------------|------|-------------|------|---------|
| override_basecaller_cfg | string | Name of the model to use for converting signal and selecting a small variant calling model. | The workflow will attempt to find the basecaller model from the headers of your input data, providing a value for this option will override the model found in the data. If the model cannot be found in the header, it must be provided with this option as the basecaller model is required for small variant calling. The basecaller model is used to automatically select the appropriate small variant calling model. The model list shows all models that are compatible for small variant calling with this workflow. You should select 'custom' to override the basecaller_cfg with clair3_model_path. |  |


