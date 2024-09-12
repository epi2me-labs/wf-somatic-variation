# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [unreleased]
### Added
- BED files for the VNTR regions in hg19/hg38 from [Severus](https://github.com/KolmogorovLab/Severus/tree/main/vntrs).
    - The appropriate file will be automatically selected for the appropriate genome, unless a user provides a custom bed with `--tr_bed`.

### Changed
- Tweaked parameters for Severus to refine SV calling.
- More informative log when a BAM called with an invalid basecaller model is provided.
- Failures of processes involved in differentially modified loci and regions detection will not cause workflow to fail.
- Updated `modkit` to v0.3.3.

### Fixed
- `-resume` failing for some `snv` processes.
- Excessive memory usage for sample_probs process when using --mod leading to exit code 137.
- Erroneous handling of inputs for the joint report.

## [v1.3.1]
### Added
- `--override_basecaller_cfg` parameter allows users to provide a basecall configuration name in cases where automatic basecall model detection fails.
- `--diff_mod` option to allow users to turn off differential modified loci (DMR) and regions (DMR) analysis with [DSS](https://bioconductor.org/packages/DSS) by setting `--diff_mod false`.

### Changed
- Updated `modkit` to v0.3.0.
- Reconciled workflow `_ingress.nf` from wf-human-variation v2.3.1.

### Fixed
- `--snv` crashing with `--include_all_ctgs true`.
- `--snv` always referring to model in `--override_basecaller_cfg` when deciding whether to call Indels.
- Automated basecaller detection not finding a basecaller model.

### Removed
- `--basecaller_cfg` as the workflow now automatically detects the basecaller model from the input data.

## [v1.3.0]
### Added
- Tumor-only mode for the base workflow, small variant calling with ClairS-TO, modified base aggregation and somatic SV calling.

### Fixed
- Parts of the documentation still referring to nanomonsv.

## Changed
- If available `basecaller_cfg` will be inferred from the `basecall_model` DS key of input read groups.
    - Providing `--basecaller_cfg` will not be required if `basecall_model` is present in the DS tag of the read groups of the input BAM.
    - `basecaller_cfg` will be ignored if a `basecall_model` is found in the input BAM.
    - The workflow will fail if the tumor and normal BAM files have not been called with the same `basecall_model`.
- Updated to Severus v1.1.

### Fixed
- Workflow crashing when the input BED file has overlapping intervals.
- Returning error in `annotate_sv` when `END` position smaller than `POS`.
- Workflow not starting in nextflow v24.04.

## [v1.2.2]
### Changed
- Update to ClairS v0.2.0.
    - The workflow now uses normal heterozygote sites to haplotag reads.
    - This behaviour can be changed using `--use_normal_hets_for_phasing` and `--use_tumor_hets_for_phasing`.
    - Moreover, it uses the indels for the phasing as default; this can be changed with `--use_het_indels_for_phasing`.
    - The workflow now uses longphase to haplotag reads; this can be changed with `--use_longphase_haplotag`.
    - The workflow now accepts a `--liquid_tumor` option, enabling presets and, where available, models specific for liquid tumors.
- Update to Clair3 v1.0.8.
    - Added `--clair3_base_err` and `--clair3_gq_bin_size` options.
- `modkit` now runs by contig.
- `modkit` bedMethyl are now in the top level output directory.

## [v1.2.1]
### Added
- A report with name `[sample name].wf-somatic-variation-report.html`, linking the individual detailed reports.

### Changed
- Use `ezcharts SeqCompare` to in QC report.
- Memory usage of alignment report reduced by using histograms.
- Retry process when `clairs.py predict` crashes with error 134.

## [v1.2.0]
### Added
- Support for input folders of BAM files for `--bam_tumor` and `--bam_normal` (instead of only allowing single BAM files).

### Changed
- ClinVar version in SnpEff container updated to version 20240307
- Update to Clair3 v1.0.7.
- Update to modkit v0.2.6.
- Improved modkit runtime by increasing the default interval size.
- Increased minimum CPU requirement for the workflow to 16.
- bedMethyl output files now follow the pattern `{{ alias }}.wf-somatic-mods.{{ type }}.bedmethyl.gz`.
- Structural variant (SV) calling is now performed with [Severus](https://github.com/KolmogorovLab/Severus) (v0.1.2).
- ARM-compatible base workflow and modified base calling.
- `minimap2` alignments will be in BAM format when `--sv` is set.

### Fixed
- Force minimap2 to clean up memory more aggressively. Empirically this reduces peak-memory use over the course of execution.
- Workflow occasionally repeating QC analyses when resuming, even if successful.
- Alignment report script using too much memory.
- HAC models not recognised as valid.
- Spurious `bamstats` crashes with implausible alignment information.

### Removed
- CRAM as supported input format.
- Reference genome and its indexes from the output directory.
- Insert classification and support for mismatching panel of control.
- Options `--qv`, `--classify_insert`, `--min_ref_support`, `--genotype_sv` and `--control_panel`, as these are no longer used by the workflow.

## [v1.1.0]
### Changed
- Updated ClairS to v0.1.7, with the new dorado 4KHz/5KHz HAC models.
- Several performance improvements which should noticeably reduce the running time of the workflow
- `makeQCreport` allows for one retry to prevent the workflow failing on report generation

## [v1.0.0]
### Added
- Tumor-only mode for the base workflow and modified base aggregation.
- `--control_panel` option to provide a non-matching control panel produced with [nanomonsv merge_control](https://github.com/friend1ws/nanomonsv#control-panel).
- Memory directives for every process.

### Changed
- Run `ClairS` `haplotype_filter` by contig.
- `--sv` will no longer emit genotypes by default, and only save the same fields from `nanomonsv`.
    - Users can still request a genotype using `--genotype_sv`.
    - `--min_ref_support` defines the minimum number of REF-supporting reads to call a heterozygote site.

### Fixed
- `--snv` calling genome-wide variants when `--bed` is specified.
- `snv:makeReport` crashing when `--annotation false`.

### Removed
- `--annotation_threads` removed as `snpEff` v5 uses only one thread.

## [v0.5.2]
### Changed
- New documentation
- Updated ClairS to v0.1.6
- ClinVar annotation of SVs has been temporarily removed due to not being correctly incorporated. SnpEff annotations are still produced as part of the final SV VCF.

### Fixed
- rVersion retried also upon success

### Removed
- Default local executor CPU and RAM limits

## [v0.5.1]
### Added
- List of reads supporting the SV events is now emitted in `{params.output}/{params.sample_name}/sv/txt`

### Fixed
- Running in `--sv` mode does not resume properly.
- `somatic_sv:report` process failing due to name collisions when running with `--annotation false`.
- Modifed base calling report showing overlapping lines in the DMR plot.

## [v0.5.0]
### Added
- Automated annotation of SNVs, small indels and SVs.
- Option to skip germline calling and variant phasing with `--germline false`.
- Workflow can now emit germline GVCFs for tumor/normal samples with `--GVCF`.
- Option `--normal_vcf` to provide a pre-computed normal VCF file.
- Add genotyping and hybrid mode.
- The workflow saves the haplotagged cram files if `--germline true`.

### Changed
- Updated `modkit` to v0.1.13, `ClairS` to v0.1.5, `Clair3` to v1.0.4, and added support for 4KHz and 5KHz dorado models.
- Runs `ClairS` `haplotype_filter` on the indels in addition to the SNVs.
- SV VCF now report a single sample with Tumor/Normal formats collected, rather than two distinct samples
    - This facilitates merging multiple samples thanks to unambiguous sample labelling
    - The original VCFs generated by nanomonSV are now saved in `${outdir}/${sample_name}/sv/vcf`
- More consistent output naming formatting.
- When `--bed` and `--<tumor|normal>_min_coverage` are specified, the workflow will process the regions with coverage above the threshold
    - The filtering will retain the union (`bedtools merge -d 0`) of the tumor and normal filtered regions

### Fixed
- snpEff crashing due to running out of memory.
- Error in `annotate_sv` when `END` position smaller than `POS`

## [v0.4.0]
### Added
- Automated annotation of SNVs and small indels.
    - Disable with `--annotation false`
- ARM-compatible base workflow and modified base calling
- Option `--qv`, to specify the expected quality value for nanomonsv

### Changed
- `Input options` and `Output options` have been combined in the `Main options` category
- Updated `DSS` to v2.38.0
- Updated `modkit` to v0.1.12
- Updated `nanomonsv` to v0.7.1

### Fixed
- The workflow is now mostly species agnostic (with the exception of `--annotation` and `--classify_insert`)

## [v0.3.0]
### Added
- Add modified base calling with `--mod`

### Changed
- Updated to nanomonsv v0.6.0.
- Refine change type counts and plot; file named [SAMPLE]_spectrum.csv renamed to [SAMPLE]_changes.csv.

## [v0.2.0]
### Added
- Coverage plots to alignment stats report

### Changed
- Improved documentation and new structure of the output sub-directories
- Variant allele frequency representation is now a scatterplot showing the relationship between normal and tumor VAF
- Created a subworkflow to call somatic SV (`somatic_sv` in `workflows/wf-somatic-sv.nf`)
- Add nanomonsv soft filtering SV when providing a bed file specifying the tandem repeat with `--tr_bed`
- Add nanomonsv insert classification with `--classify_insert`, to add RepeatMasker annotation to the SVs
- Add `report_sv`, that generates a report of the SV detected
- Enum choices are enumerated in the `--help` output
- Enum choices are enumerated as part of the error message when a user has selected an invalid choice
- Bumped minimum required Nextflow version to 22.10.8

### Fixed
- Replaced `--threads` option in fastqingress with hardcoded values to remove warning about undefined `param.threads`
- Depth plot with overlapping chromosomal coverage

## [v0.1.1]
### Fixed
- Demo data
- Fix occasional crash when candidate nested tuple is found

### Added
- Configuration for running demo data in AWS
- Fast mode

### Changed
- Updated to ClairS v0.1.1

## [v0.1.0]
### Changed
- Initialised wf-somatic-variation from wf-template
- Implemented wf-somatic-snp module, that runs ClairS in a highly parallelised way (v0.1.0).
- Implemented some accessory modules to visualise the results from ClairS (mutation counts, variant allele frequency).
- Implemented report of alignment statistics
- Update reporting script to use ezcharts
- Customizable thread count for haplotype filtering stage
- Mosdepth per-base statistics are now optional, and are not shown in the report

### Fixed
- Workflow crashing when predicting in regions without variants
- Workflow crashing when concatenating SNP and Indel VCF files
- Workflow interrupting when an empty Indel VCF file is generated
- Extremely slow reporting of alignment statistics

