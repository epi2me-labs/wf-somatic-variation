# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [unreleased]
### Added
- Automated annotation of SNVs and small indels.
    - Disable with `--annotation false`
- ARM-compatible modified base calling

### Changed
- `Input options` and `Output options` have been combined in the `Main options` category
- Updated `DSS` to v2.38.0
- Updated `modkit` to v0.1.12

### Fixed
- The workflow is now mostly species agnostic (with the exception of insert classification enabled with `--classify_insert`)

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

