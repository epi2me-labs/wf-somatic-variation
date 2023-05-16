# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [unreleased]
### Fixed
- Demo data
- Fix occasional crash when candidate nested tuple is found
### Added
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
