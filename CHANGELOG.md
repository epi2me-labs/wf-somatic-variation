# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [unreleased]
### Changed
- Initialised wf-somatic-variation from wf-template
- Implemented wf-somatic-snp module, that runs ClairS in a highly parallelised way.
- Implemented some accessory modules to visualise the results from ClairS (mutation spectra, variant allele frequency).  
- Updated to Oxford Nanopore Technologies PLC. Public License
- Implemented report of alignment statistics 

### Fixed
- Workflow crashing when predicting in regions without variants
- Workflow crashing when concatenating SNP and Indel VCF files
- Workflow interrupting when an empty Indel VCF file is generated
