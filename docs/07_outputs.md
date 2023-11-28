Outputs files may be aggregated including information for all             samples or provided per sample. Per sample files             will be prefixed with respective aliases and represented             below as {{ alias }}.

| Title | File path | Description | Per sample or aggregated |
|-------|-----------|-------------|--------------------------|
| Workflow alignment statistics report | ./{{ alias }}.wf-somatic-variation-readQC-report.html | Report of the alignment statistics for each tumor/normal paired sample. | per-sample |
| Workflow SNV report | ./{{ alias }}.wf-somatic-snv-report.html | Report of the SNV for each tumor/normal paired sample. | per-sample |
| Workflow SV report | ./{{ alias }}.wf-somatic-sv-report.html | Report of the SV for each tumor/normal paired sample. | per-sample |
| Workflow MOD report | ./{{ alias }}.wf-somatic-mod-report.html | Report of the modified bases for each tumor/normal paired sample. | per-sample |
| Somatic short variant VCF | ./{{ alias }}.wf-somatic-snv.vcf.gz | VCF file with the somatic SNVs for the sample. | per-sample |
| Somatic short variant VCF index | ./{{ alias }}.wf-somatic-snv.vcf.gz.tbi | The index of the resulting VCF file with the somatic SNVs for the sample. | per-sample |
| Somatic structural variant VCF | ./{{ alias }}.wf-somatic-sv.vcf.gz | VCF file with the somatic SVs for the sample. | per-sample |
| Somatic structural variant VCF index | ./{{ alias }}.wf-somatic-sv.vcf.gz.tbi | The index of the resulting VCF file with the somatic SVs for the sample. | per-sample |
| Modified bases BEDMethyl | ./{{ alias }}_{{ type }}.wf_mod.bedmethyl.gz | BED file with the aggregated modification counts for the tumor or normal sample. | per-sample |
| Haplotagged alignment file | ./{{ alias }}_{{ type }}.ht.{{ format }} | BAM or CRAM file with the haplotagged reads for the tumor or normal sample. | per-sample |
| Haplotagged alignment file index | ./{{ alias }}_{{ type }}.ht.{{ format }}.{{ format_index }} | The index of the resulting BAM or CRAM file with the haplotagged reads for the tumor or normal sample. | per-sample |
