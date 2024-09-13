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
| Modified bases BEDMethyl (normal) | {{ alias }}.wf-somatic-mod.normal.bedmethyl.gz | BED file with the aggregated modification counts for the normal sample. | per-sample |
| Modified bases summary (normal) | {{ alias }}.normal.mod_summary.tsv | Summary modification stats for the normal sample. | per-sample |
| Single-change BEDMethyl (normal) | {{ alias }}/mod/{{ mod type }}/bedMethyl/{{ mod type }}.{{ alias }}.wf-somatic-mod.normal.bedmethyl.gz | BED file with the aggregated modification counts for a single modification type (e.g. 5mc) for the normal sample. | per-sample |
| Single-change DSS input file (normal) | {{ alias }}/mod/{{ mod type }}/DSS/{{ mod type }}.{{ alias }}_normal.dss.tsv | Input text file for DSS for a single modification type (e.g. 5mc) for the normal sample. | per-sample |
| Modified bases BEDMethyl (tumor) | {{ alias }}.wf-somatic-mod.tumor.bedmethyl.gz | BED file with the aggregated modification counts for the tumor sample. | per-sample |
| Modified bases summary (tumor) | {{ alias }}.normal.mod_summary.tsv | Summary modification stats for the tumor sample. | per-sample |
| Single-change BEDMethyl (tumor) | {{ alias }}/mod/{{ mod type }}/bedMethyl/{{ mod type }}.{{ alias }}.wf-somatic-mod.tumor.bedmethyl.gz | BED file with the aggregated modification counts for a single modification type (e.g. 5mc) for the tumor sample. | per-sample |
| Single-change DSS input file (tumor) | {{ alias }}/mod/{{ mod type }}/DSS/{{ mod type }}.{{ alias }}_tumor.dss.tsv | Input text file for DSS for a single modification type (e.g. 5mc) for the tumor sample. | per-sample |
| Differentially modified loci (DML) per change type | {{ alias }}/mod/{{ mod type }}/DML/{{ alias }}.{{ mod type }}.dml.tsv | Differentially modified loci from DSS for a single modification type (e.g. 5mc). | per-sample |
| Differentially modified regions (DMR) per change type | {{ alias }}/mod/{{ mod type }}/DMR/{{ alias }}.{{ mod type }}.dmr.tsv | Differentially modified regions from DSS for a single modification type (e.g. 5mc). | per-sample |
| Alignment file (normal) | {{ alias }}/bam/normal/reads.bam | BAM or CRAM file with the aligned reads for the normal sample. | per-sample |
| Alignment file index (normal) | {{ alias }}/bam/normal/reads.bam.bai | The index of the resulting BAM or CRAM file with the aligned reads for the normal sample. | per-sample |
| Alignment file (tumor) | {{ alias }}/bam/tumor/reads.bam | BAM or CRAM file with the aligned reads for the tumor sample. | per-sample |
| Alignment file index (tumor) | {{ alias }}/bam/tumor/reads.bam.bai | The index of the resulting BAM or CRAM file with the aligned reads for the tumor sample. | per-sample |
| Haplotagged alignment file (normal) | {{ alias }}_normal.ht.bam | BAM or CRAM file with the haplotagged reads for the normal sample. | per-sample |
| Haplotagged alignment file index (normal) | {{ alias }}_normal.ht.bam.bai | The index of the resulting BAM or CRAM file with the haplotagged reads for the normal sample. | per-sample |
| Haplotagged alignment file (tumor) | {{ alias }}_tumor.ht.bam | BAM or CRAM file with the haplotagged reads for the tumor sample. | per-sample |
| Haplotagged alignment file index (tumor) | {{ alias }}_tumor.ht.bam.bai | The index of the resulting BAM or CRAM file with the haplotagged reads for the tumor sample. | per-sample |
