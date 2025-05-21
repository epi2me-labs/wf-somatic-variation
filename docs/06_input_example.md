<!---Example of input directory structure, delete and edit as appropriate per workflow.--->
This workflow accepts BAM files (aligned or unaligned) as input.

The `--bam_tumor` and `--bam_normal` input parameters for this workflow accept the path to a single BAM file or folder containing multiple BAM files for the tumor sample and the normal sample, respectively. The normal sample is optional. A sample name can be supplied with `--sample_name`.

```
(i)                     (ii)    
input_reads.bam     ─── input_directory
                        ├── reads0.bam
                        └── reads1.bam
```