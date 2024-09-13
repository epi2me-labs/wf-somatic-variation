+ If the workflow fails please run it with the demo data set to ensure the workflow itself is working. This will help us determine if the issue is related to the environment, input parameters or a bug.
+ See how to interpret some common nextflow exit codes [here](https://labs.epi2me.io/trouble-shooting/).
+ To run the workflow with a non-human organism, proceed as follows:
    * EPI2ME Desktop Application: disable the `Annotation` option. 
    * Command line: set `--annotation false`.
+ Short somatic Indel calling is supported only for `dna_r10` basecalling models.
+ Renaming, moving or deleting the input BAM, reference genome or the output directory from the location provided at runtime will cause IGV not to load.
+ The workflow expects either an uncompressed or bgzip-compressed reference. If the user provides a reference compressed not with bgzip, the workflow will run to completion, but won't be able to generate the necessary indexes to visualize the outputs in IGV.