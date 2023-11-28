+ If the workflow fails please run it with the demo data set to ensure the workflow itself is working. This will help us determine if the issue is related to the environment, input parameters or a bug.
+ See how to interpret some common nextflow exit codes [here](https://labs.epi2me.io/trouble-shooting/).
+ To run the workflow with a non-human organism, proceed as follows:
    * EPI2ME app: disable the `Classify insertion sequence` and `Annotation` options. 
    * Command line: set `--classify_insert false` and `--annotation false`.
+ The SV calling workflow requires a computing system supporting AVX2 instructions. Please, ensure that your system supports these before running it.
+ Short somatic Indel calling is supported only for `dna_r10` basecalling models.