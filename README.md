# Somatic variation workflow
**NOTE: this workflow is currently under active development and still doesn't call somatic variation**
This repository contains a [nextflow](https://www.nextflow.io/) workflow
to identify somatic variation in a paired control/cancer sample.
This workflow is intended to perform:
> Somatic structural variants calling (SV; in development).
> Differentialy methylated regions (DMR; in development).
> Somatic short variant calling (SNP; in development).
> Short tandem repeat detection (STR; in development).





## Introduction

This workflow enables analysis of somatic variation using the following tools:
> [Nanomonsv](https://github.com/friend1ws/nanomonsv) for the somatic SV calling
> [modbam2bed](https://github.com/epi2me-labs/modbam2bed) for the methylated regions calling, and [DSS](https://bioconductor.org/packages/release/bioc/html/DSS.html) for the differentially methylated regions detection
> [ClairS](https://github.com/HKU-BAL/ClairS) to call somatic short variants  




## Quickstart
The workflow uses [nextflow](https://www.nextflow.io/) to manage compute and 
software resources, as such nextflow will need to be installed before attempting
to run the workflow.

The workflow can currently be run using either
[Docker](https://www.docker.com/products/docker-desktop) or
[Singularity](https://docs.sylabs.io/guides/latest/user-guide/) to provide isolation of
the required software. Both methods are automated out-of-the-box provided
either Docker or Singularity is installed.

It is not required to clone or download the git repository in order to run the workflow.
For more information on running EPI2ME Labs workflows [visit our website](https://labs.epi2me.io/wfindex).

**Workflow options**

To obtain the workflow, having installed `nextflow`, users can run:

```
nextflow run epi2me-labs/wf-somatic-variation --help
```

to see the options for the workflow.




## Useful links

* [nextflow](https://www.nextflow.io/)
* [docker](https://www.docker.com/products/docker-desktop)
* [singularity](https://docs.sylabs.io/guides/latest/user-guide/)
