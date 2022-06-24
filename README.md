# ![nf-core/genomeresolvedmetagenomics](docs/images/nf-core-genomeresolvedmetagenomics_logo_light.png#gh-light-mode-only) ![nf-core/genomeresolvedmetagenomics](docs/images/nf-core-genomeresolvedmetagenomics_logo_dark.png#gh-dark-mode-only)

## Introduction

**genomeresolvedmetagenomics**  is a [Nextflow DSL2](https://www.nextflow.io/docs/latest/dsl2.html) pipeline written by Matt Olm, primarily for personal use and for use within the Sonnenburg Lab. It was built using the nf-core framework and nf-core suite of tools, but it does not conform to the nf-core overarching style and thus is not available on the nf-core website. I may update this pipeline to conform and upload it there someday

## Pipeline summary

Unlike nf-core pipelines, this pipeline is structured around `entrypoints`. This allows you to easily run some, but not other, pipelines. In pracice each `entrypoint` is a sub-workflow. 

## Quick Start

1. Install [`Nextflow`](https://www.nextflow.io/docs/latest/getstarted.html#installation) (`>=21.10.3`)

2. Install [`Docker`](https://docs.docker.com/engine/installation/) - there are other alternitives that you can use, but I prefer Docker

3. Download the pipeline and test it on a minimal dataset with a single command:

   ```console
   nextflow run nf-core/genomeresolvedmetagenomics -profile test,YOURPROFILE --outdir <OUTDIR>
   ```

4. Start running your own analysis!

   ```console
   nextflow run nf-core/genomeresolvedmetagenomics --input samplesheet.csv --outdir <OUTDIR> --genome GRCh37 -profile <docker/singularity/podman/shifter/charliecloud/conda/institute>
   ```

## Documentation

Access a list of accepted parameters by running:

   ```console
   nextflow run nf-core/genomeresolvedmetagenomics --help
   ```

## Credits

genomeresolvedmetagenomics was written by Matt Olm

## Contributions and Support

If you would like to contribute to this pipeline. please reach out to me via email or by posting an issue on this Github page 

## Citations

While this is not an nf-core pipeline it is built on the back of the nf-core framework. You can cite the `nf-core` publication as follows:

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
