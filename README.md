# genomeresolvedmetagenomics

## Introduction

**genomeresolvedmetagenomics**  is a [Nextflow DSL2](https://www.nextflow.io/docs/latest/dsl2.html) pipeline written by Matt Olm, primarily for personal use. It was built using the nf-core framework and nf-core suite of tools, but it does not conform to the nf-core overarching style and thus is not available on the nf-core website. I may update this pipeline to conform and upload it there someday.

## Pipeline summary

Unlike nf-core pipelines, this pipeline is structured around `entry` points. This allows you to easily run some, but not other, pipelines. In pracice each `entry` is a sub-workflow. Below is a quick summary of the currently available `entry` points:

### PREPROCESSREADS

Takes an input datasheet listing paired reads and runs fastqc to evaluate initial read quality, fastp to trim reads and remove adapters, bowtie2 against the human genome to remove human reads (optional), and a small script at the end to summarize the number of reads and read pairs remaining at the end.

The `--input` parameter must but a .csv file with the columns "sample, "fastq_1", and "fastq_2" (like this - https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/samplesheet/samplesheet_test_illumina_sispa.csv). The locations of the reads can be local, on S3, or on the web

### PROFILE

Takes an input datasheet listing paired reads (you probably want trimmed paired reads, like those you get from **PREPROCESSREADS**!), maps reads to the provides fasta file using bowtie2, and profiles the read mapping using inStrain.

The `--input` parameter has the same rules as above

## Quick Start

1. Install [`Nextflow`](https://www.nextflow.io/docs/latest/getstarted.html#installation) (`>=21.10.3`)

2. Install [`Docker`](https://docs.docker.com/engine/installation/) - there are other alternitives that you can use, but I prefer Docker

3. Test the pipeline on the provided dataset with a single command:

   ```console
   $ nextflow pull https://github.com/MrOlm/nf-genomeresolvedmetagenomics

   $ nextflow run https://github.com/MrOlm/nf-genomeresolvedmetagenomics -profile test,docker -entry PREPROCESSREADS --outdir testout/
   ```

4. Start running your own analysis!

Example run for PREPROCESSREADS:

   ```console
   nextflow run https://github.com/MrOlm/nf-genomeresolvedmetagenomics -entry PREPROCESSREADS --input input_samplesheet_v3.csv -with-report test_report.html --outdir results_1/ -resume
   ```

Example run for PROFILE:

   ```console
   nextflow run https://github.com/MrOlm/nf-genomeresolvedmetagenomics -entry PROFILE --input input_samplesheet_v5.csv -with-report rest_report.html --outdir results_2 -profile docker --fasta transcriptome_chopped_num.fa --stb_file transcriptome.stb -resume
   ```

## Documentation

Access a list of accepted parameters by running:

   ```console
   nextflow run https://github.com/MrOlm/nf-genomeresolvedmetagenomics --help
   ```

## Credits, Contributions, and Support

genomeresolvedmetagenomics was written by Matt Olm

If you would like to contribute to this pipeline please reach out to me via email or by posting an issue or pull request on this Github page 

While this is not an nf-core pipeline it is built on the back of the nf-core framework. You can cite the `nf-core` publication as follows:

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
