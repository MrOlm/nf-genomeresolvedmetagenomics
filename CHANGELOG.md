# nf-core/genomeresolvedmetagenomics: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project (kind of) adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v1.4.5- [1/29/24]
- Fix coverM parsing
- Add "save_bam" argumnet

## v1.4.4- [8/24/23]

- Add coverM parsing
- Don't publish .bam files to results directory

## v1.4.3- [6/23/23]

- Add aux_gtdb_classify.nf

## v1.4.2- [2/21/23]

- Add coverm to the Profile command

## v1.4.1- [2/08/23]

- Update the bowtie2 version for bowtie2_removal_align, which fixed a critical bug
- Other minor changes

## v1.4.0- [9/22/22]

- Various bugfixes and small improvements, especially to Summarize_Binning.py

## v1.3.3- [8/19/22]

- Create Summarize_Binning.py version 2, which keeps phage out of the .stb that are already in bacterial genomes, and reports bins as Euks that are Euks and bacteria 

## v1.3.2- [8/1/22]

- Add "aux_instrain_compare"

## v1.3.1- [8/1/22]

- Bugfix the "DEREPLICATE" subworkflow so it actually works with large datasets

## v1.3.0- [7/15/22]

- Add the "DEREPLICATE" subworkflow

## v1.2.0- [7/12/22]

- Add the "BIN" subworkflow

## v1.1.1- [7/5/22]

- Add a "min_length" parameter and rename scaffolds as part of the "assemble" workflow

## v1.1.0- [6/28/22]

- Add "assemble" workflow

## v1.0.2- [6/24/22]

- Fix a bug with providing genes to inStrain profile
- Add more documentation to the README

## v1.0.1 - [6/24/22]

- Add ability to pass command line arguments to inStrain in PROFILE
- Add basic teting for the pipeline

## v1.0.0 - [6/23/22]

- Initial release of genomeresolvedmetagenomics, created with the [nf-core](https://nf-co.re/) template.
