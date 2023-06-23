#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Argument parsing
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Check mandatory parameters
if (params.genomes) { genomes_input = files(params.genomes) } else { exit 1, 'Input genomes not specified' }


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Import modules
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Import modules
include { GTDBTK_CLASSIFYWF                         } from './modules/nf-core/modules/gtdbtk/classifywf/main'
include { GTDBTK_DB_PREPARATION                     } from './modules/local/gtdbtk_db_preparation'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Do the actual work
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {

    // create meta map
    def meta = [:]
    meta.id         = "GTDB"
    genome_input = [ meta, genomes_input ]

    GTDBTK_DB_PREPARATION ( params.gtdb )

    GTDBTK_CLASSIFYWF ( genome_input,
                        GTDBTK_DB_PREPARATION.out
    )
}