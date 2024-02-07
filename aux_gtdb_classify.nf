#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Argument parsing
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Validate mandatory parameters
if (!params.genomes) {
    exit 1, 'Input genomes not specified'
}

// Read genome locations from the specified file and group them into chunks of 5000
Channel
    .fromPath(params.genomes)
    .splitCsv(sep: '\n', by: 5000)
    .map { chunk -> 
        chunk.collect { it -> file(it[0]) } 
    }
    .set { genomes_chunks }

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

    // Create meta map
    def meta = [:]
    meta.id = "GTDB"

    GTDBTK_DB_PREPARATION ( params.gtdb )

    // Process each chunk of genomes
    genomes_chunks.map { chunk -> 
        [ meta.clone(), chunk ] 
    }.set { genome_input }

    genome_input.each { chunk ->
        GTDBTK_CLASSIFYWF ( chunk, GTDBTK_DB_PREPARATION.out )
    }
}