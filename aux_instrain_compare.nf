#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Argument parsing
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Set defaults
params.stb = []
params.is_compare_args = ""
if (params.stb){
    stb = file(params.stb)
}
else {
    stb = params.stb
}

// Check input path parameters to see if they exist
def checkPathParamList = [ params.input, params.stb ]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified' }


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Import modules
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Import modules
include { INPUT_CHECK_INSTRAIN               } from './subworkflows/local/input_check'
include { INSTRAINCOMPARE                    } from './modules/local/instrain_compare'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Do the actual work
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {
    ch_versions = Channel.empty() 

    // Load input table
    INPUT_CHECK_INSTRAIN()

    // Group by group
    ch_IS = INPUT_CHECK_INSTRAIN.out.ch_IS
            .map { meta, IS, bam -> [ meta.group, IS, bam ] }
            .groupTuple()
            .map { group, IS_list, bam_list -> [group, IS_list.flatten(), bam_list.flatten() ]}

    // Run
    INSTRAINCOMPARE(
        ch_IS,
        stb,
    )
}