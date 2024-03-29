/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowDereplicate.initialise(params, log)

// Check input path parameters to see if they exist
def checkPathParamList = [ params.input, params.multiqc_config, params.fasta ]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { INPUT_CHECK_DREP } from '../../subworkflows/local/input_check'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { CUSTOM_DUMPSOFTWAREVERSIONS               } from '../../modules/nf-core/modules/custom/dumpsoftwareversions/main'

//
// MODULE: Local
//
include { PREPARE_DREP                    } from '../../modules/local/drep'
include { DREP as DREP_1                  } from '../../modules/local/drep'
include { DREP as DREP_2                  } from '../../modules/local/drep'
include { DREP as DREP_3                  } from '../../modules/local/drep'
include { DREP as DREP_4                  } from '../../modules/local/drep'
include { DATABASE_CONCAT                 } from '../../modules/local/make_genomedatabase'
include { DATABASE_SCAFFDREP              } from '../../modules/local/make_genomedatabase'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow DEREPLICATE {

    ch_versions = Channel.empty()

    // **********************************
    // Read in samplesheet, validate and stage input files
    // **********************************
    INPUT_CHECK_DREP ()
    ch_drep_sf_raw = INPUT_CHECK_DREP.out.ch_drep_sf
    
    // Group by group
    ch_drep_sf = ch_drep_sf_raw
            .map { meta, drep_sf -> [ meta.group, drep_sf ] }
            .groupTuple()
            .map { group, sf_list -> [group, sf_list.flatten() ]}

    // Prepare dRep for each group
    PREPARE_DREP ( ch_drep_sf,
                params.nondatabase_extrascore,
                params.prepare_drep_params )
    ch_versions = ch_versions.mix(PREPARE_DREP.out.versions)

    // Run dRep for each group
    DREP_1 (    PREPARE_DREP.out.set1,
                params.set1_drep, 
                params.genome_extention)
    DREP_2 (    PREPARE_DREP.out.set2,
                params.set2_drep,
                params.genome_extention)
    DREP_3 (    PREPARE_DREP.out.set3,
                params.set3_drep,
                params.genome_extention)
    DREP_4 (    PREPARE_DREP.out.set4,
                params.set4_drep,
                params.genome_extention)
    ch_versions = ch_versions.mix(DREP_1.out.versions)

    // Create genome database
    ch_gdb = DREP_1.out.genomes.concat(DREP_2.out.genomes, DREP_3.out.genomes, DREP_4.out.genomes)
            .groupTuple()
            .map {group, sets ->
            if( sets.size() == 1) {
                [group, sets[0], [], [], []]
            }
            else if ( sets.size() == 2) {
                [group, sets[0], sets[1], [], []]
            }
            else if ( sets.size() == 3) {
                [group, sets[0], sets[1], [2], []]
            }
            else if ( sets.size() == 4) {
                [group, sets[0], sets[1], [2], [3]]
            }
            else {
                [group]
            }
        }

    if (!params.scaffold_level_drep){
        DATABASE_CONCAT ( ch_gdb,
                        params.genome_extention )
        ch_versions = ch_versions.mix(DATABASE_CONCAT.out.versions)
    } else {
        DATABASE_SCAFFDREP ( ch_gdb,
                        params.genome_extention )
        ch_versions = ch_versions.mix(DATABASE_SCAFFDREP.out.versions)
    }
    
    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
