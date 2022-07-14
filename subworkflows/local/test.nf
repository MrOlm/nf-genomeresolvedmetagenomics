/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowBin.initialise(params, log)

// TODO nf-core: Add all file path parameters for the pipeline to the list below
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
include { INPUT_CHECK_ASS } from '../../subworkflows/local/input_check'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { CUSTOM_DUMPSOFTWAREVERSIONS               } from '../../modules/nf-core/modules/custom/dumpsoftwareversions/main'
include { METABAT2_METABAT2                         } from '../../modules/nf-core/modules/metabat2/metabat2/main'
include { METABAT2_JGISUMMARIZEBAMCONTIGDEPTHS      } from '../../modules/nf-core/modules/metabat2/jgisummarizebamcontigdepths/main'
include { CHECKM_LINEAGEWF                          } from '../../modules/nf-core/modules/checkm/lineagewf/main'
include { GTDBTK_CLASSIFYWF                         } from '../../modules/nf-core/modules/gtdbtk/classifywf/main'
include { ABRICATE_RUN                              } from '../../modules/nf-core/modules/abricate/run/main'

//
// MODULE: Local
//
include { BOWTIE2_ASSEMBLY_BUILD                    } from '../../modules/local/bowtie2_assembly_build'
include { BOWTIE2_ASSEMBLY_ALIGN                    } from '../../modules/local/bowtie2_assembly_align'
include { COVERM                                    } from '../../modules/local/coverm'
include { GTDBTK_DB_PREPARATION                     } from '../../modules/local/gtdbtk_db_preparation'
include { CHECKVDB_PREPARATION as CHECKV_DB_PREPARATION                     } from '../../modules/local/checkV_db_preparation'
include { CHECKV                                    } from '../../modules/local/checkV'
include { VIRSORTER2                                } from '../../modules/local/virsorter2'
include { VIRSORTER2_DB_PREPARATION                 } from '../../modules/local/virsorter2'
include { EUKREP                                    } from '../../modules/local/eukrep'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow TEST {

    ch_versions = Channel.empty()

    //
    // Read in samplesheet, validate and stage input files
    //
    INPUT_CHECK_ASS ()
    ch_assemblies = INPUT_CHECK_ASS.out.ch_assemblies
    ch_assemblies_raw = ch_assemblies
        .map { meta, ass, reads -> [meta, ass] }

    

    // Run virsorter2
    VIRSORTER2_DB_PREPARATION ()
    VIRSORTER2 (ch_assemblies_raw, 
                VIRSORTER2_DB_PREPARATION.out.database)
    ch_versions = ch_versions.mix(VIRSORTER2.out.versions)
    
    // Run checkv
    // https://github.dev/replikation/What_the_Phage
    CHECKV_DB_PREPARATION ( params.checkv_db )
    CHECKV (ch_assemblies_raw, 
            CHECKV_DB_PREPARATION.out)
    ch_versions = ch_versions.mix(CHECKV.out.versions)
    

    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )

    // // *****************************
    // // Perform eukaryote binning
    // // *****************************
    // EUKREP ( ch_assemblies_raw )
    // ch_versions = ch_versions.mix(EUKREP.out.versions)

    // // ******************************
    // // Perform plasmid identification
    // // ******************************
    // ABRICATE_RUN ( ch_assemblies_raw )
    // ch_versions = ch_versions.mix(ABRICATE_RUN.out.versions)


}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
