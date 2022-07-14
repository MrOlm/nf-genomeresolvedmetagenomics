/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowAssemble.initialise(params, log)

// TODO nf-core: Add all file path parameters for the pipeline to the list below
// Check input path parameters to see if they exist
def checkPathParamList = [ params.input, params.multiqc_config, params.fasta ]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }

def assemblers = params.assemblers ? params.assemblers.split(',').collect{ it.trim().toLowerCase() } : []

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

ch_multiqc_config        = file("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config) : Channel.empty()

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { INPUT_CHECK } from '../../subworkflows/local/input_check'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { CUSTOM_DUMPSOFTWAREVERSIONS               } from '../../modules/nf-core/modules/custom/dumpsoftwareversions/main'
include { MEGAHIT                                   } from '../../modules/nf-core/modules/megahit/main'
include { SPADES                                    } from '../../modules/nf-core/modules/spades/main'
include { QUAST                                     } from '../../modules/nf-core/modules/quast/main'

//
// MODULE: Local
//
include { BOWTIE2_ASSEMBLY_BUILD                    } from '../../modules/local/bowtie2_assembly_build'
include { BOWTIE2_ASSEMBLY_ALIGN                    } from '../../modules/local/bowtie2_assembly_align'
include { COVERM                                    } from '../../modules/local/coverm'
include { ASSEMBLY_FILTER_RENAME                    } from '../../modules/local/assembly_filter_rename'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow ASSEMBLE {

    ch_versions = Channel.empty()

    //
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //
    INPUT_CHECK (
        ch_input
    )
    ch_short_reads = INPUT_CHECK.out.reads
    ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)

    // Co-assembly: prepare grouping for MEGAHIT and for pooling for SPAdes
    if (params.coassemble_group) {
        // short reads
        // group and set group as new id
        ch_short_reads_grouped = ch_short_reads
            .map { meta, reads -> [ meta.group, meta, reads ] }
            .groupTuple(by: 0)
            .map { group, metas, reads ->
                    def meta = [:]
                    meta.id          = "group-$group"
                    meta.group       = group
                    meta.single_end  = false
                    [ meta, reads.collect { it[0] }, reads.collect { it[1] } ]
            }
    } 
    else {
        ch_short_reads_grouped = ch_short_reads
            .map { meta, reads ->
                    [ meta, [reads[0], reads[1]] ] }
    }

    // Run Megahit
    ch_assemblies = Channel.empty()
    if ('megahit' in assemblers) {
        MEGAHIT ( ch_short_reads_grouped )
        ch_megahit_assemblies = MEGAHIT.out.contigs
            .map { meta, assembly ->
                def meta_new = meta.clone()
                meta_new.assembler  = "MEGAHIT"
                [ meta_new, assembly ]
            }
        ch_assemblies = ch_assemblies.mix(ch_megahit_assemblies)
        ch_versions = ch_versions.mix(MEGAHIT.out.versions)
    }

    // Run Metaspades

    // Set up pooled reads for SPAdes
    if (params.coassemble_group) {
        // short reads
        POOL_PAIRED_READS ( ch_short_reads_grouped )
        ch_short_reads_spades = POOL_PAIRED_READS.out.reads
    } else {
        ch_short_reads_spades = ch_short_reads
    }

    if ('spades' in assemblers) {
        SPADES (  
            ch_short_reads_spades.map { meta, fastq -> [ meta, fastq, [], [] ] }, 
            [] 
        )
        ch_spades_assemblies = SPADES.out.scaffolds
            .map { meta, assembly ->
                def meta_new = meta.clone()
                meta_new.assembler  = "SPAdes"
                [ meta_new, assembly ]
            }
        ch_assemblies = ch_assemblies.mix(ch_spades_assemblies)
        ch_versions = ch_versions.mix(SPADES.out.versions)
    }

    // Rename and filter assemblies
    ASSEMBLY_FILTER_RENAME (
        ch_assemblies,
        params.min_length
    )
    ch_assemblies = ASSEMBLY_FILTER_RENAME.out.assemblies
    ch_versions = ch_versions.mix(ASSEMBLY_FILTER_RENAME.out.versions)

    // Calculate assembly completeness
    if (!params.skip_assembly_completeness){
        BOWTIE2_ASSEMBLY_BUILD (
            ch_assemblies
        )

        // Figure out which reads go with which assemblies
        ch_reads_bowtie2 = ch_short_reads.map{ meta, reads -> [ meta.id, meta, reads ] }
        ch_bowtie2_input = BOWTIE2_ASSEMBLY_BUILD.out.assembly_index
            .map { meta, assembly, index -> [ meta.id, meta, assembly, index ] }
            .combine(ch_reads_bowtie2, by: 0)
            .map { id, assembly_meta, assembly, index, reads_meta, reads -> [ assembly_meta, assembly, index, reads_meta, reads ] }
        
        // Align
        BOWTIE2_ASSEMBLY_ALIGN ( ch_bowtie2_input )

        // Run coverM
        COVERM (
            BOWTIE2_ASSEMBLY_ALIGN.out.mappings.map {
                meta, ass, bam, bai -> [meta, bam]
            }
        )
        ch_versions = ch_versions.mix(COVERM.out.versions)

    }

    // Run QUAST
    if (!params.skip_quast){
        QUAST ( 
            ch_assemblies.collect{ it[1] },
            [],
            [],
            false,
            false
        )
        ch_versions = ch_versions.mix(QUAST.out.versions)
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
