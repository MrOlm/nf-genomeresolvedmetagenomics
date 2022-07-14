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
include { PRODIGAL                                  } from '../../modules/nf-core/modules/prodigal/main'
include { DIAMOND_BLASTP                            } from '../../modules/nf-core/modules/diamond/blastp/main'

//
// MODULE: Local
//
include { BOWTIE2_ASSEMBLY_BUILD                    } from '../../modules/local/bowtie2_assembly_build'
include { BOWTIE2_ASSEMBLY_ALIGN                    } from '../../modules/local/bowtie2_assembly_align'
include { COVERM                                    } from '../../modules/local/coverm'
include { GTDBTK_DB_PREPARATION                     } from '../../modules/local/gtdbtk_db_preparation'
include { CHECKVDB_PREPARATION as CHECKV_DB_PREP    } from '../../modules/local/checkV_db_preparation'
include { CHECKV                                    } from '../../modules/local/checkV'
include { VIRSORTER2                                } from '../../modules/local/virsorter2'
include { VIRSORTER2_DB_PREPARATION                 } from '../../modules/local/virsorter2'
include { EUKREP                                    } from '../../modules/local/eukrep'
include { TREP                                      } from '../../modules/local/trep'
include { PLASFLOW                                  } from '../../modules/local/plasflow'
include { VIRFINDER                                 } from '../../modules/local/virfinder'
include { PARSESTB                                  } from '../../modules/local/parse_stb'
include { SUMMARY_BINNING                           } from '../../modules/local/summarize_binning'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow BIN {

    ch_versions = Channel.empty()

    // **********************************
    // Read in samplesheet, validate and stage input files
    // **********************************
    INPUT_CHECK_ASS ()
    ch_assemblies = INPUT_CHECK_ASS.out.ch_assemblies
    ch_assemblies_raw = ch_assemblies
        .map { meta, ass, reads -> [meta, ass] }

    // Run Prodigal
    PRODIGAL ( ch_assemblies_raw,
             "gff" )

    // *************************
    // Perform bacterial binning
    // *************************

    // Make bowtie2 indicies
    BOWTIE2_ASSEMBLY_BUILD ( ch_assemblies_raw )
    ch_versions = ch_versions.mix(BOWTIE2_ASSEMBLY_BUILD.out.versions)

    // Group samples for co-mapping
    ch_reads_bowtie2 = ch_assemblies.map{ meta, ass, reads -> [ meta.group, meta, reads ] }
    if (!params.skip_crossmap) {  
        ch_bowtie2_input = BOWTIE2_ASSEMBLY_BUILD.out.assembly_index
            .map { meta, assembly, index -> [ meta.group, meta, assembly, index ] }
            .combine(ch_reads_bowtie2, by: 0)
            .map { group, assembly_meta, assembly, index, reads_meta, reads -> [ assembly_meta, assembly, index, reads_meta, reads ] }
    } 
    else {
        ch_bowtie2_input = BOWTIE2_ASSEMBLY_BUILD.out.assembly_index
            .map { meta, assembly, index -> [ meta.id, meta, assembly, index ] }
            .combine(ch_reads_bowtie2, by: 0)
            .map { id, assembly_meta, assembly, index, reads_meta, reads -> [ assembly_meta, assembly, index, reads_meta, reads ] }
    }

    // Run co-mapping 
    BOWTIE2_ASSEMBLY_ALIGN ( ch_bowtie2_input )
    ch_grouped_mappings = BOWTIE2_ASSEMBLY_ALIGN.out.mappings
        .groupTuple(by: 0)
        .map { meta, assembly, bams, bais -> [ meta, assembly.sort()[0], bams, bais ] } 
    ch_versions = ch_versions.mix(BOWTIE2_ASSEMBLY_ALIGN.out.versions)


    // Perform binning with metabat2
    ch_summarizedepth_input = ch_grouped_mappings.map { meta, assembly, bams, bais ->
                                        def meta_new = meta.clone()
                                    [ meta_new, bams, bais ]
    }
    
    METABAT2_JGISUMMARIZEBAMCONTIGDEPTHS ( ch_summarizedepth_input )
    ch_versions = ch_versions.mix(METABAT2_JGISUMMARIZEBAMCONTIGDEPTHS.out.versions.first())

    ch_metabat_depths = METABAT2_JGISUMMARIZEBAMCONTIGDEPTHS.out.depth
        .map { meta, depths ->
            def meta_new = meta.clone()
            meta_new['binner'] = 'MetaBAT2'

            [ meta_new, depths ]
        }

    ch_metabat2_input = ch_grouped_mappings
        .map { meta, assembly, bams, bais ->
            def meta_new = meta.clone()
            meta_new['binner'] = 'MetaBAT2'

            [ meta_new, assembly, bams, bais ]
        }
        .join( ch_metabat_depths, by: 0 )
        .map { meta, assembly, bams, bais, depths ->
            [ meta, assembly, depths ]
        }

    METABAT2_METABAT2 ( ch_metabat2_input )
    ch_versions = ch_versions.mix(METABAT2_METABAT2.out.versions.first())

    PARSESTB ( METABAT2_METABAT2.out.fasta )
    ch_versions = ch_versions.mix(PARSESTB.out.versions)

    // Run checkM on metabat2 bins
    if (!params.skip_checkm){
        CHECKM_LINEAGEWF ( METABAT2_METABAT2.out.fasta,
                        ".fa"
        )
        ch_versions = ch_versions.mix(CHECKM_LINEAGEWF.out.versions)
    }

    // Run GTDB on metabat2 bins
    if (!params.skip_gtdb){
        GTDBTK_DB_PREPARATION ( params.gtdb )
        GTDBTK_CLASSIFYWF ( METABAT2_METABAT2.out.fasta,
                            GTDBTK_DB_PREPARATION.out
        )
        ch_versions = ch_versions.mix(GTDBTK_CLASSIFYWF.out.versions)
    }

    // *****************************
    // Perform bacteriophage binning
    // *****************************
    // https://github.dev/replikation/What_the_Phage

    // Run virsorter2
    if (params.run_virsorter){
        VIRSORTER2_DB_PREPARATION ()
        VIRSORTER2 (ch_assemblies_raw, 
                    VIRSORTER2_DB_PREPARATION.out.database)
        ch_versions = ch_versions.mix(VIRSORTER2.out.versions)
    } 
    
    // Run checkv
    // https://github.dev/replikation/What_the_Phage
    if (!params.skip_checkv){
        CHECKV_DB_PREP ( params.checkv_db )
        CHECKV (ch_assemblies_raw, 
                CHECKV_DB_PREP.out)
        ch_versions = ch_versions.mix(CHECKV.out.versions)
    }

    // Run virfinder
    if (params.run_virfinder){
        VIRFINDER ( ch_assemblies_raw )
        ch_versions = ch_versions.mix(VIRFINDER.out.versions)
    }

    // *****************************
    // Perform eukaryote binning
    // *****************************
    ch_eukrep_meta = ch_assemblies_raw
            .map { meta, ass -> [meta.id, meta]}
    ch_eukrep_ass_in = ch_assemblies_raw
            .map { meta, ass -> [meta.id, ass]}
    ch_eukrep_stbs =  PARSESTB.out.stb
            .map { meta, stb -> [meta.id, stb]}

    ch_euk_input = ch_eukrep_meta
            .combine( ch_eukrep_ass_in, by: 0)
            .combine( ch_eukrep_stbs, by: 0)
            .map {id, meta, ass, stb -> [meta, ass, stb]}

    EUKREP ( ch_euk_input )
    ch_versions = ch_versions.mix(EUKREP.out.versions)

    // ******************************
    // Perform plasmid identification
    // ******************************
    ABRICATE_RUN ( ch_assemblies_raw )
    ch_versions = ch_versions.mix(ABRICATE_RUN.out.versions)

    // THIS TOOL DOESN'T WORK AND IS DEPREACTED BY THE AUTHOR
    // PLASFLOW ( ch_assemblies_raw )
    // ch_versions = ch_versions.mix(PLASFLOW.out.versions)

    // ******************************
    // Preform gene annotation
    // ******************************
        
    DIAMOND_BLASTP (
        PRODIGAL.out.amino_acid_fasta,
        params.trep_diamond,
        'txt',
        []
    )
    ch_versions = ch_versions.mix(DIAMOND_BLASTP.out.versions)

    // Run tRep
    ch_annoinput_stbs =  PARSESTB.out.stb
            .map { meta, stb -> [meta.id, stb]}
    ch_annoinput_genes = PRODIGAL.out.amino_acid_fasta
            .map { meta, genes -> [meta.id, genes]}
    ch_annoinput_diamond = DIAMOND_BLASTP.out.txt
            .map { meta, diam -> [meta.id, diam]}
    ch_annoinput_meta = DIAMOND_BLASTP.out.txt
            .map { meta, diam -> [meta.id, meta]}

    ch_trep_input = ch_annoinput_meta
            .combine( ch_annoinput_diamond, by: 0)
            .combine( ch_annoinput_stbs, by: 0)
            .combine( ch_annoinput_genes, by: 0)
            .map {id, meta, dia, stb, genes -> [meta, dia, stb, genes]}

    TREP (
        ch_trep_input,
        params.trep_ttable
    )
    ch_versions = ch_versions.mix(TREP.out.versions)


    // ******************************
    // Finish
    // ******************************
    
    ch_sum_meta     = ch_assemblies_raw.map { meta, item -> [meta.id, meta]}

    // Run summarize
    ch_sum_fasta    = ch_assemblies_raw
            .map { meta, item -> [meta.id, item]}
    ch_sum_stbs     =  PARSESTB.out.stb
            .map { meta, item -> [meta.id, item]}
    ch_sum_euk      =  EUKREP.out.euk_bins
            .map { meta, item -> [meta.id, item]}
    ch_sum_gtdb     = GTDBTK_CLASSIFYWF.out.summary
            .map { meta, item -> [meta.id, item]}
    ch_sum_chkm     = CHECKM_LINEAGEWF.out.checkm_tsv
            .map { meta, item -> [meta.id, item]}
    ch_sum_chv      = CHECKV.out.checkV_results_ch
            .map { meta, item -> [meta.id, item]}
    ch_sum_abr      = ABRICATE_RUN.out.report
            .map { meta, item -> [meta.id, item]}
    ch_sum_tScaff   = TREP.out.trep_scaff_tax
            .map { meta, item -> [meta.id, item]}
    ch_sum_tGenome  = TREP.out.trep_genome_tax
            .map { meta, item -> [meta.id, item]}

    ch_sum_input = ch_sum_meta
            .combine( ch_sum_fasta, by: 0)
            .combine( ch_sum_stbs, by: 0)
            .combine( ch_sum_euk, by: 0)
            .combine( ch_sum_gtdb, by: 0)
            .combine( ch_sum_chkm, by: 0)
            .combine( ch_sum_chv, by: 0)
            .combine( ch_sum_abr, by: 0)
            .combine( ch_sum_tScaff, by: 0)
            .combine( ch_sum_tGenome, by: 0)
            .map {id, meta, fasta, stb, euk, gtdb, chkm, chv, abr, tScaff, tGenome -> 
                     [meta, fasta, stb, euk, gtdb, chkm, chv, abr, tScaff, tGenome]}

    // ch_sum_input = ch_sum_meta
    //         .combine( ch_sum_fasta, by: 0)
    //         .combine( ch_sum_stbs, by: 0)
    //         // .combine( ch_sum_euk, by: 0)
    //         // .combine( ch_sum_gtdb, by: 0)
    //         // .combine( ch_sum_chkm, by: 0)
    //         // .combine( ch_sum_chv, by: 0)
    //         // .combine( ch_sum_abr, by: 0)
    //         // .combine( ch_sum_tScaff, by: 0)
    //         // .combine( ch_sum_tGenome, by: 0)
    //         .map {id, meta, fasta, stb -> [meta, fasta, stb]}

    SUMMARY_BINNING (
        ch_sum_input
    )
    
    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )


}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
