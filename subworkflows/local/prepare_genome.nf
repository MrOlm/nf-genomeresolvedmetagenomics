//
// Uncompress and prepare reference genome files
//

include { GUNZIP        } from '../../modules/nf-core/modules/gunzip/main'
include { UNTAR         } from '../../modules/nf-core/modules/untar/main'
include { BOWTIE2_BUILD } from '../../modules/nf-core/modules/bowtie2/build/main'

workflow PREPARE_GENOME {
    main:

    ch_versions = Channel.empty()

    //
    // Uncompress genome fasta file if required
    //
    if (params.fasta.endsWith('.gz')) {
        GUNZIP (
            [ [:], params.fasta ]
        )
        ch_fasta    = GUNZIP.out.gunzip.map { it[1] }
        ch_versions = ch_versions.mix(GUNZIP.out.versions)
    } 
    else {
        ch_fasta = file(params.fasta)
    }

    //
    // Prepare bowtie2 index
    //
    ch_bowtie2_index = Channel.empty()
    if (params.bowtie2_index) {
        if (params.bowtie2_index.endsWith('.tar.gz')) {
            UNTAR_BOWTIE2_INDEX (
                params.bowtie2_index
            )
            ch_bowtie2_index = UNTAR_BOWTIE2_INDEX.out.untar
            ch_versions      = ch_versions.mix(UNTAR_BOWTIE2_INDEX.out.versions)
        } else {
            ch_bowtie2_index = file(params.bowtie2_index)
        }
    } else {
        BOWTIE2_BUILD (
            ch_fasta
        )
        ch_bowtie2_index = BOWTIE2_BUILD.out.index
        ch_versions      = ch_versions.mix(BOWTIE2_BUILD.out.versions)
    }

    emit:
    fasta                = ch_fasta                // path: genome.fasta
    bowtie2_index        = ch_bowtie2_index        // path: bowtie2/index/
    versions             = ch_versions             // channel: [ versions.yml ]
}
