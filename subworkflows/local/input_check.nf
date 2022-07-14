//
// Check input samplesheet and get read channels
//

include { SAMPLESHEET_CHECK     } from '../../modules/local/samplesheet_check'

workflow INPUT_CHECK {
    take:
    samplesheet // file: /path/to/samplesheet.csv

    main:
    SAMPLESHEET_CHECK ( samplesheet )
        .csv
        .splitCsv ( header:true, sep:',' )
        .map { create_fastq_channel(it) }
        .set { reads }

    emit:
    reads                                     // channel: [ val(meta), [ reads ] ]
    versions = SAMPLESHEET_CHECK.out.versions // channel: [ versions.yml ]
}

// Function to get list of [ meta, [ fastq_1, fastq_2 ] ]
def create_fastq_channel(LinkedHashMap row) {
    // create meta map
    def meta = [:]
    meta.id         = row.sample
    meta.single_end = row.single_end.toBoolean()

    // add path(s) of the fastq file(s) to the meta map
    def fastq_meta = []
    if (!file(row.fastq_1).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> Read 1 FastQ file does not exist!\n${row.fastq_1}"
    }
    if (meta.single_end) {
        fastq_meta = [ meta, [ file(row.fastq_1) ] ]
    } else {
        if (!file(row.fastq_2).exists()) {
            exit 1, "ERROR: Please check input samplesheet -> Read 2 FastQ file does not exist!\n${row.fastq_2}"
        }
        fastq_meta = [ meta, [ file(row.fastq_1), file(row.fastq_2) ] ]
    }
    return fastq_meta
}

def hasExtension(it, extension) {
    it.toString().toLowerCase().endsWith(extension.toLowerCase())
}

workflow INPUT_CHECK_ASS {
    main:
    if(hasExtension(params.input, "csv")){
        // extracts read files from samplesheet CSV and distribute into channels
        ch_input_rows = Channel
            .from(file(params.input))
            .splitCsv(header: true)
            .map { row ->
                    if (row.size() == 5) {
                        def id = row.sample
                        def group = row.group
                        def sr1 = row.fastq_1 ? file(row.fastq_1, checkIfExists: true) : false
                        def sr2 = row.fastq_2 ? file(row.fastq_2, checkIfExists: true) : false
                        def ass = row.assembly ? file(row.assembly, checkIfExists: true) : false

                        // Check if given combination is valid
                        if (!ass) exit 1, "Invalid input samplesheet: assembly can not be empty."
                        return [ id, group, sr1, sr2, ass ]
                    } else {
                        exit 1, "Input samplesheet contains row with ${row.size()} column(s). Expects 5."
                    }
                }
        // make the metamap
        ch_ass = ch_input_rows
            .map { id, group, sr1, sr2, ass ->
                        def meta = [:]
                        meta.id           = id
                        meta.group        = group
                        meta.crossmap     = !params.skip_crossmap
                        if (!params.skip_crossmap)
                            return [ meta, ass, [sr1, sr2] ]
                        else
                            return [ meta, ass ]
                }

    } else {
        exit 1, "Input samplesheet must be a .csv file"
    }

    // Ensure sample IDs are unique
    ch_input_rows
        .map { id, group, sr1, sr2, lr -> id }
        .toList()
        .map { ids -> if( ids.size() != ids.unique().size() ) {exit 1, "ERROR: input samplesheet contains duplicated sample IDs!" } }

    emit:
    ch_assemblies = ch_ass
}

workflow INPUT_CHECK_DREP {
    main:
    if(hasExtension(params.input, "csv")){
        // extracts read files from samplesheet CSV and distribute into channels
        ch_input_rows = Channel
            .from(file(params.input))
            .splitCsv(header: true)
            .map { row ->
                    if (row.size() == 2) {
                        def id = row.sample
                        def drep_sf = row.drep_sf
                        def group = '0'

                        // Check if given combination is valid
                        if (!drep_sf) exit 1, "Invalid input samplesheet: assembly can not be empty."
                        return [ id, drep_sf, group ]

                    } else if (row.size() == 3) {
                        def id = row.sample
                        def drep_sf = row.drep_sf
                        def group = row.group

                        // Check if given combination is valid
                        if (!drep_sf) exit 1, "Invalid input samplesheet: assembly can not be empty."
                        return [ id, drep_sf, group ]
                    } else {
                        exit 1, "Input samplesheet contains row with ${row.size()} column(s). Expects 2 or 3."
                    }
                }
        // make the metamap
        ch_ass = ch_input_rows
            .map { id, drep_sf, group ->
                        def meta = [:]
                        meta.id           = id
                        meta.group        = group
                        return [ meta, drep_sf ]
            }

    } else {
        exit 1, "Input samplesheet must be a .csv file"
    }

    // Ensure sample IDs are unique
    ch_input_rows
        .map { id, drep_sf, group -> id }
        .toList()
        .map { ids -> if( ids.size() != ids.unique().size() ) {exit 1, "ERROR: input samplesheet contains duplicated sample IDs!" } }

    emit:
    ch_drep_sf = ch_ass
}

