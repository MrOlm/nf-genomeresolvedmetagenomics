process SUMMARY_BINNING {
    tag "$meta.id"
    label 'process_low'

    container "quay.io/biocontainers/drep:3.3.0--pyhdfd78af_0"

    input:
    tuple val(meta), path(fasta), path(stb), path(euk), path(gtdb), path(chkm), path(chv), path(abr), path(tScaff), path(tGenome)
    // tuple val(meta), path(fasta), path(stb)

    output:
    tuple val(meta), path ("${meta.id}.sum.AllGenomes.csv")     , emit: genome_info
    tuple val(meta), path ("${meta.id}.sum.PhagePlasmid.csv")   , emit: pp_info
    tuple val(meta), path ("${meta.id}.sum.stb")   , emit: stb
    tuple val(meta), path ("dRep_sf_${meta.id}/")               , emit: drep_sf
    path "versions.yml"                                 , emit: versions

    script:
    def name =  "${meta.id}"
    """
    summarize_binning.py \\
        -f ${fasta} \\
        -s ${stb} \\
        -e ${euk} \\
        -g ${gtdb} \\
        -chm ${chkm} \\
        -chv ${chv} \\
        -a ${abr} \\
        -ts ${tScaff} \\
        -tg ${tGenome} \\
        -o ${meta.id} \\

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        SummarizeBinning: \$(summarize_binning.py --version | awk '{print \$2}')
    END_VERSIONS
    """
}
