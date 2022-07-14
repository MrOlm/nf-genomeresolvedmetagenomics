process PARSESTB {
    tag "$meta.id"
    label 'process_low'

    container "quay.io/biocontainers/drep:3.3.0--pyhdfd78af_0"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*.stb")                      , emit: stb
    path "versions.yml"                                 , emit: versions

    script:
    def name =  "${meta.id}"
    """
    parse_stb.py -f ${fasta} -o ${meta.id}.stb --reverse

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        parse_stb: \$(dRep --version | awk '{print \$2}')
    END_VERSIONS
    """
}
