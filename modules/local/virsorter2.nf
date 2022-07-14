process VIRSORTER2_DB_PREPARATION {
    container "quay.io/biocontainers/gnu-wget:1.18--h7132678_6"
    label 'process_low'

    output:
    path "db/*"            , emit: database

    script:
    """
    wget --no-check-certificate https://osf.io/v46sc/download -O db.tgz
    tar -zxvf db.tgz
    chmod -R a+rX db
    rm db.tgz
    """
}

process VIRSORTER2 {
    tag "$meta.id"
    label 'process_high'
    label 'error_ignore'
    label 'process_long'

    container "staphb/virsorter2:2.1"

    input:
        tuple val(meta), path(fasta)
        path("db/*")

    output:
        tuple val(meta), path("virsorter2_*.out/final-viral-score.tsv"), emit: virsorter_scores
        tuple val(meta), path("virsorter2_*.out")                      , emit: virsorter_out  
        path "versions.yml"                                            , emit: versions

    script:
    def name =  "${meta.id}"
    """
    virsorter run \
        --db-dir db \
        -w virsorter2_${name}.out \
        -i ${fasta} \
        -j ${task.cpus} \
        --provirus-off \

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        virsorter: \$(virsorter --version | tail -n 1 | awk '{print \$3}')
    END_VERSIONS
    """
}
