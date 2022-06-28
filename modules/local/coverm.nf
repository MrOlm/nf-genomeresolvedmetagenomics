process COVERM {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::coverm=0.6.1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/instrain':
        'quay.io/biocontainers/coverm:0.6.1--h1535e20_4' }"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("${meta.id}.coverm.tsv"), path("${meta.id}.coverm.log") , emit: coverm
    path "versions.yml"                                                       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    coverm contig \\
        -b $bam \\
        -t $task.cpus \\
        -m mean covered_bases length count --output-format sparse --no-zeros \\
        $args \\
        -o ${meta.id}.coverm.tsv \\
        2> ${meta.id}.coverm.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        coverm: \$(echo \$(coverm --version 2>&1) | sed 's/^.*coverm //; s/Using.*\$//' ))
    END_VERSIONS
    """
}
