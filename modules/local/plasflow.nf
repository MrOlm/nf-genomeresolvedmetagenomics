// https://github.com/replikation/docker_pipelines/blob/e1c3ca8a11b2e3ec305e15cf0b767ecbf110b87f/modules/plasflow.nf
process PLASFLOW {
    tag "$meta.id"
    label 'process_high'
    label 'error_ignore'

    container "quay.io/biocontainers/plasflow:1.1.0--py35_0"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), file("${meta.id}_chromosomes.fasta"), file("${meta.id}_plasmids.fasta"), file("${meta.id}_unclassified.fasta")
    path "versions.yml"                                 , emit: versions

    script:
    def name =  "${meta.id}"
    """
    PlasFlow.py --input ${fasta} --output ${name} --threshold 0.7

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        plasflow: 1.1.0
    END_VERSIONS
    """
}
