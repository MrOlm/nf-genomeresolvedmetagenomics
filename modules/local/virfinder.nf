process VIRFINDER {
    tag "$meta.id"
    label 'process_low'
    label 'process_long'

    container "multifractal/virfinder:0.2"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("${meta.id}_*.list.gz")              , emit: results
    path "versions.yml"                                     , emit: versions

    script:
    def name =  "${meta.id}"
    """
    fasta=$fasta
    if [[ "\${fasta: -3}" == ".gz" ]]
    then
        gzip -f -d $fasta
        fasta=\${fasta::-3}
    fi
    
    virfinder_execute.R \$fasta 
    cp results.txt ${name}_virfinder.list
    gzip *.list

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        VirFinder: \$(VirFinder --version | awk '{print \$2}')
    END_VERSIONS
    """
}
