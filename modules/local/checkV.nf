process CHECKV {
    tag "$meta.id"
    label 'process_high'

    container "quay.io/biocontainers/checkv:0.9.0--pyhdfd78af_0"

    input:
    tuple val(meta), path(fasta)
    path "database/*"

    output:
    tuple val(meta), path("${meta.id}_quality_summary.tsv"), emit: sample_quality_ch optional true
    tuple val(meta), path("${meta.id}_results/")           , emit: checkV_results_ch optional true
    path "versions.yml"                                 , emit: versions

    script:
    def name =  "${meta.id}"
    """
    fasta=$fasta
    if [[ "\${fasta: -3}" == ".gz" ]]
    then
        gzip -d $fasta
        fasta=\${fasta::-3}
    fi

    checkv end_to_end \$fasta -d database/ -t ${task.cpus} ${name}_results  

    cp ${name}_results/quality_summary.tsv ${name}_quality_summary.tsv 

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        checkv: \$(checkv -h | awk '{print substr(\$2,2,length(\$2)-2)}' | head -n 1)
    END_VERSIONS
    """
}
