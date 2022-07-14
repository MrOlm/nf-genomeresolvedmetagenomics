process ASSEMBLY_FILTER_RENAME {
    label "process_low"

    conda (params.enable_conda ? "bioconda::bbmap=38.96" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bbmap:38.96' :
        'quay.io/biocontainers/bbmap:38.96--h5c4e2a8_1' }"

    input:
    tuple val(meta), path(assembly)
    val min_length

    output:
    tuple val(meta), path("${meta.id}_${meta.assembler}.fasta.gz")       , emit: assemblies 
    path "versions.yml"                                                  , emit: versions

    script: 
    def min_length = min_length ?: '1000'
    def name = "${meta.id}_${meta.assembler}"
    """
    reformat.sh in=${assembly} out=${name}.min.fasta.gz minlength=${min_length}

    zcat ${name}.min.fasta.gz |  awk 'BEGIN{c=1} {if ( substr(\$1,0,1) == ">" ) {print ">${meta.id}_scaffold_" c " " substr(\$1,2,10000000); c = c + 1} else print \$0}' | gzip > ${name}.fasta.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        reformat.sh: \$(echo \$(reformat.sh --version 2>&1) | awk '{print \$10}')
    END_VERSIONS
    """
}

