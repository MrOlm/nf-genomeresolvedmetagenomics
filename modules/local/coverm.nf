process COVERM {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::coverm=0.6.1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/instrain':
        'mattolm/epiricial-detection-docker:3.0' }"

    input:
    tuple val(meta), path(bam)
    path genome_fasta
    path stb_file

    output:
    tuple val(meta), path("${meta.id}.coverm.tsv"), path("${meta.id}.coverm.log"), path("${meta.id}.coverm.parsed.tsv") , emit: coverm
    path "versions.yml"                                                       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def fasta_stb_present = (genome_fasta && !genome_fasta.trim().empty) && (stb_file && !stb_file.trim().empty) ? 1 : 0

    """
    coverm contig \\
        -b $bam \\
        -t $task.cpus \\
        -m mean covered_bases length count --output-format sparse --no-zeros \\
        $args \\
        -o ${meta.id}.coverm.tsv \\
        2> ${meta.id}.coverm.log

    if [ "$fasta_stb_present" == "true" ]; then
        echo parse_coverm.py -c ${meta.id}.coverm.tsv -f $genome_fasta -s $stb_file -o ${meta.id}.coverm.parsed.tsv
        parse_coverm.py -c ${meta.id}.coverm.tsv -f $genome_fasta -s $stb_file -o ${meta.id}.coverm.parsed.tsv
    else
        echo 'Skipping parse_coverm.py due to missing inputs'
        touch ${meta.id}.coverm.parsed.tsv
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        coverm: \$(echo \$(coverm --version 2>&1) | sed 's/^.*coverm //; s/Using.*\$//' ))
    END_VERSIONS
    """
}
