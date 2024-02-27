process STROBEALIGN {
    tag "$meta.id"
    label 'process_high'

    conda (params.enable_conda ? 'bioconda::bowtie2=2.4.4' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bowtie2:2.4.4--py39hbb4e92a_0' :
        'mattolm/strobealign-samtools' }"

    input:
    tuple val(meta), path(reads)
    path  fasta

    output:
    tuple val(meta), path("*.bam")    , emit: bam
    tuple val(meta), path("*.log")    , emit: log
    path  "versions.yml"              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    reads_args = "${reads[0]} ${reads[1]}"
    """
    strobealign \\
        $fasta \\
        $reads_args \\
        -t $task.cpus \\
        -U \\
        2> ${prefix}.strobealign.log \\
        | samtools sort --threads $task.cpus -o ${prefix}.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        strobealign: \$(echo \$(strobealign --version 2>&1) | sed 's/^.*bowtie2-align-s version //; s/ .*\$//')
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
