process DREP {
    tag "$group"
    label 'process_medium'

    container "quay.io/biocontainers/drep:3.3.0--pyhdfd78af_0"

    input:
    tuple val(group), path(genomes), path(genome_info), path(extra_score)
    val(parameters)

    output:
    // tuple val(group), path ("${group}/dereplicated_genomes/*"),  emit: genomes
    // tuple val(group), path ("${group}/data_tables/*"),           emit: tables, optional:true
    path "versions.yml"                                 ,           emit: versions

    script:
    def args = task.ext.args ?: ''
    """
    # Figure out if there's just 1 genome
    genomes=( ${genomes} )
    if [ \${#genomes[@]} -eq 1 ]
    then
        echo "only 1 genome"
        mkdir -p ${group}/dereplicated_genomes/
        cp $genomes ${group}/dereplicated_genomes/
    else
        dRep dereplicate ${group} -g ${genomes} --genomeInfo ${genome_info} -extraW ${extra_score} -p $task.cpus ${parameters} ${args}
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dRep: \$(dRep -h | awk '{if (NR==2) print substr(\$3, 2, 1000)}')
    END_VERSIONS
    """
}

process PREPARE_DREP {
    tag "$group"
    label 'process_low'

    container "quay.io/biocontainers/drep:3.3.0--pyhdfd78af_0"

    input:
    tuple val(group), path(fasta)

    output:
    tuple val(group), path ("${group}.set1_bins/*"), path ("${group}.set1_genomeInfo.csv"), path ("${group}.set1_extraScore.tsv"), emit: set1
    tuple val(group), path ("${group}.set2_bins/*"), path ("${group}.set2_genomeInfo.csv"), path ("${group}.set2_extraScore.tsv"), emit: set2, optional: true
    tuple val(group), path ("${group}.set3_bins/*"), path ("${group}.set3_genomeInfo.csv"), path ("${group}.set3_extraScore.tsv"), emit: set3, optional: true
    tuple val(group), path ("${group}.set4_bins/*"), path ("${group}.set4_genomeInfo.csv"), path ("${group}.set4_extraScore.tsv"), emit: set4, optional: true
    path "versions.yml"                                 , emit: versions

    script:
    """
    prepare_drep.py -f $fasta -o $group

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dRep: \$(dRep -h | awk '{if (NR==2) print substr(\$3, 2, 1000)}')
    END_VERSIONS
    """
}

// process PARSESTB {
//     tag "$meta.id"
//     label 'process_low'

//     container "quay.io/biocontainers/drep:3.3.0--pyhdfd78af_0"

//     input:
//     tuple val(meta), path(fasta)

//     output:
//     tuple val(meta), path("*.stb")                      , emit: stb
//     path "versions.yml"                                 , emit: versions

//     script:
//     def name =  "${meta.id}"
//     """
//     parse_stb.py -f ${fasta} -o ${meta.id}.stb --reverse

//     cat <<-END_VERSIONS > versions.yml
//     "${task.process}":
//         parse_stb: \$(dRep --version | awk '{print \$2}')
//     END_VERSIONS
//     """
// }
