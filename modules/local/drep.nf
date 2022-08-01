process DREP {
    tag "$group.$set"
    label 'process_high'
    label 'process_high_memory'
    label 'process_long'

    container "sonnenburglab/drep:3.4.0"

    input:
    tuple val(group), val(set), path(genomes, stageAs: 'genomes-tmp'), path(genome_info), path(extra_score)
    val(parameters)
    val(extention)

    output:
    tuple val(group), path ("${group}.${set}/dereplicated_genomes/")  ,  emit: genomes
    tuple val(group), path ("${group}.${set}/data_tables/")           ,  emit: tables, optional:true
    path "versions.yml"                                         ,  emit: versions

    script:
    def args = task.ext.args ?: ''
    def extention = extention ?: '.fa.gz'
    """
    echo "start"
    
    # Make a text file of all genomes
    find genomes-tmp/ -name "*" -type f > genomes.txt
    
    # Figure out if there's just 1 genome
    if [ \$(cat genomes.txt | wc -l) -eq 1 ]
    then
        echo "only 1 genome"
        mkdir -p ${group}.${set}/dereplicated_genomes/
        cp \$(find genomes-tmp/ -name "*" -type f) ${group}.${set}/dereplicated_genomes/
    else
        dRep dereplicate ${group}.${set} -g genomes.txt --genomeInfo ${genome_info} -extraW ${extra_score} -p $task.cpus ${parameters} ${args}
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dRep: \$(dRep -h | awk '{if (NR==2) print substr(\$3, 2, 1000)}')
    END_VERSIONS
    """
}

// process DREP {
//     tag "$group.$set"
//     label 'process_high'
//     beforeScript 'ulimit -Ss unlimited'

//     container "quay.io/biocontainers/drep:3.3.0--pyhdfd78af_0"

//     input:
//     tuple val(group), val(set), path(genomes), path(genome_info), path(extra_score)
//     val(parameters)
//     val(extention)

//     output:
//     tuple val(group), path ("${group}.${set}/dereplicated_genomes/*")  ,  emit: genomes
//     tuple val(group), path ("${group}.${set}/data_tables/*")           ,  emit: tables, optional:true
//     path "versions.yml"                                         ,  emit: versions

//     script:
//     def args = task.ext.args ?: ''
//     def extention = extention ?: '.fa.gz'
//     """
//     echo "start"
    
//     # Make a text file of all genomes
//     find -name "*${extention}" -print > genomes.txt
//     genomes=\$(find -name "*${extention}" -print)
    
//     # Figure out if there's just 1 genome
//     if [ \$(cat genomes.txt | wc -l) -eq 1 ]
//     then
//         echo "only 1 genome"
//         mkdir -p ${group}.${set}/dereplicated_genomes/
//         cp \$genomes ${group}.${set}/dereplicated_genomes/
//     else
//         dRep dereplicate ${group}.${set} -g genomes.txt --genomeInfo ${genome_info} -extraW ${extra_score} -p $task.cpus ${parameters} ${args}
//     fi

//     cat <<-END_VERSIONS > versions.yml
//     "${task.process}":
//         dRep: \$(dRep -h | awk '{if (NR==2) print substr(\$3, 2, 1000)}')
//     END_VERSIONS
//     """
// }

process PREPARE_DREP {
    tag "$group"
    label 'process_low'

    container "quay.io/biocontainers/drep:3.3.0--pyhdfd78af_0"

    input:
    tuple val(group), path(fasta, stageAs: 'fastas/*')
    val(es)

    output:
    tuple val(group), val('set1'), path ("${group}.set1_bins/"), path ("${group}.set1_genomeInfo.csv"), path ("${group}.set1_extraScore.tsv"), emit: set1
    tuple val(group), val('set2'), path ("${group}.set2_bins/"), path ("${group}.set2_genomeInfo.csv"), path ("${group}.set2_extraScore.tsv"), emit: set2, optional: true
    tuple val(group), val('set3'), path ("${group}.set3_bins/"), path ("${group}.set3_genomeInfo.csv"), path ("${group}.set3_extraScore.tsv"), emit: set3, optional: true
    tuple val(group), val('set4'), path ("${group}.set4_bins/"), path ("${group}.set4_genomeInfo.csv"), path ("${group}.set4_extraScore.tsv"), emit: set4, optional: true
    path "versions.yml"                                 , emit: versions

    script:
    def es = es ?: '0'
    """
    prepare_drep.py -f fastas/* -o $group -e $es

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
