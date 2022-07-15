process DATABASE_SCAFFDREP {
    tag "$group"
    label 'process_low'
    beforeScript 'ulimit -Ss unlimited'

    container "quay.io/biocontainers/drep:3.3.0--pyhdfd78af_0"

    input:
    tuple val(group), path(gs1), path(gs2), path(gs3), path(gs4)
    val(extention)

    output:
    path "${group}.database.stb"                                ,  emit: stb
    path "${group}.database.fa.gz"                              ,  emit: fasta
    path "${group}.SL_out/DereplicationInfo.csv"                ,  emit: drep_info
    path "versions.yml"                                         ,  emit: versions

    script:
    def args = task.ext.args ?: ''
    """    
    # Make an .stb file
    find -name "*${extention}" -print > genomes.txt
    parse_stb_2.py -f genomes.txt --reverse -o ${group}.database.stb

    # Concatonate fasta sets
    cat $gs1 > set1.fa.gz
    cat $gs2 > set2.fa.gz
    cat $gs3 > set3.fa.gz
    cat $gs4 > set4.fa.gz

    # Have to uncompress for nucmer to work
    if [ -s set1.fa.gz ]; then
        gzip -d set1.fa.gz
    else
        touch set1.fa
    fi
     
    if [ -s set2.fa.gz ]; then
        gzip -d set2.fa.gz
    else
        touch set2.fa
    fi
    
    if [ -s set3.fa.gz ]; then
        gzip -d set3.fa.gz
    else
        touch set3.fa
    fi

    if [ -s set4.fa.gz ]; then
        gzip -d set4.fa.gz
    else
        touch set4.fa
    fi

    # Run program
    ScaffoldLevel_dRep2.py -f set1.fa set2.fa set3.fa set4.fa -o ${group}.SL_out -p $task.cpus ${args}

    # Make final file
    gzip ${group}.SL_out -c > ${group}.database.fa.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dRep: \$(dRep -h | awk '{if (NR==2) print substr(\$3, 2, 1000)}')
    END_VERSIONS
    """
}

process DATABASE_CONCAT {
    tag "$group"
    label 'process_low'
    beforeScript 'ulimit -Ss unlimited'

    container "quay.io/biocontainers/drep:3.3.0--pyhdfd78af_0"

    input:
    tuple val(group), path(gs1), path(gs2), path(gs3), path(gs4)
    val(extention)

    output:
    path "${group}.database.stb"                                ,  emit: stb
    path "${group}.database.${extention}"                              ,  emit: fasta
    path "versions.yml"                                         ,  emit: versions

    script:
    def args = task.ext.args ?: ''
    def extention = extention ?: '.fa.gz'
    """
    # Make a text file of all genomes
    find -name "*${extention}" -print > genomes.txt

    # Make an .stb file
    parse_stb_2.py -f genomes.txt --reverse -o ${group}.database.stb

    # Make an .fasta file
    { xargs cat < genomes.txt ; } > ${group}.database.${extention}
    # gzip ${group}.database.

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        Database concat: 1.0.0
    END_VERSIONS
    """
}