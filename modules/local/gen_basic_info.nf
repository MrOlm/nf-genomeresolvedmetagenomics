process BASICINFO_MERGE {
    label "process_low"

    conda (params.enable_conda ? "conda-forge::python=3.8.3" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.8.3' :
        'quay.io/biocontainers/python:3.8.3' }"

    input:
    path basic_info_list

    output:
    path '*.csv'       , emit: csv

    script: 
    """
    # Make into an array
    IFS=' ' read -r -a read_array <<< "${basic_info_list}"

    # Write the header
    cat "\${read_array[0]}" | head -n 1 > basic_info.csv

    # Loop array
    for i in "\${read_array[@]}"
    do
        cat "\$i" | tail -n +2 >> basic_info.csv
    done
    """
}

process BASICINFO {
    tag "$meta.id"
    label "process_low"

    conda (params.enable_conda ? "conda-forge::python=3.8.3" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.8.3' :
        'quay.io/biocontainers/python:3.8.3' }"

    input:
    tuple val(meta), path(final_reads)
    tuple val(metaO), path(ori_reads)

    output:
    path '*.csv'       , emit: csv

    script: 
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # Write the header
    echo "sample,read1,read2,reads,bases,RL,ori_reads,ori_bases" > ${meta.id}_basic_info.csv

    # Calculate reads and bases for the final samples
    reads_bases=\$(zcat ${final_reads[0]} | awk '{if (NR%4 == 2) tBases += length(\$0); if (NR%4 == 2) tReads += 1} END {printf tReads * 2 "," tBases ","; printf "%.1f", tBases / (tReads * 2) }')

    # Calculate reads and bases for the ori samples
    reads_bases_ori=\$(zcat ${ori_reads[0]} | awk '{if (NR%4 == 2) tBases += length(\$0); if (NR%4 == 2) tReads += 1} END {printf tReads * 2 "," tBases}')
    
    echo "${meta.id},${final_reads[0]},${final_reads[1]},\$reads_bases,\$reads_bases_ori" >> ${meta.id}_basic_info.csv
    """
}

// process BASICINFO_TOGETHER {
//     label "process_low"

//     conda (params.enable_conda ? "conda-forge::python=3.8.3" : null)
//     container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
//         'https://depot.galaxyproject.org/singularity/python:3.8.3' :
//         'quay.io/biocontainers/python:3.8.3' }"

//     input:
//     val meta_list
//     path read_list
//     // tuple val(meta), path(final_reads)
//     // tuple val(metaO), path(ori_reads)

//     output:
//     path '*.csv'       , emit: csv

//     script: 
//     """
//     echo ${meta_list} > testtttt.csv

//     echo ${read_list} >> testtttt.csv

//     IFS=' ' read -r -a read_array <<< "${read_list}"

//     echo \$read_array >> testtttt.csv

//     echo \$read_array >> testtttt.csv

//     echo "break" >> testtttt.csv

//     echo ${read_list[0]} >> testtttt.csv

//     echo ${read_list[1]} >> testtttt.csv

//     echo "\${read_array}" >> testtttt.csv

//     echo "\${read_array[1]}" >> testtttt.csv


//     """
//     // echo "\${array[1]}" >> testtttt.csv
//     // for f in \$(${final_reads});
//     //     do echo \$f >> testtttt.csv;
//     // done;
// }




