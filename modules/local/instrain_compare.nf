process INSTRAINCOMPARE {
    tag "$group"
    label 'process_medium'
    label 'process_vlong'
    label 'process_high_memory'
    disk '1 TB'

    // MO - UPDATE TO LASTEST inSTRAIN (v1.6.1) WHEN YOU CAN
    conda (params.enable_conda ? "bioconda::instrain=1.6.1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/instrain':
        'quay.io/biocontainers/instrain:1.6.1--pyhdfd78af_0' }"

    input:
    tuple val(group), path(IS_list), path(bam_list)
    path(stb_file)
    // val(group)
    // path IS_list
    // path bam_list

    output:
    tuple val(group), path("*.RC"), emit: profile
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${group}"
    def stb_args = stb_file ? "-s ${stb_file}": ''
    """
    echo $bam_list

    pip install instrain --upgrade
    
    inStrain \\
        compare \\
        -i $IS_list \\
        -o ${prefix}.RC \\
        -p 4 \\
        -b $bam_list \\
        $stb_args \\
        $args \\

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        instrain: \$(echo \$(inStrain profile --version 2>&1) | sed 's/^.*inStrain //; s/Using.*\$//' ))
    END_VERSIONS
    """
}
