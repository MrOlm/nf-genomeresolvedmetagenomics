process INSTRAIN {
    resourceLabels project: 'testingTags', owner: 'mattolm', pipeline: 'profile', stage: 'inStrain_profile'
    tag "$meta.id"
    label 'process_medium' 
    label 'error_ignore'
    label 'process_long'
    label 'process_high_memory'

    // MO - UPDATE TO LASTEST inSTRAIN (v1.6.1) WHEN YOU CAN
    conda (params.enable_conda ? "bioconda::instrain=1.6.1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/instrain':
        'quay.io/biocontainers/instrain:1.6.1--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(bam)
    path genome_fasta
    path genes_fasta
    path stb_file

    output:
    tuple val(meta), path("*.IS"), emit: profile
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def genes_args = genes_fasta != [] ? "-g ${genes_fasta}": ''
    def stb_args = stb_file != [] ? "-s ${stb_file}": ''

    // TODO nf-core: Where possible, a command MUST be provided to obtain the version number of the software e.g. 1.10
    //               If the software is unable to output a version number on the command-line then it can be manually specified
    //               e.g. https://github.com/nf-core/modules/blob/master/modules/homer/annotatepeaks/main.nf
    //               Each software used MUST provide the software name and version number in the YAML version file (versions.yml)
    // TODO nf-core: It MUST be possible to pass additional parameters to the tool as a command-line string via the "task.ext.args" directive
    // TODO nf-core: If the tool supports multi-threading then you MUST provide the appropriate parameter
    //               using the Nextflow "task" variable e.g. "--threads $task.cpus"
    // TODO nf-core: Please replace the example samtools command below with your module's command
    // TODO nf-core: Please indent the command appropriately (4 spaces!!) to help with readability ;)
    """
    inStrain \\
        profile \\
        $bam \\
        $genome_fasta \\
        -o ${meta.id}.IS \\
        -p $task.cpus \\
        $genes_args \\
        $stb_args \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        instrain: \$(echo \$(inStrain profile --version 2>&1) | sed 's/^.*inStrain //; s/Using.*\$//' ))
    END_VERSIONS
    """
}

// process INSTRAIN_QUICKPROFILE {
//     tag "$meta.id"
//     label 'process_medium' 
//     label 'error_ignore'

//     // MO - UPDATE TO LASTEST inSTRAIN (v1.6.1) WHEN YOU CAN
//     conda (params.enable_conda ? "bioconda::instrain=1.6.1" : null)
//     container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
//         'https://depot.galaxyproject.org/singularity/instrain':
//         'quay.io/biocontainers/instrain:1.6.1--pyhdfd78af_0' }"

//     input:
//     tuple val(meta), path(bam)
//     path genome_fasta
//     path stb_file

//     output:
//     tuple val(meta), path("*.QP"), emit: profile
//     path "versions.yml"           , emit: versions

//     when:
//     task.ext.when == null || task.ext.when

//     script:
//     def args = task.ext.args ?: ''
//     def prefix = task.ext.prefix ?: "${meta.id}"
//     def stb_args = stb_file != [] ? "-s ${stb_file}": ''

//     // TODO nf-core: Where possible, a command MUST be provided to obtain the version number of the software e.g. 1.10
//     //               If the software is unable to output a version number on the command-line then it can be manually specified
//     //               e.g. https://github.com/nf-core/modules/blob/master/modules/homer/annotatepeaks/main.nf
//     //               Each software used MUST provide the software name and version number in the YAML version file (versions.yml)
//     // TODO nf-core: It MUST be possible to pass additional parameters to the tool as a command-line string via the "task.ext.args" directive
//     // TODO nf-core: If the tool supports multi-threading then you MUST provide the appropriate parameter
//     //               using the Nextflow "task" variable e.g. "--threads $task.cpus"
//     // TODO nf-core: Please replace the example samtools command below with your module's command
//     // TODO nf-core: Please indent the command appropriately (4 spaces!!) to help with readability ;)
//     """
//     inStrain \\
//         quick_profile \\
//         $bam \\
//         $genome_fasta \\
//         -o ${meta.id}.QP \\
//         -p $task.cpus \\
//         $stb_args \\
//         $args

//     cat <<-END_VERSIONS > versions.yml
//     "${task.process}":
//         instrain: \$(echo \$(inStrain quick_profile --version 2>&1) | sed 's/^.*inStrain //; s/Using.*\$//' ))
//     END_VERSIONS
//     """
// }
