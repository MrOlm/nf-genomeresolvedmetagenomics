process GTDBTK_CLASSIFYWF {
    tag "${meta.id}"
    label 'process_medium'
    label 'process_high_memory'

    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    conda (params.enable_conda ? "bioconda::gtdbtk=2.1.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gtdbtk:2.1.0--pyhdfd78af_5' :
        'quay.io/biocontainers/gtdbtk:2.1.0--pyhdfd78af_5' }"

    input:
    tuple val(meta), path("bins/*")
    tuple val(db_name), path("database/*")

    output:
    tuple val(meta), path ("gtdbtk.${meta.id}.*.summary.tsv")        , optional:true, emit: summary
    tuple val(meta), path ("gtdbtk.${meta.id}.*.classify.tree.gz")  , optional:true, emit: tree
    tuple val(meta), path ("gtdbtk.${meta.id}.*.markers_summary.tsv"), optional:true, emit: markers
    tuple val(meta), path ("gtdbtk.${meta.id}.*.msa.fasta.gz")       , optional:true, emit: msa
    tuple val(meta), path ("gtdbtk.${meta.id}.*.user_msa.fasta")     , optional:true, emit: user_msa
    tuple val(meta), path ("gtdbtk.${meta.id}.*.filtered.tsv")       , optional:true, emit: filtered
    tuple val(meta), path ("gtdbtk.${meta.id}.log")                  ,                emit: log
    tuple val(meta), path ("gtdbtk.${meta.id}.warnings.log")         , optional:true, emit: warnings
    tuple val(meta), path ("gtdbtk.${meta.id}.failed_genomes.tsv")   , optional:true, emit: failed
    path "versions.yml"                                            ,                emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def pplacer_scratch = "--scratch_dir pplacer_tmp"
    def VERSION = '2.1.0' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    export GTDBTK_DATA_PATH="\${PWD}/database"
    if [ "${pplacer_scratch}" != "" ] ; then
        mkdir pplacer_tmp
    fi

    mkdir classify
    touch classify/gtdbtk.N5_271_010G1.bac120.summary.tsv

    pattern="bins/*.gz"
    if compgen -G "\$pattern" >/dev/null; then
        gzip -d bins/*.gz
    fi

    gtdbtk classify_wf \\
        $args \\
        --genome_dir bins \\
        --prefix "gtdbtk.${meta.id}" \\
        --out_dir "\${PWD}" \\
        --cpus $task.cpus \\
        --pplacer_cpus $task.cpus \\
        $pplacer_scratch \\
        --extension fa \\

    if [ -f "gtdbtk.${meta.id}".*.classify.tree ]; then
        gzip "gtdbtk.${meta.id}".*.classify.tree "gtdbtk.${meta.id}".*.msa.fasta
    fi

    mv gtdbtk.log "gtdbtk.${meta.id}.log"
    mv gtdbtk.warnings.log "gtdbtk.${meta.id}.warnings.log"


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gtdbtk: \$(echo \$(gtdbtk --version -v 2>&1) | sed "s/gtdbtk: version //; s/ Copyright.*//")
    END_VERSIONS
    """

    stub:
    """
    touch gtdbtk.${meta.id}.stub.summary.tsv
    touch gtdbtk.${meta.id}.stub.classify.tree.gz
    touch gtdbtk.${meta.id}.stub.markers_summary.tsv
    touch gtdbtk.${meta.id}.stub.msa.fasta.gz
    touch gtdbtk.${meta.id}.stub.user_msa.fasta
    touch gtdbtk.${meta.id}.stub.filtered.tsv
    touch gtdbtk.${meta.id}.log
    touch gtdbtk.${meta.id}.warnings.log
    touch gtdbtk.${meta.id}.failed_genomes.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gtdbtk: $VERSION
    END_VERSIONS
    """
}
