process CHECKVDB_PREPARATION {
    tag "${database}"

    container "quay.io/biocontainers/diamond:2.0.15--hb97b32f_0"

    input:
    path(database)

    output:
    path("database/*")

    script:
    """
    mkdir database
    tar -xzf ${database} -C database --strip 1

    if [[ ! -f "database/genome_db/checkv_reps.dmnd" ]]; then
        diamond makedb --in database/genome_db/checkv_reps.faa --db database/genome_db/checkv_reps
    fi

    rm ${database}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tar: \$(tar --version 2>&1 | sed -n 1p | sed 's/tar (GNU tar) //')
    END_VERSIONS
    """
}
