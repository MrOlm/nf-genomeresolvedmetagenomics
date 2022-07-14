process TREP {
    tag "$meta.id"
    label 'process_low'

    container "sonnenburglab/annotation:0.5.3"

    input:
    tuple val(meta), path(blast_results), path(stb_file),  path(faa_file)
    path(translation_file)

    output:
    tuple val(meta), path ("${meta.id}.trep.tax_fullScaffoldTaxonomy.tsv.gz")          , emit: trep_scaff_tax
    tuple val(meta), path ("${meta.id}.trep.tax_fullGeneTaxonomy.tsv.gz")              , emit: trep_gene_tax
    tuple val(meta), path ("${meta.id}.trep.tax_fullGenomeTaxonomy.tsv.gz")            , emit: trep_genome_tax, optional: true
    tuple val(meta), path ("${meta.id}.trep.function_geneFunctionalAnnotation.tsv.gz") , emit: trep_function
    path "versions.yml"                                             , emit: versions

    script:
    def name =  "${meta.id}"
    def stb_args = stb_file ? "-stb ${stb_file}": ''
    def gene_args = faa_file ? "-a ${faa_file}": ''
    """
    stb="$stb_args"
    if [[ "\${stb: -3}" == ".gz" ]]
    then
        stb_file=\$(echo \$stb | awk '{print \$2}')
        gzip -d \$stb_file
        stb="-stb \${stb_file::-3}"
    fi

    aa="$gene_args"
    if [[ "\${aa: -3}" == ".gz" ]]
    then
        faa_file=\$(echo \$aa | awk '{print \$2}')
        gzip -d \$faa_file
        aa="-a \${faa_file::-3}"
    fi

    tax_collector.py  -b $blast_results -o ${meta.id}.trep.tax \$stb \$aa
    functional_tax.py -b $blast_results -o ${meta.id}.trep.function -d ${translation_file}

    gzip *.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tRep: \$(make_Tdb.py --version | awk '{print \$2}')
    END_VERSIONS
    """
}
