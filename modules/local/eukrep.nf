process EUKREP {
    tag "$meta.id"
    label 'process_low'

    container "sonnenburglab/eukrep:0.6.7"

    input:
    tuple val(meta), path(fasta), path(stb)

    output:
    tuple val(meta), path ("${meta.id}.eukrep")                            , emit: euk_scaffolds
    tuple val(meta), path ("${meta.id}.eukrep.csv")                        , emit: euk_bins
    path "versions.yml"                                                  , emit: versions

    script:
    def name =  "${meta.id}"
    """
    fasta=$fasta
    if [[ "\${fasta: -3}" == ".gz" ]]
    then
        gzip -d $fasta
        fasta=\${fasta::-3}
    fi
    
    EukRep -i \$fasta -o ${meta.id}.eukrep

    parse_eukrep.py -f \$fasta -e ${meta.id}.eukrep -s ${stb} -o ${meta.id}.eukrep.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        EukRep: \$(EukRep --version | awk '{print \$2}')
    END_VERSIONS
    """
}
