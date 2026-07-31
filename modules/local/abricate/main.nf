// Local module - based on nf-core/modules ABRICATE_RUN (master branch),
// adapted to CholeraSeq's classic versions.yml convention. databasedir is
// optional - pass [] to use ABRicate's bundled default databases.
process ABRICATE_RUN {
    tag "$meta.id"
    label 'process_medium'

    container "quay.io/biocontainers/abricate:1.0.1--ha8f3691_1"

    input:
    tuple val(meta), path(assembly)
    path databasedir

    output:
    tuple val(meta), path("*.txt"), emit: report
    path 'versions.yml'           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args    = task.ext.args   ?: ''
    def prefix  = task.ext.prefix ?: "${meta.id}"
    def datadir = databasedir ? "--datadir ${databasedir}" : ''
    """
    if [[ "${assembly}" != "${prefix}.fasta" ]]; then
        ln -s ${assembly} ${prefix}.fasta
    fi

    abricate \\
        ${prefix}.fasta \\
        ${args} \\
        ${datadir} \\
        --threads ${task.cpus} \\
        > ${prefix}.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        abricate: \$(abricate --version | sed 's/^.* //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.txt
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        abricate: 1.0.1
    END_VERSIONS
    """
}
