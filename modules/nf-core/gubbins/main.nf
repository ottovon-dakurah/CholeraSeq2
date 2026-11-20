process GUBBINS {
    label 'process_medium'
    tag "${meta.id}"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gubbins:3.4.1--py310pl5321h577a1d6_0' :
        'biocontainers/gubbins:3.4.1--py310pl5321h577a1d6_0' }"

    input:
    tuple val(meta), path(alignment)

    output:
    tuple val(meta), path(alignment), path("*.gff")    , emit: fasta_gff
    path "*.fasta"                                     , emit: fasta
    path "*.gff"                                       , emit: gff
    path "*.vcf"                                       , emit: vcf
    path "*.csv"                                       , emit: stats
    path "*.phylip"                                    , emit: phylip
    path "*.recombination_predictions.embl"            , emit: embl_predicted
    path "*.branch_base_reconstruction.embl"           , emit: embl_branch
    path "*.final_tree.tre"                            , emit: tree
    path "*.node_labelled.final_tree.tre"              , emit: tree_labelled
    path "versions.yml"                                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    mkdir numba_cache_dir
    export NUMBA_CACHE_DIR='./numba_cache_dir'

    # Gubbins' ancestral reconstruction (pyjar) hands the phylogeny to worker
    # processes, which pickles it. Pickling recurses once per node, so a large or
    # unbalanced tree exceeds Python's default 1000-frame limit and the run dies
    # after the tree has already been built, hours in:
    #   RecursionError: maximum recursion depth exceeded while pickling an object
    # Python imports any `sitecustomize` module found on PYTHONPATH at interpreter
    # start-up, before Gubbins runs, so raising the limit there fixes this without
    # patching Gubbins itself. Deeper recursion also needs a larger thread stack,
    # set before the worker pool creates its threads; the stack size is a virtual
    # reservation and costs no resident memory.
    mkdir -p pyfix
    printf '%s\\n' \\
        'import sys, threading' \\
        'sys.setrecursionlimit(200000)' \\
        'threading.stack_size(512 * 1024 * 1024)' \\
        > pyfix/sitecustomize.py
    export PYTHONPATH="\$PWD/pyfix\${PYTHONPATH:+:\$PYTHONPATH}"

    # pyjar's ancestral reconstruction runs under Python multiprocessing, whose
    # Manager process creates its temp directory under TMPDIR. Nextflow points
    # TMPDIR at the task work directory, and where that is a network filesystem a
    # file deleted while still open is silly-renamed to .nfsXXXXXX and cannot be
    # unlinked, so the Manager's shutdown cleanup raises
    #   OSError: [Errno 16] Device or resource busy: '.nfsXXXXXX'
    # The Manager then exits while its workers are still running and the task
    # hangs indefinitely instead of failing: no exit status is written, one core
    # idles, and only the wall-clock limit ever ends it. Keeping the Manager's
    # temp directory on container-local storage avoids the rename entirely.
    export TMPDIR=/tmp

    run_gubbins.py \\
        --threads $task.cpus \\
        $args \\
        $alignment

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gubbins: 3.4.1
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    """
    mkdir numba_cache_dir
    export NUMBA_CACHE_DIR='./numba_cache_dir'

    touch ${alignment.baseName}.fasta
    touch ${alignment.baseName}.gff
    touch ${alignment.baseName}.vcf
    touch ${alignment.baseName}.csv
    touch ${alignment.baseName}.phylip
    touch ${alignment.baseName}.recombination_predictions.embl
    touch ${alignment.baseName}.branch_base_reconstruction.embl
    touch ${alignment.baseName}.final_tree.tre
    touch ${alignment.baseName}.node_labelled.final_tree.tre

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gubbins: 3.4.1
    END_VERSIONS
    """
}
