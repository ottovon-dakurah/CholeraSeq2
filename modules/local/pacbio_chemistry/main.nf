// Detect PacBio chemistry (HiFi vs CLR) from the reads themselves.
//
// Submitter labels, instrument model, read length and reported base quality all
// fail to predict which Flye mode will work: public runs are frequently
// mislabelled, many report no instrument model, and some carry placeholder
// quality strings. What does separate the two chemistries is the error rate of
// the reads, which shows up directly as the divergence between overlapping
// reads of the same sample.
//
// minimap2's all-versus-all mode reports that divergence per overlap in the
// dv:f: tag. HiFi reads (>99% accurate) overlap each other at a few tenths of a
// percent; CLR reads differ by around 8-12%. The two groups are separated by
// roughly half an order of magnitude, so a single threshold classifies reliably.
//
// Only a subsample of reads is used, so this costs a minute or two per sample
// against the hours a wrongly-assembled sample wastes.
process PACBIO_CHEMISTRY {
    tag "$meta.id"
    label 'process_medium'

    container "quay.io/biocontainers/minimap2:2.24--h7132678_1"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path(reads), env('CHEMISTRY'), emit: typed
    tuple val(meta), path("*.chemistry.txt"), emit: report
    path 'versions.yml', emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix    = task.ext.prefix ?: "${meta.id}"
    def nreads    = params.pacbio_detect_reads    ?: 20000
    def threshold = params.pacbio_detect_threshold ?: 0.03
    def fallback  = (params.pacbio_mode == 'hifi') ? 'pacbio_hifi' : 'pacbio_clr'
    """
    # Whole FASTQ records only: 4 lines per read. head closes the pipe once it has
    # enough, which kills zcat with SIGPIPE; pipefail would turn that into a failed
    # task, so the subsampling runs with pipefail off.
    set +o pipefail
    zcat -f ${reads} | head -n \$(( ${nreads} * 4 )) > subsample.fastq
    set -o pipefail

    minimap2 -x ava-pb -t ${task.cpus} subsample.fastq subsample.fastq 2> minimap2.log \\
        | grep -o 'dv:f:[0-9.]*' | cut -d: -f3 | sort -n > divergences.txt || true

    N=\$(wc -l < divergences.txt)
    if [ "\$N" -ge 100 ]; then
        MEDIAN=\$(awk '{a[NR]=\$1} END {printf "%.4f", a[int(NR/2)+1]}' divergences.txt)
        CHEMISTRY=\$(awk -v m="\$MEDIAN" -v t=${threshold} 'BEGIN {print (m < t) ? "pacbio_hifi" : "pacbio_clr"}')
        BASIS="median read-to-read divergence \$MEDIAN over \$N overlaps (threshold ${threshold})"
    else
        # Too few overlaps to judge (very low coverage, or a tiny archive record).
        # Fall back to --pacbio_mode: CLR mode tolerates accurate reads, while
        # HiFi mode fails outright on error-prone reads, so this errs safely.
        MEDIAN=NA
        CHEMISTRY=${fallback}
        BASIS="only \$N overlaps found - too few to judge; defaulted to \$CHEMISTRY"
    fi
    export CHEMISTRY

    printf '%s\\t%s\\t%s\\t%s\\n' "${meta.id}" "\$CHEMISTRY" "\$MEDIAN" "\$BASIS" > ${prefix}.chemistry.txt
    echo "PACBIO_CHEMISTRY: ${meta.id} -> \$CHEMISTRY (\$BASIS)" >&2

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        minimap2: \$(minimap2 --version)
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    CHEMISTRY=pacbio_clr
    export CHEMISTRY
    printf '%s\\tpacbio_clr\\tNA\\tstub\\n' "${meta.id}" > ${prefix}.chemistry.txt
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        minimap2: 2.24-r1122
    END_VERSIONS
    """
}
