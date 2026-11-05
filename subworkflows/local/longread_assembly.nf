include { NANOPLOT     } from '../../modules/local/nanoplot/main.nf'
include { FILTLONG     } from '../../modules/local/filtlong/main.nf'
include { FLYE as FLYE_ONT          } from '../../modules/local/flye/main.nf'
include { FLYE as FLYE_PACBIO_HIFI  } from '../../modules/local/flye/main.nf'
include { FLYE as FLYE_PACBIO_CLR   } from '../../modules/local/flye/main.nf'
include { PACBIO_CHEMISTRY          } from '../../modules/local/pacbio_chemistry/main.nf'

// QC + assembly for ONT/PacBio long-read samples. Once Flye produces contigs,
// they're tagged is_contig:true and emitted in the same shape as any other
// contig sample - VARIANT_CALLING_WF's existing Snippy --ctgs mode,
// ASSEMBLY_TYPING_AMR_WF's contig-skip-assembly branch, and CLUSTERING_WF all
// already handle is_contig:true samples without modification, so long reads
// join those paths "for free" once assembled here.
//
// FLYE's mode is a broadcast value, not per-sample, so ONT and PacBio
// samples are split and run as two separate FLYE calls (aliased FLYE_ONT/
// FLYE_PACBIO, same underlying module) rather than one mixed call - applying
// the wrong mode to the wrong platform would silently produce a bad assembly
// rather than an obvious error.
workflow LONGREAD_ASSEMBLY_WF {

    take:
        longread_ch   // channel: [ val(meta), [ reads ] ]  meta.platform: 'nanopore' | 'pacbio'

    main:
        ch_versions = Channel.empty()

        // NanoPlot/Filtlong take a single fastq, not a list - reads are always
        // single-end for long-read platforms (enforced in check_samplesheet.py)
        longread_single_ch = longread_ch.map { meta, reads -> [ meta, reads instanceof List ? reads[0] : reads ] }

        NANOPLOT ( longread_single_ch )
        ch_versions = ch_versions.mix(NANOPLOT.out.versions)

        FILTLONG ( longread_single_ch )
        ch_versions = ch_versions.mix(FILTLONG.out.versions)

        // Chemistry can be set per sample: platform 'pacbio_hifi' / 'pacbio_clr'
        // are explicit. Plain 'pacbio' is detected from the reads when
        // --pacbio_mode is 'auto' (the default), and otherwise falls back to the
        // run-wide --pacbio_mode. Detection is worth it because submitter labels,
        // instrument model, read length and reported quality all fail to predict
        // which Flye mode works, while read-to-read divergence separates the two
        // chemistries cleanly (see modules/local/pacbio_chemistry/main.nf).
        FILTLONG.out.reads
            .branch { meta, _reads ->
                detect: meta.platform == 'pacbio' && params.pacbio_mode.toString().toLowerCase() == 'auto'
                other:  true
            }
            .set { pre_ch }

        PACBIO_CHEMISTRY ( pre_ch.detect )
        ch_versions = ch_versions.mix(PACBIO_CHEMISTRY.out.versions)

        // Detected samples carry their chemistry as the third element; put it on
        // meta.platform so both routes below look identical from here on.
        detected_ch = PACBIO_CHEMISTRY.out.typed.map { meta, reads, chem -> [ meta + [platform: chem.trim()], reads ] }

        detected_ch.mix(pre_ch.other)
            .branch { meta, _reads ->
                def chem = meta.platform == 'pacbio' ? (params.pacbio_mode == 'hifi' ? 'pacbio_hifi' : 'pacbio_clr') : meta.platform
                ont:          chem == 'nanopore'
                pacbio_hifi:  chem == 'pacbio_hifi'
                pacbio_clr:   chem == 'pacbio_clr'
                unrecognized: true
            }
            .set { filtered_ch }

        // A sample matching none of the labels above would otherwise vanish
        // silently after NanoPlot/Filtlong - e.g. if a newer input stage sends
        // 'pacbio_hifi'/'pacbio_clr' samples to an older copy of this file that
        // only knows 'nanopore'/'pacbio', branch has no label for them and
        // drops them with no error. Fail loudly instead.
        filtered_ch.unrecognized.subscribe { meta, reads ->
            error "LONGREAD_ASSEMBLY_WF: sample '${meta.id}' has platform '${meta.platform}', which has no assembly route. Expected nanopore, pacbio, pacbio_hifi or pacbio_clr - check that all pipeline files are from the same version."
        }

        FLYE_ONT         ( filtered_ch.ont,         '--nano-raw' )
        FLYE_PACBIO_HIFI ( filtered_ch.pacbio_hifi, '--pacbio-hifi' )
        FLYE_PACBIO_CLR  ( filtered_ch.pacbio_clr,  '--pacbio-raw' )
        ch_versions = ch_versions.mix(FLYE_ONT.out.versions).mix(FLYE_PACBIO_HIFI.out.versions).mix(FLYE_PACBIO_CLR.out.versions)

        // Gunzip so output matches the plain (uncompressed) fasta convention
        // already used elsewhere for is_contig:true samples supplied directly
        // by the user (e.g. via --contigs_dir) - keeps every is_contig:true
        // sample in a consistent format regardless of how it got assembled,
        // rather than mixing gzipped and plain fasta under the same meta flag.
        assembled_ch = FLYE_ONT.out.fasta.mix(FLYE_PACBIO_HIFI.out.fasta).mix(FLYE_PACBIO_CLR.out.fasta)
        GUNZIP_LONGREAD_ASSEMBLY ( assembled_ch )

        contigs = GUNZIP_LONGREAD_ASSEMBLY.out.fasta
            .map { meta, fasta -> [ meta + [ is_contig: true ], [ fasta ] ] }

    emit:
        contigs             // channel: [ val(meta), [ fasta ] ]  meta.is_contig == true
        versions = ch_versions
}

// Uses ubuntu:22.04 (definitely-real, always-available base image) rather
// than a specific gzip-utility image, since all this needs is `zcat` - not
// worth pulling in another single-purpose container for one command.
process GUNZIP_LONGREAD_ASSEMBLY {
    tag "$meta.id"
    label 'process_single'
    container "docker.io/ubuntu:22.04"

    input:
    tuple val(meta), path(gz_fasta)

    output:
    tuple val(meta), path("*.fasta"), emit: fasta

    script:
    """
    zcat ${gz_fasta} > ${meta.id}.assembly.fasta
    """
}
