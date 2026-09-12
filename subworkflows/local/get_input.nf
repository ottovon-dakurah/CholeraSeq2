include { INPUT_CHECK  } from './input_check'
include { SRA_DOWNLOAD as SRA_DOWNLOAD_ILLUMINA } from '../../modules/local/sra_download/main.nf'
include { SRA_DOWNLOAD as SRA_DOWNLOAD_ONT      } from '../../modules/local/sra_download/main.nf'
include { SRA_DOWNLOAD as SRA_DOWNLOAD_PACBIO   } from '../../modules/local/sra_download/main.nf'

// Single entry point for getting input reads/contigs into the pipeline, in any
// of the four modes Cholera_genomics supported:
//   --input        a samplesheet CSV (existing CholeraSeq behaviour, unchanged)
//   --contigs_dir  directory of pre-assembled *_contigs.fasta files
//   --reads_dir    directory of fastq files - paired (_R1/_R2 or _1/_2 naming)
//                   or single-end (no such suffix)
//   --sra_list     text file, one SRA/ENA accession per line (Illumina)
// --platform (default 'illumina') tags every sample discovered via
// --reads_dir as 'illumina', 'nanopore', or 'pacbio' - one directory is
// assumed to be one platform per run.
//
// For accession lists specifically, THREE separate, platform-specific lists
// can be combined in a single run: --sra_list (Illumina), --ont_sra_list
// (Nanopore), --pacbio_sra_list (PacBio) - each downloaded independently via
// its own SRA_DOWNLOAD_<PLATFORM> alias (same underlying module, called
// three times - a process can be invoked more than once in one workflow).
// Any combination of the three can be supplied together in one run, e.g.
// --sra_list illumina.txt --ont_sra_list ont.txt. A single list can't mix
// platforms within itself - every accession in a given list is tagged with
// that list's platform uniformly.
//
// Samplesheet rows (--input) remain the most flexible option: they can set
// platform per-ROW via the optional 'platform' column, so a single
// samplesheet can mix Illumina/Nanopore/PacBio at the individual-sample
// level rather than one list per platform.
//
// --input is mutually exclusive with the other modes; --contigs_dir/
// --reads_dir/--sra_list/--ont_sra_list/--pacbio_sra_list can be freely
// combined with each other in a single run.
workflow GET_INPUT_WF {

    main:
        ch_versions = Channel.empty()
        reads_ch    = Channel.empty()

        def any_dir_or_list_mode = params.reads_dir || params.contigs_dir || params.sra_list || params.ont_sra_list || params.pacbio_sra_list

        if (params.input && any_dir_or_list_mode) {
            error("--input cannot be combined with --reads_dir/--contigs_dir/--sra_list/--ont_sra_list/--pacbio_sra_list. Use a samplesheet OR the directory/list flags, not both.")
        }

        if (params.input) {
            //
            // Samplesheet mode (existing behaviour, including sra_id column support)
            //
            INPUT_CHECK ( file(params.input) )
            reads_ch    = INPUT_CHECK.out.reads
            ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)

        } else if (any_dir_or_list_mode) {
            //
            // Directory / accession-list auto-discovery mode - no samplesheet
            //
            if (params.contigs_dir) {
                contigs_ch = Channel
                    .fromPath("${params.contigs_dir}/*_contigs.fasta", checkIfExists: false)
                    .map { fasta ->
                        def id = fasta.getBaseName().replaceAll(/_contigs$/, '')
                        [ [ id: id, single_end: true, is_contig: true, platform: 'na' ], [ fasta ] ]
                    }
                reads_ch = reads_ch.mix(contigs_ch)
            }

            if (params.reads_dir) {
                // Support both naming conventions: <id>_R1/_R2.fastq.gz and
                // <id>_1/_2.fastq.gz. If a directory happens to contain both
                // patterns for the same sample ID, it will appear twice - an
                // edge case left unhandled here rather than over-engineered.
                fastq_r_ch   = Channel.fromFilePairs("${params.reads_dir}/*_R{1,2}.fastq.gz", checkIfExists: false)
                fastq_num_ch = Channel.fromFilePairs("${params.reads_dir}/*_{1,2}.fastq.gz",   checkIfExists: false)
                fastq_pairs_ch = fastq_r_ch.mix(fastq_num_ch)
                    .map { id, pair -> [ [ id: id, single_end: false, is_contig: false, platform: params.platform ], pair ] }
                reads_ch = reads_ch.mix(fastq_pairs_ch)

                // Single-end: any *.fastq.gz file NOT matching the paired
                // suffix patterns above (e.g. <id>.fastq.gz with no _1/_2/_R1/_R2).
                // Filtered in Groovy rather than a separate glob, since a glob
                // for "everything except X" isn't directly expressible.
                fastq_se_ch = Channel
                    .fromPath("${params.reads_dir}/*.fastq.gz", checkIfExists: false)
                    .filter { fq -> !(fq.name ==~ /.*_R?[12]\.fastq\.gz$/) }
                    .map { fq ->
                        def id = fq.name.replaceAll(/\.fastq\.gz$/, '')
                        [ [ id: id, single_end: true, is_contig: false, platform: params.platform ], [ fq ] ]
                    }
                reads_ch = reads_ch.mix(fastq_se_ch)
            }

            //
            // Illumina accession list
            //
            if (params.sra_list) {
                illumina_accessions_ch = Channel
                    .fromPath(params.sra_list, checkIfExists: true)
                    .splitText()
                    .map { it.trim() }
                    .filter { it }

                SRA_DOWNLOAD_ILLUMINA (
                    illumina_accessions_ch.map { acc -> [ id: acc, sra_id: acc, platform: 'illumina' ] }
                )
                ch_versions = ch_versions.mix(SRA_DOWNLOAD_ILLUMINA.out.versions)

                illumina_sra_ch = SRA_DOWNLOAD_ILLUMINA.out.reads
                    .map { meta, reads ->
                        def reads_list = reads instanceof List ? reads : [ reads ]
                        def full_meta  = meta + [ single_end: reads_list.size() == 1, is_contig: false ]
                        [ full_meta, reads_list ]
                    }
                reads_ch = reads_ch.mix(illumina_sra_ch)
            }

            //
            // ONT (Nanopore) accession list
            //
            if (params.ont_sra_list) {
                ont_accessions_ch = Channel
                    .fromPath(params.ont_sra_list, checkIfExists: true)
                    .splitText()
                    .map { it.trim() }
                    .filter { it }

                SRA_DOWNLOAD_ONT (
                    ont_accessions_ch.map { acc -> [ id: acc, sra_id: acc, platform: 'nanopore' ] }
                )
                ch_versions = ch_versions.mix(SRA_DOWNLOAD_ONT.out.versions)

                ont_sra_ch = SRA_DOWNLOAD_ONT.out.reads
                    .map { meta, reads ->
                        def reads_list = reads instanceof List ? reads : [ reads ]
                        // ONT reads are always single-end - forced true here
                        // regardless of how many files were produced, since
                        // long-read platforms have no paired-end concept.
                        def full_meta  = meta + [ single_end: true, is_contig: false ]
                        [ full_meta, reads_list ]
                    }
                reads_ch = reads_ch.mix(ont_sra_ch)
            }

            //
            // PacBio accession list
            //
            if (params.pacbio_sra_list) {
                pacbio_accessions_ch = Channel
                    .fromPath(params.pacbio_sra_list, checkIfExists: true)
                    .splitText()
                    .map { it.trim() }
                    .filter { it }

                SRA_DOWNLOAD_PACBIO (
                    pacbio_accessions_ch.map { acc -> [ id: acc, sra_id: acc, platform: 'pacbio' ] }
                )
                ch_versions = ch_versions.mix(SRA_DOWNLOAD_PACBIO.out.versions)

                pacbio_sra_ch = SRA_DOWNLOAD_PACBIO.out.reads
                    .map { meta, reads ->
                        def reads_list = reads instanceof List ? reads : [ reads ]
                        def full_meta  = meta + [ single_end: true, is_contig: false ]
                        [ full_meta, reads_list ]
                    }
                reads_ch = reads_ch.mix(pacbio_sra_ch)
            }

        } else {
            error("No input specified. Use --input <samplesheet.csv>, or one or more of --reads_dir/--contigs_dir/--sra_list/--ont_sra_list/--pacbio_sra_list.")
        }

    emit:
        reads    = reads_ch   // channel: [ val(meta), [ reads ] ]
        versions = ch_versions
}
