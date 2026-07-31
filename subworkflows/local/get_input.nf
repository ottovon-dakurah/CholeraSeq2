include { INPUT_CHECK  } from './input_check'
include { SRA_DOWNLOAD } from '../../modules/local/sra_download/main.nf'

// Single entry point for getting input reads/contigs into the pipeline, in any
// of the four modes Cholera_genomics supported:
//   --input        a samplesheet CSV (existing CholeraSeq behaviour, unchanged)
//   --contigs_dir  directory of pre-assembled *_contigs.fasta files
//   --reads_dir    directory of fastq files - paired (_R1/_R2 or _1/_2 naming)
//                   or single-end (no such suffix)
//   --sra_list     text file, one SRA/ENA accession per line
// --input is mutually exclusive with the other three; --contigs_dir/--reads_dir/
// --sra_list can be freely combined with each other in a single run.
workflow GET_INPUT_WF {

    main:
        ch_versions = Channel.empty()
        reads_ch    = Channel.empty()

        if (params.input && (params.reads_dir || params.contigs_dir || params.sra_list)) {
            error("--input cannot be combined with --reads_dir/--contigs_dir/--sra_list. Use a samplesheet OR the directory/list flags, not both.")
        }

        if (params.input) {
            //
            // Samplesheet mode (existing behaviour, including sra_id column support)
            //
            INPUT_CHECK ( file(params.input) )
            reads_ch    = INPUT_CHECK.out.reads
            ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)

        } else if (params.reads_dir || params.contigs_dir || params.sra_list) {
            //
            // Directory / accession-list auto-discovery mode - no samplesheet
            //
            if (params.contigs_dir) {
                contigs_ch = Channel
                    .fromPath("${params.contigs_dir}/*_contigs.fasta", checkIfExists: false)
                    .map { fasta ->
                        def id = fasta.getBaseName().replaceAll(/_contigs$/, '')
                        [ [ id: id, single_end: true, is_contig: true ], [ fasta ] ]
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
                    .map { id, pair -> [ [ id: id, single_end: false, is_contig: false ], pair ] }
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
                        [ [ id: id, single_end: true, is_contig: false ], [ fq ] ]
                    }
                reads_ch = reads_ch.mix(fastq_se_ch)
            }

            if (params.sra_list) {
                sra_accessions_ch = Channel
                    .fromPath(params.sra_list, checkIfExists: true)
                    .splitText()
                    .map { it.trim() }
                    .filter { it }

                SRA_DOWNLOAD (
                    sra_accessions_ch.map { acc -> [ id: acc, sra_id: acc ] }
                )
                ch_versions = ch_versions.mix(SRA_DOWNLOAD.out.versions)

                sra_ch = SRA_DOWNLOAD.out.reads
                    .map { meta, reads ->
                        def reads_list = reads instanceof List ? reads : [ reads ]
                        def full_meta  = meta + [ single_end: reads_list.size() == 1, is_contig: false ]
                        [ full_meta, reads_list ]
                    }
                reads_ch = reads_ch.mix(sra_ch)
            }

        } else {
            error("No input specified. Use --input <samplesheet.csv>, or one or more of --reads_dir/--contigs_dir/--sra_list.")
        }

    emit:
        reads    = reads_ch   // channel: [ val(meta), [ reads ] ]
        versions = ch_versions
}
