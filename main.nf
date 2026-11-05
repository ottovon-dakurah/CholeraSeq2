#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CERI-KRISP/CholeraSeq
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/CERI-KRISP/CholeraSeq----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOW FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PARAMETER TYPES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Values given on the command line arrive as text: without a declared type,
// `--max_cpus 8` is the string "8" (so arithmetic on it fails) and
// `--skip_fastbaps false` is the string "false", which is truthy and therefore
// does NOT switch the step off. Declaring the type here makes Nextflow convert
// command-line values, so numbers and booleans can be set with plain flags
// instead of only through a params file.
//
// Types only: every default stays in nextflow.config, which remains the single
// source of truth. String parameters need no declaration, as text is already
// what they expect. Decimals must be Float; Double and BigDecimal reject
// command-line values.
params {
    // switches
    use_global_alignment: Boolean
    skip_fastbaps: Boolean
    skip_clustering: Boolean
    skip_assembly: Boolean
    skip_mlst: Boolean
    skip_amr: Boolean
    plaintext_email: Boolean
    monochrome_logs: Boolean
    help: Boolean
    version: Boolean
    validate_params: Boolean
    show_hidden_params: Boolean

    // whole numbers
    max_cpus: Integer
    pacbio_detect_reads: Integer
    long_read_min_length: Integer
    long_read_keep_percent: Integer
    min_trim_quality: Integer
    min_trim_length: Integer
    min_mapping_quality: Integer
    min_base_quality: Integer
    min_site_coverage: Integer
    snippy_minqual: Integer
    max_missing_percentage: Integer
    abricate_minid: Integer
    abricate_mincov: Integer
    min_partition_size: Integer

    // decimals
    pacbio_detect_threshold: Float
    min_allele_fraction: Float
    min_parsimony_coverage: Float
}

include { validateParameters; paramsHelp; paramsSummaryLog; paramsSummaryMap } from 'plugin/nf-schema'
include { CHOLERASEQ } from './workflows/choleraseq'

//
// WORKFLOW: Run main CERI-KRISP/CholeraSeq analysis pipeline
//
workflow CERI_KRISP {
    CHOLERASEQ ()
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN ALL WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// WORKFLOW: Execute a single named workflow for the pipeline
// See: https://github.com/nf-core/rnaseq/issues/619
//
workflow {
    // VALIDATE & PRINT PARAMETER SUMMARY
    // NOTE: moved here from top-level script scope - bare statements outside
    // a process/workflow/function are not allowed under Nextflow's strict
    // script syntax (default since Nextflow 26.04).
    // Help, parameter summary and validation now come from nf-schema, driven by
    // nextflow_schema.json. --help is handled by the plugin itself (validation.help).
    if (params.help) {
        log.info WorkflowMain.logoAndCitation(workflow, params)
        log.info paramsHelp()
        System.exit(0)
    }
    log.info WorkflowMain.logoAndCitation(workflow, params)
    log.info paramsSummaryLog(workflow)
    validateParameters()

    // Remaining checks that are not expressible in the schema
    WorkflowMain.initialise(workflow, params, log)

    CERI_KRISP ()

    // COMPLETION EMAIL AND SUMMARY
    // NOTE: this must live at this exact nesting level (main.nf's entry workflow),
    // not inside CHOLERASEQ - registering workflow.onComplete deeper than this
    // throws a NullPointerException at runtime because `params` doesn't reliably
    // resolve inside a closure invoked later (at session shutdown) from that far
    // down the workflow-call chain. summary_params is recomputed here rather than
    // forwarded from CHOLERASEQ, since it only depends on workflow/params, both of
    // which are already known-good at this level.
    //
    // multiqc_report is built as a direct file path rather than forwarded through
    // CHOLERASEQ/CERI_KRISP's `emit:` - attempting that hits an unresolved Nextflow
    // bug ("Missing workflow output parameter") when a workflow output is nested
    // two `emit:` hops deep (see nextflow-io/nextflow#6204 for the same failure
    // mode). Since this value is only used to attach the report to a completion
    // email (a path only reached if params.email/email_on_fail is set), and
    // MULTIQC's publishDir is deterministic, constructing the path directly avoids
    // the bug entirely with no loss of functionality.
    //
    // `wf`/`prms`/`lg`/`pd` capture the corresponding implicit bindings BEFORE
    // entering the closure below. Referencing bare implicit identifiers (workflow,
    // params, log, projectDir) *inside* `workflow.onComplete {}` is ambiguous -
    // Groovy sets the closure's delegate to the workflow object itself, so those
    // names can resolve against that delegate instead of the outer script binding,
    // silently returning null. This affects every implicit binding used inside the
    // closure, not just `workflow` (confirmed: first hit workflow.revision, then
    // hit params.containsKey() the same way after aliasing only `workflow`) -
    // aliasing all four up front avoids fixing this one variable at a time.
    def wf   = workflow
    def prms = params
    def lg   = log
    def pd   = projectDir
    workflow.onComplete {
        def summary_params  = paramsSummaryMap(wf)
        def multiqc_report  = ["${prms.outdir}/multiqc/multiqc_report.html"]
        if (prms.email || prms.email_on_fail) {
            NfcoreTemplate.email(wf, prms, summary_params, pd, lg, multiqc_report)
        }
        NfcoreTemplate.summary(wf, prms, lg)
        if (prms.hook_url) {
            NfcoreTemplate.IM_notification(wf, prms, summary_params, pd, lg)
        }
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
