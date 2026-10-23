#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VIBTOBIlab/MethylSeq
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/VIBTOBIlab/MethylSeq
    Website: https://vibtobilab.github.io/MethylSeq/
----------------------------------------------------------------------------------------
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { METHYLSEQ               } from './workflows/methylseq'
include { PIPELINE_INITIALISATION } from './subworkflows/local/utils_nfcore_methylseq_pipeline'
include { PIPELINE_COMPLETION     } from './subworkflows/local/utils_nfcore_methylseq_pipeline'
include { getGenomeAttribute      } from './subworkflows/local/utils_nfcore_methylseq_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PIPELINE PARAMETERS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

params {

    // Input/output options
    // Path to comma-separated file containing information about the samples in the experiment.
    input: String?
    // The output directory where the results will be saved. You have to use absolute paths to storage on Cloud infrastructure.
    outdir: String?
    // Email address for completion summary.
    email: String?
    // MultiQC report title. Printed as page header, used for filename if not otherwise specified.
    multiqc_title: String?

    // Save intermediate files
    // Save reference(s) to results directory
    save_reference: Boolean
    // Save aligned intermediates to results directory
    save_align_intermeds: Boolean
    // Bismark only - Save unmapped reads to FastQ files
    unmapped: Boolean
    // Save trimmed reads to results directory.
    save_trimmed: Boolean

    // Reference genome options
    // Name of iGenomes reference.
    genome: String?
    // Path to FASTA genome file
    fasta: String?
    // Path to Fasta index file.
    fasta_index: String?
    // Path to a directory containing a Bismark reference index.
    bismark_index: String?
    // bwameth index filename base
    bwa_meth_index: String?

    // Alignment options
    // Alignment tool to use.
    aligner: String = 'bismark'
    // Output information for all cytosine contexts.
    comprehensive: Boolean

    // Special library types
    // Preset for working with PBAT libraries.
    pbat: Boolean
    // Turn on if dealing with MspI digested material.
    rrbs: Boolean
    // Run bismark in SLAM-seq mode.
    slamseq: Boolean
    // Preset for EM-seq libraries.
    em_seq: Boolean
    // Trimming preset for single-cell bisulfite libraries.
    single_cell: Boolean
    // Trimming preset for the Accel kit.
    accel: Boolean
    // Trimming preset for the CEGX bisulfite kit.
    cegx: Boolean
    // Trimming preset for the Epignome kit.
    epignome: Boolean
    // Trimming preset for the Zymo kit.
    zymo: Boolean

    // Adapter Trimming
    // Trim bases from the 5' end of read 1 (or single-end reads).
    clip_r1: Integer = 0
    // Trim bases from the 5' end of read 2 (paired-end only).
    clip_r2: Integer = 0
    // Trim bases from the 3' end of read 1 AFTER adapter/quality trimming.
    three_prime_clip_r1: Integer = 0
    // Trim bases from the 3' end of read 2 AFTER adapter/quality trimming
    three_prime_clip_r2: Integer = 0
    // Trim bases below this quality value from the 3' end of the read, ignoring high-quality G bases
    nextseq_trim: Integer = 0

    // Bismark options
    // Run alignment against all four possible strands.
    non_directional: Boolean
    // Output stranded cytosine report, following Bismark's bismark_methylation_extractor step.
    cytosine_report: Boolean
    // Turn on to relax stringency for alignment (set allowed penalty with --num_mismatches).
    relax_mismatches: Boolean
    // 0.6 will allow a penalty of bp * -0.6 - for 100bp reads (bismark default is 0.2)
    num_mismatches: Float = 0.6
    // Specify a minimum read coverage to report a methylation call
    meth_cutoff: Integer?
    // Ignore read 2 methylation when it overlaps read 1
    no_overlap: Boolean = true
    // Ignore methylation in first n bases of 5' end of R1
    ignore_r1: Integer = 0
    // Ignore methylation in first n bases of 5' end of R2
    ignore_r2: Integer = 2
    // Ignore methylation in last n bases of 3' end of R1
    ignore_3prime_r1: Integer = 0
    // Ignore methylation in last n bases of 3' end of R2
    ignore_3prime_r2: Integer = 2
    // Supply a .gtf file containing known splice sites (bismark_hisat only).
    known_splices: String?
    // Allow soft-clipping of reads (potentially useful for single-cell experiments).
    local_alignment: Boolean
    // The minimum insert size for valid paired-end alignments.
    minins: Integer?
    // The maximum insert size for valid paired-end alignments.
    maxins: Integer?
    // Sample is NOMe-seq or NMT-seq. Runs coverage2cytosine.
    nomeseq: Boolean
    // Filtering incomplete bisulfite conversion from Bismark BAM files.
    filter_non_conversion: Boolean
    // Instead of filtering on an absolute count of methylated cytosines in non-CG context this option allows you to define an overall percentage of methylation in non-CG context (both CHH and CHG) which, if reached or exceeded, results in the read or read pair being filtered out.
    percentage_cutoff: Integer = 90
    // At least this number of cytosines in non-CG context (CHH or CHG) have to be seen in a read (irrespective of their methylation state) before the --percentage_cutoff filter kicks in.
    minimum_count: Integer = 3

    // Optical duplicates removal options
    // Mark and remove optical duplicates using PicardMarkDuplicates.
    remove_optic_duplicates: Boolean = true
    // Type of sequencing machine used.
    sequencer: String = 'NovaSeq'

    // methurator options
    // Maximum extrapolation factor (t) value.
    t_max: Integer = 10
    // Minimum number of counts to define a new CpG
    minimum_coverage: String = '1'
    // Run methurator sequencing saturation curve (def. true).
    run_methurator: Boolean = true

    // bwa-meth options
    // Specify a minimum read coverage for MethylDackel to report a methylation call.
    min_depth: Integer = 0
    // MethylDackel - ignore SAM flags
    ignore_flags: Boolean
    // Save files for use with methylKit
    methyl_kit: Boolean

    // Qualimap Options
    // A GFF or BED file containing the target regions which will be passed to Qualimap/Bamqc.
    bamqc_regions_file: String?

    // Skip pipeline steps
    // Skip read trimming.
    skip_trimming: Boolean
    // Skip deduplication step after alignment.
    skip_deduplication: Boolean
    // Skip MultiQC
    skip_multiqc: Boolean
    // Skip preseq step.
    skip_preseq: Boolean

    // Generic options
    // Display the help message.
    help: String?
    // Display version and exit.
    version: Boolean
    // Email address for completion summary, only when pipeline fails.
    email_on_fail: String?
    // Send plain-text email instead of HTML.
    plaintext_email: Boolean
    // File size limit when attaching MultiQC reports to summary emails.
    max_multiqc_email_size: String = '25.MB'
    // Do not use coloured log outputs.
    monochrome_logs: Boolean
    // Custom config file to supply to MultiQC.
    multiqc_config: String?
    // Custom logo file to supply to MultiQC. File name must also be set in the MultiQC config file
    multiqc_logo: String?
    // Custom MultiQC yaml file containing HTML including a methods description.
    multiqc_methods_description: String?
    // Boolean whether to validate parameters against the schema at runtime
    validate_params: Boolean = true
    // Display the full detailed help message.
    help_full: Boolean
    // Display hidden parameters in the help message (only works when --help or --help_full are provided).
    show_hidden: Boolean
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOWS FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// WORKFLOW: Run main analysis pipeline
//
workflow VIBTOBILAB_METHYLSEQ {

    take:
    samplesheet    // channel: samplesheet read in from --input
    fasta          //  string: path to the genome FASTA file
    fasta_index    //  string: path to the genome FASTA index
    bismark_index  //  string: path to the Bismark index
    bwa_meth_index //  string: path to the bwa-meth index

    main:

    //
    // WORKFLOW: Run pipeline
    //
    METHYLSEQ (
        samplesheet,
        fasta,
        fasta_index,
        bismark_index,
        bwa_meth_index,
        params.multiqc_config,
        params.multiqc_logo,
        params.multiqc_methods_description,
        params.outdir,
    )

    emit:
    multiqc_report = METHYLSEQ.out.multiqc_report // channel: /path/to/multiqc_report.html
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {

    main:
    //
    // SUBWORKFLOW: Run initialisation tasks
    //
    PIPELINE_INITIALISATION (
        params.version,
        params.validate_params,
        params.monochrome_logs,
        args,
        params.outdir,
        params.input,
        params.help,
        params.help_full,
        params.show_hidden
    )

    //
    // WORKFLOW: Run main workflow
    //
    //
    // Fall back to the iGenomes references of --genome when not given explicitly
    //
    VIBTOBILAB_METHYLSEQ (
        PIPELINE_INITIALISATION.out.samplesheet,
        params.fasta          ?: getGenomeAttribute('fasta'),
        params.fasta_index    ?: getGenomeAttribute('fasta_index'),
        params.bismark_index  ?: getGenomeAttribute('bismark'),
        params.bwa_meth_index ?: getGenomeAttribute('bwa_meth')
    )

    //
    // SUBWORKFLOW: Run completion tasks
    //
    PIPELINE_COMPLETION (
        params.email,
        params.email_on_fail,
        params.plaintext_email,
        params.outdir,
        params.monochrome_logs,
        VIBTOBILAB_METHYLSEQ.out.multiqc_report
    )
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
