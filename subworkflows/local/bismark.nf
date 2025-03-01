/*
 * bismark subworkflow
 */
include { BISMARK_ALIGN                                                               } from '../../modules/nf-core/bismark/align/main'
include { SAMTOOLS_SORT as SAMTOOLS_SORT_DEDUPLICATED                                 } from '../../modules/nf-core/samtools/sort/main'
include { SAMTOOLS_SORT as SAMTOOLS_SORT_ALIGNED                                      } from '../../modules/nf-core/samtools/sort/main'
include { SAMTOOLS_QUERYNAME_SORT                                                     } from '../../modules/nf-core/samtools/queryname_sort/main'
include { SAMTOOLS_INDEX                                                              } from '../../modules/nf-core/samtools/index/main'
include { PICARD_MARKDUPLICATES                                                       } from '../../modules/nf-core/picard/markduplicates/main'
include { BISMARK_DEDUPLICATE                                                         } from '../../modules/nf-core/bismark/deduplicate/main'
include { SEQ_SATURATION_METHEXTRACTOR                                                } from '../../modules/local/seqsaturation_methextractor/main'
include { BISMARK_METHYLATIONEXTRACTOR                                                } from '../../modules/nf-core/bismark/methylationextractor/main'
include { BISMARK_COVERAGE2CYTOSINE                                                   } from '../../modules/nf-core/bismark/coverage2cytosine/main'
include { BISMARK_REPORT                                                              } from '../../modules/nf-core/bismark/report/main'
include { BISMARK_SUMMARY                                                             } from '../../modules/nf-core/bismark/summary/main'
include { BISMARK_FILTER_NON_CONVERSION                                               } from '../../modules/nf-core/bismark/filter_non_conversion/main'
include { SEQ_SATURATION                                                              } from '../../modules/local/seqsaturation/main'
include { PLOT_SEQ_SATURATION                                                         } from '../../modules/local/plot_seqsaturation/main'


workflow BISMARK {
    take:
    reads              // channel: [ val(meta), [ reads ] ]
    bismark_index      // channel: /path/to/BismarkIndex/
    skip_deduplication // boolean: whether to deduplicate alignments
    cytosine_report    // boolean: whether the run coverage2cytosine
    fasta              // channel: /path/to/fasta
    fasta_index        // channel: /path/to/fasta_index

    main:
    versions = Channel.empty()
    picard_metrics = Channel.empty()


    /*
     * Align with bismark
     */
    BISMARK_ALIGN (
        reads,
        bismark_index
    )
    alignments = BISMARK_ALIGN.out.bam
    versions = versions.mix(BISMARK_ALIGN.out.versions)


    /*
     * If seq saturation curve specified, it will generate the 
     * necessary files and plot
     */
    if (params.sequencing_curve | params.rrbs) {

        percentages_ch = Channel.fromList(params.downsampling_percentages.split(",").toList())
        downsample_input = BISMARK_ALIGN.out.bam.combine(percentages_ch)

        SEQ_SATURATION(
            downsample_input
        )
        versions = versions.mix(SEQ_SATURATION.out.versions)

        SEQ_SATURATION_METHEXTRACTOR(
            SEQ_SATURATION.out.bam,
            bismark_index

        )
        versions = versions.mix(SEQ_SATURATION_METHEXTRACTOR.out.versions)

        bam_res = SEQ_SATURATION.out.csv.
            collectFile(name: 'downsampling_reads_results.csv',
                        storeDir: "${params.outdir}/${params.aligner}/sequencing_saturation_curve/")
        cov_res = SEQ_SATURATION_METHEXTRACTOR.out.csv.
            collectFile(name: 'downsampling_cpgs_results.csv',
                        storeDir: "${params.outdir}/${params.aligner}/sequencing_saturation_curve/")

        PLOT_SEQ_SATURATION(
            bam_res,
            cov_res
        )
        
    }


    /*
     * If filter_non_conversion flag has been specified, 
     * run this step (filtering out non-bisulfite converted reads)
     */
    if (params.filter_non_conversion) {
        BISMARK_FILTER_NON_CONVERSION (
            alignments
        )
        alignments = BISMARK_FILTER_NON_CONVERSION.out.filter_bam
        versions = versions.mix(BISMARK_FILTER_NON_CONVERSION.out.versions)
    }

    /*
     * Sort the raw, aligned BAM file
     */
    SAMTOOLS_SORT_ALIGNED (
        BISMARK_ALIGN.out.bam
    )
    //alignments = SAMTOOLS_SORT_ALIGNED.out.bam
    versions = versions.mix(SAMTOOLS_SORT_ALIGNED.out.versions)

    /*
     * Remove optical duplicates if specified
     */
    if (params.remove_optic_duplicates) {
        // Sort by queryname for Picard module
        SAMTOOLS_QUERYNAME_SORT (
            alignments
        )
        PICARD_MARKDUPLICATES (
            alignments,
            fasta,
            fasta_index
        )
        alignments = PICARD_MARKDUPLICATES.out.bam
        picard_metrics = PICARD_MARKDUPLICATES.out.metrics
        versions = versions.mix(PICARD_MARKDUPLICATES.out.versions)
    }

    if (skip_deduplication) {
        alignment_reports = BISMARK_ALIGN.out.report.map{ meta, report -> [ meta, report, [] ] }
    } else {
        /*
        * Run deduplicate_bismark
        */
        BISMARK_DEDUPLICATE( alignments )
        alignments = BISMARK_DEDUPLICATE.out.bam
        alignment_reports = BISMARK_ALIGN.out.report.join(BISMARK_DEDUPLICATE.out.report)
        versions = versions.mix(BISMARK_DEDUPLICATE.out.versions)
    }

    /*
     * Run bismark_methylation_extractor
     */
    BISMARK_METHYLATIONEXTRACTOR (
        alignments,
        bismark_index
    )
    versions = versions.mix(BISMARK_METHYLATIONEXTRACTOR.out.versions)

    /*
     * Run coverage2cytosine
     */
    if (cytosine_report) {
        BISMARK_COVERAGE2CYTOSINE (
            BISMARK_METHYLATIONEXTRACTOR.out.coverage,
            bismark_index
        )
        versions = versions.mix(BISMARK_COVERAGE2CYTOSINE.out.versions)
    }

    /*
     * Generate bismark sample reports
     */
    BISMARK_REPORT (
        alignment_reports
            .join(BISMARK_METHYLATIONEXTRACTOR.out.report)
            .join(BISMARK_METHYLATIONEXTRACTOR.out.mbias)
    )
    versions = versions.mix(BISMARK_REPORT.out.versions)

    /*
     * Generate bismark summary report
     */
    BISMARK_SUMMARY (
        BISMARK_ALIGN.out.bam.collect{ it[1].name }.ifEmpty([]),
        alignment_reports.collect{ it[1] }.ifEmpty([]),
        alignment_reports.collect{ it[2] }.ifEmpty([]),
        BISMARK_METHYLATIONEXTRACTOR.out.report.collect{ it[1] }.ifEmpty([]),
        BISMARK_METHYLATIONEXTRACTOR.out.mbias.collect{ it[1] }.ifEmpty([])
    )
    versions = versions.mix(BISMARK_SUMMARY.out.versions)

    /*
     * Sort deduplicated output BAM
     */
    SAMTOOLS_SORT_DEDUPLICATED (
        alignments
    )
    versions = versions.mix(SAMTOOLS_SORT_DEDUPLICATED.out.versions)

    /*
     * Run samtools index on deduplicated alignment
     */
    SAMTOOLS_INDEX (SAMTOOLS_SORT_DEDUPLICATED.out.bam)
    versions = versions.mix(SAMTOOLS_INDEX.out.versions)

    /*
     * Collect MultiQC inputs
     */
    BISMARK_SUMMARY.out.summary.ifEmpty([])
        .mix(picard_metrics.collect{ it[1] })
        .mix(alignment_reports.collect{ it[1] })
        .mix(alignment_reports.collect{ it[2] })
        .mix(BISMARK_METHYLATIONEXTRACTOR.out.report.collect{ it[1] })
        .mix(BISMARK_METHYLATIONEXTRACTOR.out.mbias.collect{ it[1] })
        .mix(BISMARK_REPORT.out.report.collect{ it[1] })
        .set{ multiqc_files }

    emit:
    bam        = SAMTOOLS_SORT_ALIGNED.out.bam            // channel: [ val(meta), [ bam ] ] ## sorted, non-deduplicated (raw) BAM from aligner
    dedup      = SAMTOOLS_SORT_DEDUPLICATED.out.bam       // channel: [ val(meta), [ bam ] ] ## sorted, possibly deduplicated BAM
    mqc        = multiqc_files                            // path: *{html,txt}
    versions                                              // path: *.version.txt
}
