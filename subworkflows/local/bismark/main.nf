/*
 * bismark subworkflow
 */
include { BISMARK_ALIGN                                                               } from '../../../modules/nf-core/bismark/align/main'
include { SAMTOOLS_SORT as SAMTOOLS_SORT_DEDUPLICATED                                 } from '../../../modules/nf-core/samtools/sort/main'
include { SAMTOOLS_SORT as SAMTOOLS_SORT_ALIGNED                                      } from '../../../modules/nf-core/samtools/sort/main'
include { SAMTOOLS_QUERYNAME_SORT                                                     } from '../../../modules/local/samtools/queryname_sort/main'
include { SAMTOOLS_INDEX                                                              } from '../../../modules/nf-core/samtools/index/main'
include { PICARD_MARKOPTICALDUPLICATES                                                } from '../../../modules/local/picard/markopticalduplicates/main'
include { BISMARK_DEDUPLICATE                                                         } from '../../../modules/nf-core/bismark/deduplicate/main'
include { BISMARK_METHYLATIONEXTRACTOR                                                } from '../../../modules/nf-core/bismark/methylationextractor/main'
include { BISMARK_COVERAGE2CYTOSINE                                                   } from '../../../modules/nf-core/bismark/coverage2cytosine/main'
include { BISMARK_REPORT                                                              } from '../../../modules/nf-core/bismark/report/main'
include { BISMARK_SUMMARY                                                             } from '../../../modules/nf-core/bismark/summary/main'
include { BISMARK_FILTER_NON_CONVERSION                                               } from '../../../modules/local/bismark/filter_non_conversion/main'


workflow BISMARK {
    take:
    reads              // channel: [ val(meta), [ reads ] ]
    bismark_index      // channel: /path/to/BismarkIndex/
    skip_deduplication // boolean: whether to deduplicate alignments
    cytosine_report    // boolean: whether the run coverage2cytosine
    fasta              // channel: /path/to/fasta

    main:
    picard_metrics = channel.empty()

    def ch_fasta = fasta.map { fasta_file -> [ [:], fasta_file ] }
    def ch_bismark_index = bismark_index.map { index -> [ [:], index ] }

    /*
     * Align with bismark
     */
    BISMARK_ALIGN (
        reads,
        ch_fasta,
        ch_bismark_index
    )
    alignments = BISMARK_ALIGN.out.bam

    /*
     * If filter_non_conversion flag has been specified,
     * run this step (filtering out non-bisulfite converted reads)
     */
    if (params.filter_non_conversion) {
        BISMARK_FILTER_NON_CONVERSION (
            alignments
        )
        alignments = BISMARK_FILTER_NON_CONVERSION.out.filter_bam
    }

    /*
     * Sort the raw, aligned BAM file
     */
    SAMTOOLS_SORT_ALIGNED (
        BISMARK_ALIGN.out.bam,
        [ [:], [], [] ],
        ''
    )

    /*
     * Remove optical duplicates if specified
     */
    if (params.remove_optic_duplicates) {
        // Sort by queryname for Picard module
        SAMTOOLS_QUERYNAME_SORT (
            alignments
        )
        PICARD_MARKOPTICALDUPLICATES (
            alignments,
            fasta
        )
        alignments = PICARD_MARKOPTICALDUPLICATES.out.bam
        picard_metrics = PICARD_MARKOPTICALDUPLICATES.out.metrics
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
    }

    /*
     * Run bismark_methylation_extractor
     */
    BISMARK_METHYLATIONEXTRACTOR (
        alignments,
        ch_bismark_index
    )

    /*
     * Run coverage2cytosine
     */
    if (cytosine_report) {
        BISMARK_COVERAGE2CYTOSINE (
            BISMARK_METHYLATIONEXTRACTOR.out.coverage,
            ch_fasta,
            ch_bismark_index
        )
    }

    /*
     * Generate bismark sample reports
     */
    BISMARK_REPORT (
        alignment_reports
            .join(BISMARK_METHYLATIONEXTRACTOR.out.report)
            .join(BISMARK_METHYLATIONEXTRACTOR.out.mbias)
    )

    /*
     * Generate bismark summary report
     */
    BISMARK_SUMMARY (
        BISMARK_ALIGN.out.bam.collect{ _meta, bam -> bam.name }.ifEmpty([]),
        alignment_reports.collect{ _meta, align_report, _dedup_report -> align_report }.ifEmpty([]),
        alignment_reports.collect{ _meta, _align_report, dedup_report -> dedup_report }.ifEmpty([]),
        BISMARK_METHYLATIONEXTRACTOR.out.report.collect{ _meta, report -> report }.ifEmpty([]),
        BISMARK_METHYLATIONEXTRACTOR.out.mbias.collect{ _meta, mbias -> mbias }.ifEmpty([])
    )

    /*
     * Sort potentially deduplicated output BAM
     */
    SAMTOOLS_SORT_DEDUPLICATED (
        alignments,
        [ [:], [], [] ],
        ''
    )

    /*
     * Run samtools index on potentially deduplicated alignment
     */
    SAMTOOLS_INDEX (SAMTOOLS_SORT_DEDUPLICATED.out.bam)

    /*
     * Collect MultiQC inputs
     */
    BISMARK_SUMMARY.out.summary.ifEmpty([])
        .mix(picard_metrics.collect{ _meta, metrics -> metrics })
        .mix(alignment_reports.collect{ _meta, align_report, _dedup_report -> align_report })
        .mix(alignment_reports.collect{ _meta, _align_report, dedup_report -> dedup_report })
        .mix(BISMARK_METHYLATIONEXTRACTOR.out.report.collect{ _meta, report -> report })
        .mix(BISMARK_METHYLATIONEXTRACTOR.out.mbias.collect{ _meta, mbias -> mbias })
        .mix(BISMARK_REPORT.out.report.collect{ _meta, report -> report })
        .set{ multiqc_files }

    emit:
    bam        = SAMTOOLS_SORT_ALIGNED.out.bam            // channel: [ val(meta), [ bam ] ] ## sorted, non-deduplicated (raw) BAM from aligner
    dedup      = SAMTOOLS_SORT_DEDUPLICATED.out.bam       // channel: [ val(meta), [ bam ] ] ## sorted, **potentially** deduplicated BAM
    dedup_bai  = SAMTOOLS_INDEX.out.index                 // channel: [ val(meta), [ bai ] ] ## dedup bam index
    mqc        = multiqc_files                            // path: *{html,txt}
}
