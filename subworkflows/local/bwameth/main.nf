/*
 * bwameth subworkflow
 */
include { SAMTOOLS_STATS                                } from '../../../modules/nf-core/samtools/stats/main'
include { SAMTOOLS_FLAGSTAT                             } from '../../../modules/nf-core/samtools/flagstat/main'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_ALIGNMENTS   } from '../../../modules/nf-core/samtools/index/main'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_DEDUPLICATED } from '../../../modules/nf-core/samtools/index/main'
include { SAMTOOLS_SORT                                 } from '../../../modules/nf-core/samtools/sort/main'
include { PICARD_MARKDUPLICATES                         } from '../../../modules/nf-core/picard/markduplicates/main'
include { BWAMETH_ALIGN                                 } from '../../../modules/nf-core/bwameth/align/main'
include { METHYLDACKEL_EXTRACT                          } from '../../../modules/nf-core/methyldackel/extract/main'
include { METHYLDACKEL_MBIAS                            } from '../../../modules/nf-core/methyldackel/mbias/main'

workflow BWAMETH {
    take:
    reads               // channel: [ val(meta), [ reads ] ]
    bwameth_index       // channel: /path/to/bwa_index/
    fasta               // channel: /path/to/genome.fa
    fasta_index         // channel: /path/to/genome.fa.fai
    skip_deduplication  // boolean: whether to deduplicate alignments

    main:
    def ch_fasta = fasta.map { fasta_file -> [ [:], fasta_file ] }
    // combine() returns a queue channel: first() turns it back into a value channel used by every sample
    def ch_fasta_fai = fasta.combine(fasta_index).map { fasta_file, fai -> [ [:], fasta_file, fai ] }.first()

    /*
     * Align with bwameth
     */
    BWAMETH_ALIGN (
        reads,
        ch_fasta,
        bwameth_index.map { index -> [ [:], index ] }
    )

    /*
     * Sort raw output BAM
     */
    SAMTOOLS_SORT (
        BWAMETH_ALIGN.out.bam,
        [ [:], [], [] ],
        ''
    )

    /*
     * Run samtools index on alignment
     */
    SAMTOOLS_INDEX_ALIGNMENTS (SAMTOOLS_SORT.out.bam)

    /*
     * Run samtools flagstat and samtools stats
     */
    SAMTOOLS_FLAGSTAT( BWAMETH_ALIGN.out.bam.join(SAMTOOLS_INDEX_ALIGNMENTS.out.index) )
    SAMTOOLS_STATS( BWAMETH_ALIGN.out.bam.join(SAMTOOLS_INDEX_ALIGNMENTS.out.index), [ [:], [], [] ] )

    if (skip_deduplication) {
        alignments = SAMTOOLS_SORT.out.bam
        bam_index = SAMTOOLS_INDEX_ALIGNMENTS.out.index
        picard_metrics = channel.empty()
    } else {
        /*
        * Run Picard MarkDuplicates
        */
        PICARD_MARKDUPLICATES (
            SAMTOOLS_SORT.out.bam,
            ch_fasta_fai
        )
        /*
         * Run samtools index on deduplicated alignment
        */
        SAMTOOLS_INDEX_DEDUPLICATED (PICARD_MARKDUPLICATES.out.bam)

        alignments = PICARD_MARKDUPLICATES.out.bam
        bam_index = SAMTOOLS_INDEX_DEDUPLICATED.out.index
        picard_metrics = PICARD_MARKDUPLICATES.out.metrics
    }


    /*
     * Extract per-base methylation and plot methylation bias
     */

    METHYLDACKEL_EXTRACT(
        alignments.join(bam_index),
        ch_fasta_fai
    )
    METHYLDACKEL_MBIAS(
        alignments.join(bam_index),
        ch_fasta_fai
    )

    /*
     * Collect MultiQC inputs
     */
    picard_metrics.collect{ _meta, metrics -> metrics }
        .mix(SAMTOOLS_FLAGSTAT.out.flagstat.collect{ _meta, flagstat -> flagstat })
        .mix(SAMTOOLS_STATS.out.stats.collect{ _meta, stats -> stats })
        .mix(METHYLDACKEL_EXTRACT.out.bedgraph.collect{ _meta, bedgraph -> bedgraph })
        .mix(METHYLDACKEL_MBIAS.out.txt.collect{ _meta, txt -> txt })
        .set{ multiqc_files }

    emit:
    bam                  = SAMTOOLS_SORT.out.bam  // channel: [ val(meta), [ bam ] ] ## sorted, non-deduplicated (raw) BAM from aligner
    dedup                = alignments             // channel: [ val(meta), [ bam ] ] ## sorted, **potentially** deduplicated BAM
    dedup_bai            = bam_index              // channel: [ val(meta), [ bai ] ] ## dedup bam index
    mqc                  = multiqc_files          // path: *{html,txt}
}
