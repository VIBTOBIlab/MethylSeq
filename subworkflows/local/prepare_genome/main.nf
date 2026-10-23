//
// Prepare reference genome files
//

include { UNTAR                       } from '../../../modules/nf-core/untar/main'
include { BISMARK_GENOMEPREPARATION   } from '../../../modules/nf-core/bismark/genomepreparation/main'
include { BWAMETH_INDEX               } from '../../../modules/nf-core/bwameth/index/main'
include { SAMTOOLS_FAIDX              } from '../../../modules/nf-core/samtools/faidx/main'

workflow PREPARE_GENOME {

    take:
    fasta          // string: path to the genome FASTA file
    fasta_index    // string: path to the genome FASTA index
    bismark_index  // string: path to the Bismark index
    bwa_meth_index // string: path to the bwa-meth index

    main:
    ch_fasta         = channel.empty()
    ch_bismark_index = channel.empty()
    ch_bwameth_index = channel.empty()
    ch_fasta_index   = channel.empty()

    // FASTA, if supplied or from genome
    if (fasta) {
        ch_fasta = channel.value(file(fasta))
    }

    // Aligner: bismark or bismark_hisat
    if( params.aligner =~ /bismark/ ){

        /*
         * Generate bismark index if not supplied
         */
        if (bismark_index) {
            if (bismark_index.endsWith('.gz')) {
                ch_bismark_index = UNTAR ( [ [:], file(bismark_index) ] ).untar.map { _meta, index -> index }
            } else {
                ch_bismark_index = channel.value(file(bismark_index))
            }
        } else {
            BISMARK_GENOMEPREPARATION(ch_fasta.map { fasta_file -> [ [:], fasta_file ] })
            ch_bismark_index = BISMARK_GENOMEPREPARATION.out.index.map { _meta, index -> index }
        }

    }
    // Aligner: bwameth
    else if ( params.aligner == 'bwameth' ){

        /*
         * Generate bwameth index if not supplied
         */
        if (bwa_meth_index) {
            if (bwa_meth_index.endsWith('.tar.gz')) {
                ch_bwameth_index = UNTAR ( [ [:], file(bwa_meth_index) ] ).untar.map { _meta, index -> index }
            } else {
                ch_bwameth_index = channel.value(file(bwa_meth_index))
            }
        } else {
            BWAMETH_INDEX(ch_fasta.map { fasta_file -> [ [:], fasta_file ] }, false)
            ch_bwameth_index = BWAMETH_INDEX.out.index.map { _meta, index -> index }
        }
    }

    /*
     * Generate fasta index if not supplied
     */
    if (fasta_index) {
        ch_fasta_index = channel.value(file(fasta_index))
    } else {
        SAMTOOLS_FAIDX(ch_fasta.map { fasta_file -> [ [:], fasta_file, [] ] }, false)
        ch_fasta_index = SAMTOOLS_FAIDX.out.fai.map { _meta, fai -> fai }
    }

    emit:
    fasta         = ch_fasta                  // channel: path(genome.fasta)
    bismark_index = ch_bismark_index          // channel: path(genome.fasta)
    bwameth_index = ch_bwameth_index          // channel: path(genome.fasta)
    fasta_index   = ch_fasta_index            // channel: path(genome.fasta)

}
