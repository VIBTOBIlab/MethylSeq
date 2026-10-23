/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { PREPARE_GENOME         } from '../subworkflows/local/prepare_genome'
include { BISMARK                } from '../subworkflows/local/bismark'
include { BWAMETH                } from '../subworkflows/local/bwameth'

//
// MODULE: Installed directly from nf-core/modules
//
include { CAT_FASTQ              } from '../modules/nf-core/cat/fastq/main'
include { FASTQC                 } from '../modules/nf-core/fastqc/main'
include { MULTIQC                } from '../modules/nf-core/multiqc/main'
include { TRIMGALORE             } from '../modules/nf-core/trimgalore/main'
include { QUALIMAP_BAMQC         } from '../modules/nf-core/qualimap/bamqc/main'
include { PRESEQ_LCEXTRAP        } from '../modules/nf-core/preseq/lcextrap/main'
include { METHURATOR_GTESTIMATOR } from '../modules/nf-core/methurator/gtestimator/main'
include { METHURATOR_PLOT        } from '../modules/nf-core/methurator/plot/main'

//
// FUNCTIONS
//
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_methylseq_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow METHYLSEQ {

    take:
    ch_samplesheet // channel: [ val(meta), [ fastq ] ] read in from --input
    fasta          //  string: path to the genome FASTA file
    fasta_index    //  string: path to the genome FASTA index
    bismark_index  //  string: path to the Bismark index
    bwa_meth_index //  string: path to the bwa-meth index
    multiqc_config
    multiqc_logo
    multiqc_methods_description
    outdir

    main:

    def ch_versions = channel.empty()
    def ch_multiqc_files = channel.empty()

    //
    // SUBWORKFLOW: Prepare any required reference genome indices
    //
    PREPARE_GENOME(
        fasta,
        fasta_index,
        bismark_index,
        bwa_meth_index
    )

    //
    // Branch samples with one or multiple FastQ files per read
    //
    def ch_fastq = ch_samplesheet
        .branch { meta, fastq ->
            single: fastq.size() == 1
                return [ meta, fastq.flatten() ]
            multiple: fastq.size() > 1
                return [ meta, fastq.flatten() ]
        }

    //
    // MODULE: Concatenate FastQ files from same sample if required
    //
    CAT_FASTQ (
        ch_fastq.multiple
    )
    def ch_cat_fastq = CAT_FASTQ.out.reads.mix(ch_fastq.single)

    //
    // MODULE: Run FastQC
    //
    FASTQC (
        ch_cat_fastq
    )

    //
    // MODULE: Run TrimGalore!
    //
    def reads = ch_cat_fastq
    if (!params.skip_trimming) {
        TRIMGALORE(ch_cat_fastq)
        reads = TRIMGALORE.out.reads
    }

    //
    // SUBWORKFLOW: Align reads, deduplicate and extract methylation
    //
    def ch_bam = channel.empty()
    def ch_dedup = channel.empty()
    def ch_bai = channel.empty()
    def ch_aligner_mqc = channel.empty()

    // Aligner: bismark or bismark_hisat
    if (params.aligner =~ /bismark/) {

        BISMARK (
            reads,
            PREPARE_GENOME.out.bismark_index,
            params.skip_deduplication || params.rrbs,
            params.cytosine_report || params.nomeseq,
            PREPARE_GENOME.out.fasta
        )
        ch_bam = BISMARK.out.bam
        ch_dedup = BISMARK.out.dedup
        ch_bai = BISMARK.out.dedup_bai
        ch_aligner_mqc = BISMARK.out.mqc
    }
    // Aligner: bwameth
    else if (params.aligner == 'bwameth') {

        BWAMETH (
            reads,
            PREPARE_GENOME.out.bwameth_index,
            PREPARE_GENOME.out.fasta,
            PREPARE_GENOME.out.fasta_index,
            params.skip_deduplication || params.rrbs,
        )
        ch_bam = BWAMETH.out.bam
        ch_dedup = BWAMETH.out.dedup
        ch_bai = BWAMETH.out.dedup_bai
        ch_aligner_mqc = BWAMETH.out.mqc
    }

    //
    // MODULE: Qualimap BamQC
    //
    QUALIMAP_BAMQC (
        ch_dedup,
        params.bamqc_regions_file ? channel.fromPath( params.bamqc_regions_file, checkIfExists: true ).toList() : []
    )

    //
    // MODULE: Run Preseq
    //
    if (!params.skip_preseq && !params.rrbs) {
        PRESEQ_LCEXTRAP (
            ch_bam
        )
    }

    //
    // MODULE: Run methurator
    //
    if (params.run_methurator) {

        METHURATOR_GTESTIMATOR (
            ch_dedup
                .join(ch_bai)
                .combine(PREPARE_GENOME.out.fasta)
                .combine(PREPARE_GENOME.out.fasta_index)
        )

        METHURATOR_PLOT (
            METHURATOR_GTESTIMATOR.out.summary_report
        )
    }

    //
    // Collate and save software versions
    //
    def topic_versions = channel.topic("versions")
        .distinct()
        .branch { entry ->
            versions_file: entry instanceof Path
            versions_tuple: true
        }

    def topic_versions_string = topic_versions.versions_tuple
        .map { process, tool, version ->
            [ process[process.lastIndexOf(':')+1..-1], "  ${tool}: ${version}" ]
        }
        .groupTuple(by:0)
        .map { process, tool_versions ->
            tool_versions.unique().sort()
            "${process}:\n${tool_versions.join('\n')}"
        }

    def ch_collated_versions = softwareVersionsToYAML(ch_versions.mix(topic_versions.versions_file))
        .mix(topic_versions_string)
        .collectFile(
            storeDir: "${outdir}/pipeline_info",
            name:  'methylseq_software_'  + 'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        )

    //
    // MODULE: MultiQC
    //
    def ch_multiqc_report = channel.empty()
    if (!params.skip_multiqc) {
        ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
        def ch_summary_params = paramsSummaryMap(parameters_schema: "nextflow_schema.json")
        def ch_workflow_summary = channel.value(paramsSummaryMultiqc(ch_summary_params))
        ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
        def ch_multiqc_custom_methods_description = multiqc_methods_description
            ? file(multiqc_methods_description, checkIfExists: true)
            : file("${projectDir}/assets/methods_description_template.yml", checkIfExists: true)
        def ch_methods_description = channel.value(methodsDescriptionText(ch_multiqc_custom_methods_description))
        ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml', sort: true))
        ch_multiqc_files = ch_multiqc_files.mix(QUALIMAP_BAMQC.out.results.collect{ _meta, results -> results }.ifEmpty([]))
        ch_multiqc_files = ch_multiqc_files.mix(ch_aligner_mqc.ifEmpty([]))
        if (!params.skip_trimming) {
            ch_multiqc_files = ch_multiqc_files.mix(TRIMGALORE.out.log.collect{ _meta, log -> log })
        }
        if (!params.skip_preseq && !params.rrbs) {
            ch_multiqc_files = ch_multiqc_files.mix(PRESEQ_LCEXTRAP.out.log.collect{ _meta, log -> log }.ifEmpty([]))
        }
        ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{ _meta, zip -> zip }.ifEmpty([]))

        // The pipeline MultiQC config is always used, a custom config is applied on top of it
        def multiqc_configs = [ file("${projectDir}/assets/multiqc_config.yml", checkIfExists: true) ]
        if (multiqc_config) {
            multiqc_configs << file(multiqc_config, checkIfExists: true)
        }

        MULTIQC(
            ch_multiqc_files.flatten().collect().map { files ->
                [
                    [id: 'methylseq'],
                    files,
                    multiqc_configs,
                    multiqc_logo ? file(multiqc_logo, checkIfExists: true) : [],
                    [],
                    [],
                ]
            }
        )
        ch_multiqc_report = MULTIQC.out.report.map { _meta, report -> [report] }.toList()
    }

    emit:
    multiqc_report = ch_multiqc_report // channel: /path/to/multiqc_report.html
    versions       = ch_versions       // channel: [ path(versions.yml) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
