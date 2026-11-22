process PICARD_MARKOPTICALDUPLICATES {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/picard:3.0.0--hdfd78af_1' :
        'biocontainers/picard:3.0.0--hdfd78af_1' }"

    input:
    tuple val(meta), path(bam)
    path fasta

    output:
    tuple val(meta), path("*.bam")        , emit: bam
    tuple val(meta), path("*.bai")        , optional:true, emit: bai
    tuple val(meta), path("*.metrics.txt"), emit: metrics
    tuple val("${task.process}"), val('picard'), eval("picard MarkDuplicates --version 2>&1 | sed -n 's/.*Version://p'"), topic: versions, emit: versions_picard

    script:
    def args = task.ext.args ?: ''
    def avail_mem = 3072
    if (!task.memory) {
        log.info '[Picard MarkDuplicates] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = (task.memory.mega*0.8).intValue()
    }
    def opt_dupl_dist = 100
    if (params.sequencer=="NovaSeq") {
        opt_dupl_dist = 12000
    } else if (params.sequencer=="HiSeq") {
        opt_dupl_dist = 2500
    }
    """
    picard \\
        -Xmx${avail_mem}M \\
        MarkDuplicates \\
        $args \\
        --INPUT $bam \\
        --OUTPUT ${meta.id}.markdup.bam \\
        --REFERENCE_SEQUENCE $fasta \\
        --METRICS_FILE ${meta.id}.MarkDuplicates.metrics.txt \\
        --OPTICAL_DUPLICATE_PIXEL_DISTANCE ${opt_dupl_dist} \\
        --ASSUME_SORT_ORDER queryname \\
        --READ_NAME_REGEX '[a-zA-Z0-9]+:[0-9]+:[a-zA-Z0-9]+:[0-9]+:([0-9]+):([0-9]+):([0-9]+)_[0-9]+:[a-zA-Z0-9]+:[0-9]+:[a-zA-Z0-9]+[+][a-zA-Z0-9]+' \\
        --COMPRESSION_LEVEL 5
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bam
    touch ${prefix}.bam.bai
    touch ${prefix}.MarkDuplicates.metrics.txt
    """
}
