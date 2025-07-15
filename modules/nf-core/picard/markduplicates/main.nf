process PICARD_MARKDUPLICATES {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::picard=3.0.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/picard:3.0.0--hdfd78af_1' :
        'biocontainers/picard:3.0.0--hdfd78af_1' }"

    input:
    tuple val(meta), path(bam)
    path fasta
    path fai

    output:
    tuple val(meta), path("*.bam")        , emit: bam
    tuple val(meta), path("*.bai")        , optional:true, emit: bai
    tuple val(meta), path("*.metrics.txt"), emit: metrics
    path  "versions.yml"                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def avail_mem = 3072
    if (!task.memory) {
        log.info '[Picard MarkDuplicates] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = (task.memory.mega*0.8).intValue()
    }
    if (params.sequencer=="NovaSeq") {
        opt_dupl_dist = 12000
    } else if (params.sequencer=="HiSeq") {
        opt_dupl_dist = 2500
    } else { opt_dupl_dist = 100 }
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
        --COMPRESSION_LEVEL 5 \\

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        picard: \$(echo \$(picard MarkDuplicates --version 2>&1) | grep -o 'Version:.*' | cut -f2- -d:)
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bam
    touch ${prefix}.bam.bai
    touch ${prefix}.MarkDuplicates.metrics.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        picard: \$(echo \$(picard MarkDuplicates --version 2>&1) | grep -o 'Version:.*' | cut -f2- -d:)
    END_VERSIONS
    """
}
