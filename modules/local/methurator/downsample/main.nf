process METHURATOR_DOWNSAMPLE {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::methurator=0.1.8"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/methurator:0.1.8--pyhdfd78af_0' :
        'biocontainers/methurator:0.1.8--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(bam)
    path fasta

    output:
    tuple val(meta), path("methurator_*.yml")    , emit: summary_report
    path  "versions.yml"              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    methurator downsample \\
        $bam \\
        --fasta $fasta \\
        -mc ${params.min_counts} \\
        -ds ${params.downsampling_percentages} \\
        -@ ${task.cpus} \\
        -o .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        methurator: \$(echo \$(methurator --version 2>&1) | sed 's/^.*methurator Version: v//; s/Copyright.*\$//')
    END_VERSIONS
    """
}
