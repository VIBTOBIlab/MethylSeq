process METHURATOR_GTESTIMATOR {
    tag "$meta.id"
    label 'process_medium'
    errorStrategy 'ignore'
    conda "bioconda::methurator=2.0.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/methurator:2.0.0--pyhdfd78af_0' :
        'biocontainers/methurator:2.0.0--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(bam)
    path fasta

    output:
    tuple val(meta), path("methurator_*.yml")    , emit: summary_report
    path  "versions.yml"              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = ''
    if(params.rrbs) { args += '--rrbs'}
    """
    methurator gt-estimator \\
        $bam \\
        --fasta $fasta \\
        -mc ${params.min_counts} \\
        --t-max ${params.tmax} \\
        -@ ${task.cpus} \\
        --compute_ci \\
        -o . $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        methurator: "\$(methurator --version 2>&1 | sed -E 's/.*version[[:space:]]+([0-9.]+).*/\\1/')"
    END_VERSIONS
    """
}
