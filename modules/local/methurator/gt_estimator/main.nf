process METHURATOR_GTESTIMATOR {
    tag "$meta.id"
    label 'process_medium'
    errorStrategy 'ignore'

    conda "bioconda::methurator=2.1.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/methurator:2.1.1--pyhdfd78af_0' :
        'biocontainers/methurator:2.1.1--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(bam)
    tuple val(meta), path(bai)
    path fasta
    path fai

    output:
    tuple val(meta), path("methurator_summary.yml") , emit: summary_report
    path  "versions.yml"                            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    if (params.rrbs) {
        args += '--rrbs'
    }
    """
    methurator gt-estimator \\
        ${bam} \\
        ${args} \\
        --fasta ${fasta} \\
        --minimum-coverage ${params.minimum_coverage} \\
        --t-max ${params.t_max} \\
        -@ ${task.cpus} \\
        --outdir . \\
        --compute_ci

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        methurator: "\$(methurator --version --version | sed 's/.* //')"
    END_VERSIONS
    """
}
