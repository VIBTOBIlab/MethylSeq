process METHURATOR_PLOT {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::methurator=0.1.8"
    container 'biocontainers/methurator:0.1.8--pyhdfd78af_0'
    label 'error_ignore'

    input:
    tuple val(meta), path(summary_report)

    output:
    tuple val(meta), path("plots/*.html")  , emit: plots
    path "versions.yml"                           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    methurator plot --summary $summary_report -o .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        methurator: \$(echo \$(methurator --version 2>&1) | sed 's/^.*methurator Version: v//; s/Copyright.*\$//')
    END_VERSIONS
    """
}
