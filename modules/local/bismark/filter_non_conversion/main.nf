process BISMARK_FILTER_NON_CONVERSION {
    tag "$meta.id"
    label 'process_high'

    conda "bioconda::bismark=0.25.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bismark:0.25.0--hdfd78af_0' :
        'biocontainers/bismark:0.25.0--hdfd78af_0' }"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.nonCG_filtered.bam"), emit: filter_bam
    path "*.nonCG_removed_seqs.bam"
    path "*.non-conversion_filtering.txt", emit: report
    path "versions.yml" , emit: versions

    script:
    def args = ''
    if (!meta.single_end) {
        args += '-p'
    }
    """
    filter_non_conversion \\
    $args \\
    --percentage_cutoff ${params.percentage_cutoff} \\
    --minimum_count ${params.minimum_count} \\
    ${bam} \\

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bismark: \$(echo \$(bismark -v 2>&1) | sed 's/^.*Bismark Version: v//; s/Copyright.*\$//')
    END_VERSIONS
    """
}
