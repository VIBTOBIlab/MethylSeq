process BISMARK_FILTER_NON_CONVERSION {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/10/106c5858e0e1021f5eb80619e7ec56f02e5fc857d2bdbb00570d30a5add9b5e9/data' :
        'community.wave.seqera.io/library/bismark:3.4.0--4721154cfbeba5e3' }"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.nonCG_filtered.bam"), emit: filter_bam
    path "*.nonCG_removed_seqs.bam"
    path "*.non-conversion_filtering.txt", emit: report
    tuple val("${task.process}"), val('bismark'), eval("bismark --version 2>&1 | grep -Eo '[0-9]+\\.[0-9]+\\.[0-9]+'"), emit: versions_bismark, topic: versions

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
    ${bam}
    """
}
