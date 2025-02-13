process SEQ_SATURATION {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::samtools=1.17"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.17--h00cdaf9_0' :
        'biocontainers/samtools:1.17--h00cdaf9_0' }"

    input:
    tuple val(meta), path(bam), val(percentage)

    output:
    tuple val(meta), path("*.bam"), val(percentage), emit: bam
    path "${meta.id}_${percentage}.readcount.csv", emit: csv
    path  "versions.yml"          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    if (percentage.toFloat() !=1 ) { 
        args += " samtools view -s ${percentage} -b ${bam} > ${meta.id}_downsampled_${percentage}.bam" }
    else { args += " ln -s ${bam} ${meta.id}_downsampled_${percentage}.bam"}
    if ("$bam" == "${prefix}.bam") error "Input and output names are the same, use \"task.ext.prefix\" to disambiguate!"
    """
    $args
    READ_COUNT=\$(samtools view -c ${meta.id}_downsampled_${percentage}.bam)
    echo "${meta.id},${percentage},\$READ_COUNT" > ${meta.id}_${percentage}.readcount.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}