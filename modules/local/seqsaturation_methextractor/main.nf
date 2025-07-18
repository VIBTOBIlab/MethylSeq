process SEQ_SATURATION_METHEXTRACTOR {
    tag "$meta.id"
    label 'process_high'

    conda "bioconda::bismark=0.24.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bismark:0.24.0--hdfd78af_0' :
        'biocontainers/bismark:0.24.0--hdfd78af_0' }"

    input:
    tuple val(meta), path(bam), val(percentage)
    path index

    output:
    tuple val(meta), path("*.cov.gz"), val(percentage)          , emit: coverage
    path "*.csv"                                                , emit: csv
    path "versions.yml"                                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    // Assign sensible numbers for multicore and buffer_size based on bismark docs
    if(!args.contains('--multicore') && task.cpus >= 6){
        args += " --multicore ${(task.cpus / 3) as int}"
    }
    // Only set buffer_size when there are more than 6.GB of memory available
    if(!args.contains('--buffer_size') && task.memory?.giga > 6){
        args += " --buffer_size ${task.memory.giga - 2}G"
    }

    def seqtype  = meta.single_end ? '-s' : '-p'
    """
    bismark_methylation_extractor \\
        $bam \\
        --bedGraph \\
        --counts \\
        --gzip \\
        $seqtype \\
        $args \\

    IFS=',' read -ra numbers <<< "${params.min_counts}"
    echo "${meta.id}_${percentage}.cpgcount.csv"
    for i in "\${numbers[@]}"; do
        CpG_COUNT=\$(zcat ${meta.id}_downsampled_${percentage}.bismark.cov.gz | awk -v OFS='\\t' -v i="\$i" '\$5 + \$6 >= i' | wc -l)
        echo "${meta.id},${percentage},\${i},\$CpG_COUNT" >> "${meta.id}_${percentage}.cpgcount.csv"
    done
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bismark: \$(echo \$(bismark -v 2>&1) | sed 's/^.*Bismark Version: v//; s/Copyright.*\$//')
    END_VERSIONS
    """
}
