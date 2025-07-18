process PLOT_SEQ_SATURATION {

    container "docker.io/egiuili/python3_pandas:v1"

    label 'process_low'    
    
    input:
    path reads
    path cpgs
    //val counts

    output:
    path ("*.png")

    script:
    """
    plot_reads_vs_cpgs.py \
    --cpgs_file ${cpgs} \
    --read_file ${reads} \
    --percentages ${params.downsampling_percentages} \
    """
}

