process bam2wig {

    tag "${id}"
    label 'big_cpu'

    publishDir( 
        "${params.outputs}/bigwig", 
        mode: 'copy',
    )

    input:
    tuple val( id ), path( bamfile )

    output:
    tuple val( id ), path( "${id}.bw" )

    script:
    """
    samtools index "${bamfile}"
    bamCoverage \
        -b "${bamfile}" \
        --outFileName "${id}.bw" \
        --outFileFormat bigwig \
        --binSize 1 \
        --numberOfProcessors ${task.cpus}
    
    """
}


process plot_peaks {

    tag "${id}"
    label 'big_cpu'

    publishDir( 
        "${params.outputs}/coverage", 
        mode: 'copy',
        saveAs: { "${id}.${it}" },
    )

    input:
    tuple val( id ), path( bigwigs, stageAs: 'dir??/*' ), path( ann_bed )

    output:
    tuple val( id ), path( "matrix.mat.gz" ), emit: matrix
    tuple val( id ), path( "peak-heatmap.png" ), emit: plot

    script:
    """
    export MPLCONFIGDIR=mpltemp
    mkdir "\$MPLCONFIGDIR"

    computeMatrix reference-point \
        --verbose \
        -S dir??/*.bw \
        --regionsFileName "${ann_bed}" \
        --referencePoint TSS \
        -a 1000 \
        -b 1000 \
        --binSize 10 \
        --numberOfProcessors ${task.cpus} \
        -o matrix.mat.gz

    plotHeatmap \
        --verbose \
        --matrixFile matrix.mat.gz \
        --dpi 300 \
        --colorMap magma \
        --refPointLabel "Gene Start" \
        --legendLocation none \
        --outFileName peak-heatmap.png
    
    """

}