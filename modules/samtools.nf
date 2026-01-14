process SAMtools_stats {

    tag "${id}"
    label 'big_cpu'

    publishDir( 
        "${params.outputs}/samtools", 
        mode: 'copy',
    )

    input:
    tuple val( id ), path( bamfile )

    output:
    tuple val( id ), path( "${id}.stats.txt" )

    script:
    """
    samtools stats -@ ${task.cpus} "${bamfile}" > "${id}.stats.txt"
    """

}

process plot_bamstats {

    tag "${id}"

    publishDir( 
        "${params.outputs}/samtools", 
        mode: 'copy',
        saveAs: { "${id}.${it}" },
    )

    input:
    tuple val( id ), path( txt )

    output:
    tuple val( id ), path( "stats-*.png" )

    script:
    """
    plot-bamstats -p stats "${txt}"
    """

}



process SAMtools_flagstat {

    tag "${id}"
    label 'big_cpu'

    publishDir( 
        "${params.outputs}/samtools", 
        mode: 'copy',
    )

    input:
    tuple val( id ), path( bamfile )

    output:
    tuple val( id ), path( "${id}.flagstat.tsv" )

    script:
    """
    #samtools collate -@ ${task.cpus} -O -u "${bamfile}" \
    #| samtools fixmate -@ ${task.cpus} -m -u - - \
    #| samtools sort -@ ${task.cpus} -m ${Math.round(task.memory.getGiga() * 0.8)}G -u - \
    #| samtools markdup -@ ${task.cpus} - markdup.bam
    #samtools flagstat -@ ${task.cpus} -O tsv markdup.bam > flagstat.tsv
    samtools flagstat -@ ${task.cpus} -O tsv "${bamfile}" > "${id}.flagstat.tsv"
    """

}

process SAMtools_coverage {

    tag "${id}"

    publishDir( 
        "${params.outputs}/samtools", 
        mode: 'copy',
    )

    input:
    tuple val( id ), path( bamfile )

    output:
    tuple val( id ), path( "${id}.coverage.tsv" )

    script:
    """
    samtools coverage "${bamfile}" -o "${id}.coverage.tsv"
    """

}

process sort_and_index_bam {

    tag "${id}"
    label 'big_cpu'

    publishDir( 
        "${params.outputs}/mapped", 
        mode: 'copy',
    )

    input:
    tuple val( id ), path( bamfile )

    output:
    tuple val( id ), path( "${id}.sorted.bam" )

    script:
    """
    # filter out multimappers
    samtools view -@ ${task.cpus} -F260 -b -q 25 "${bamfile}" -o "${id}.filtered.bam"
    samtools sort -@ ${task.cpus} -m ${Math.round(task.memory.getGiga() * 0.8)}G "${id}.filtered.bam" -o "${id}.sorted.bam"
    #samtools index -@ ${task.cpus} "${id}.sorted.bam" -o "${id}.sorted.bai"
    #rm filtered.bam
    """

}


