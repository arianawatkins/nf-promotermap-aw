process macs3 {

    tag "${id}"
    label 'med_mem'

    publishDir( 
        "${params.outputs}/peaks", 
        mode: 'copy',
        saveAs: { "${id}.${it}" },
    )

    input:
    tuple val( id ), path( bamfile ), path( bamfile_ctrl, stageAs: "control/*" )

    output:
    tuple val( id ), path( "peaks.bed" ), emit: peaks
    tuple val( id ), path( "summits.bed" ), emit: summits
    tuple val( id ), path( "macs3/*_peaks.tsv" ), emit: excel
    tuple val( id ), path( "macs3/*.bed" ), emit: bed
    tuple val( id ), path( "macs3/*.bdg" ), emit: bg
    tuple val( id ), path( "macs3/*_cutoff_analysis.txt" ), emit: cutoffs
    tuple val( id ), path( "peak_length_histogram.png" ), emit: histogram 

    script:
    """
    samtools merge -@ ${task.cpus} -o merged.bam ${bamfile}
    samtools view -bS merged.bam -h -F 20 -o f.bam
    samtools view -bS merged.bam -h -f 16 -o r.bam

    samtools view -bS "${bamfile_ctrl}" -h -F 20 -o f_ctrl.bam
    samtools view -bS "${bamfile_ctrl}" -h -f 16 -o r_ctrl.bam

    samtools view -H "${bamfile_ctrl}" \
    | awk '\
        \$1 == "@SQ" { 
            for(i=1; i<=NF; i++) { 
                if(\$i ~ /^LN:/){
                    split(\$i, a, ":"); 
                    len+=a[2]
                }
            } 
        }
        END { print len }
    ' \
    > genome-size.txt

    for f in {f,r}.bam
    do
        strand=\$(basename \$f .bam)
        ctrl=\$(basename \$f .bam)_ctrl.bam
        macs3 callpeak \
            --treatment "\$f" \
            --control "\$ctrl" \
            --format BAMPE \
            --call-summits \
            -p 0.01 \
            --min-length 30 \
            --cutoff-analysis \
            --outdir macs3 \
            --name "${id}.\$strand" \
            --bdg \
            --trackline \
            --gsize \$(cat genome-size.txt) \
            --buffer-size ${Math.round(task.memory.toMega() * 0.8 / 800)}
    done

    for f in macs3/*.xls
    do
        # these aren't actually xls files
        mv \$f \$(dirname \$f)/\$(basename \$f .xls).tsv
    done

    summit_BED=(macs3/*.bed)
    head -n1 "\${summit_BED[0]}" \
    | sed 's/${id}\\.f/${id}/g' \
    | cat - \
        <(
            cat \
                <(awk -v OFS='\\t' 'NR>1 { print \$0, "+" }' "\${summit_BED[0]}") \
                <(awk -v OFS='\\t' 'NR>1 { print \$0, "-" }' "\${summit_BED[1]}") \
            | sort -k1,1 -k2,2n
        ) \
    > "summits.bed"

    peak_tsv=(macs3/*_peaks.tsv)
    cat \
        <(echo 'track name="${id} (peaks)" description="Peaks for ${id} (Made with MACS v3, 09/03/25)" visibility=1') \
        <(
            cat "\${peak_tsv[@]}" \
            | grep -v -e '^#' -e '^\$' -e '^chr\\s' \
            | awk -v OFS='\\t' '
                {
                    chr=\$1; start=\$2; end=\$3; name=\$10;
                    # infer strand from name tag
                    strand=".";
                    if (name ~ /\\.f_/) strand="+";
                    else if (name ~ /\\.r_/) strand="-";
                    # convert start to 0-based; end stays as-is (end-exclusive in BED)
                    print chr, start-1, end, name, 0, strand
                }
            ' \
            | sort -k1,1 -k2,2n \
        ) \
    > "peaks.bed"

 # Calculate peak lengths and generate histogram of number of peaks vs peak length
    awk 'NR>1 {print \$3-\$2}' peaks.bed > peak_lengths.txt

    python3 - <<EOF
import matplotlib.pyplot as plt
import numpy as np

# Load peak lengths
lengths = np.loadtxt('peak_lengths.txt')

# Create adaptive bins for axis
bins = np.arange(min(lengths), max(lengths) + 2) - 0.5

# Calculate mean and median
mean_len = round(np.mean(lengths))
median_len = round(np.median(lengths))

# Print stats
print(f"Mean peak length: {mean_len} bp")
print(f"Median peak length: {median_len} bp")

# Plot histogram of number of peaks vs peak length
plt.figure(figsize=(8,6))
plt.hist(lengths, bins=bins, color='cornflowerblue')
plt.xlabel('Peak Length (bp)', fontsize=12)
plt.ylabel('Number of Peaks', fontsize=12)
plt.title('Peak Length Distribution (${id})', fontsize=14)
plt.grid(alpha=0.3)

# Add median and mean lines
plt.axvline(median_len, color='midnightblue', linestyle='--', label=f'Median = {median_len} bp')
plt.axvline(mean_len, color='midnightblue', linestyle='-', label=f'Mean = {mean_len} bp')
plt.legend()

# Save figure
plt.tight_layout()
plt.savefig("peak_length_histogram.png", dpi=300)
EOF

    rm merged.bam

    """
}
