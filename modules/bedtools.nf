process gff2bed {

    tag "${id}"

    publishDir( 
        "${params.outputs}/genome", 
        mode: 'copy',
        saveAs: { "${id}.${it}" },
    )

    input:
    tuple val( id ), path( gff )

    output:
    tuple val( id ), path( "genes.bed" )

    script:
    """
    awk '\$3 == "gene"' "${gff}" \
    | gff2bed --do-not-sort \
    | awk -v OFS='\\t' '
        { 
            delete c
            split(\$10, a, ";"); 
            for (i in a) {
                split(a[i], b, "="); 
                c[b[1]] = b[2]
            }; 
            if ("old_locus_tag" in c) { 
                key = "old_locus_tag" 
            } else if ("locus_tag" in c) {
                key = "locus_tag"
            } else {
                key = "ID"
            }; 
            \$4 = (key == "" ? "." : c[key])
            print \$0;
        }
    ' \
    | sort -k1,1 -k2,2n \
    > genes.bed

    """
}

process annotate_nearest_gene {

    tag "${id}"

    publishDir( 
        "${params.outputs}/peaks", 
        mode: 'copy',
        saveAs: { "${id}.${it}" },
    )

    input:
    tuple val( id ), path( summits ), path( ann_bed )

    output:
    tuple val( id ), path( "summits-ann.tsv" )

    script:
    """
    grep -v '^track' "${summits}" \
    | sort -k1,1 -k2,2n \
    > summits.sorted.bed

    printf 'chr\\tsummit_start\\tsummit_end\\tpeak_name\\tscore\\tstrand\\tann_chr\\tann_start\\tann_end\\tann_name\\tann_score\\tann_strand\\tann_source\\tann_type\\tann_phase\\tann_attributes\\tann_distance\\n' \
    > summits-ann0.tsv
    bedtools closest -a summits.sorted.bed -b "${ann_bed}" \
        -s -D b -t first -iu \
    >> summits-ann0.tsv

    python -c '
    import pandas as pd
    from collections import defaultdict

    df = pd.read_csv("summits-ann0.tsv", sep="\\t")
    attrs = {
        "ann_" + a.split("=")[0]: [] 
        for v in df["ann_attributes"] 
        for a in v.split(";")
    }

    for v in df["ann_attributes"]:
        these_attrs = {a.split("=")[0]: a.split("=")[1] for a in v.split(";") if "=" in a}
        for key in attrs:
            attrs[key].append(these_attrs.get(key.replace("ann_", "")))

    (
        df
        .assign(**attrs)
        .drop(columns="ann_attributes")
        .to_csv("summits-ann.tsv", sep="\\t", index=False)
    )
    
    '

    """
}


process extract_peak_sequences {

    tag "${id}"

    publishDir( 
        "${params.outputs}/peaks", 
        mode: 'copy',
        saveAs: { "${id}.${it}" },
    )


    // Need to add the input val( orientation ) here - need to figure out what is going in 
    input:
    tuple val( id ), path( peaks ), path( genome ) //, val( orientation )

    output:
    tuple val( id ), path( "peak-seqs.fasta" ), emit: fasta
    tuple val( id ), path( "peak-seqs.tsv" ), emit: tsv

    script:
    """
    bedtools getfasta -s -fi "${genome}" -bed "${peaks}" -name > peak-seqs.fasta
    printf 'peak_id\\tpeak_name\\tchr\\tpeak_start\\tpeak_end\\tpeak_sequence\\n' \
    | cat - <(
        bedtools getfasta -s -fi "${genome}" -bed "${peaks}" -name -tab \
        | awk -F'\\t' -v OFS='\\t' '
            { 
                split(\$1, a, ":"); 
                split(a[4], c0, "-"); 
                split(c0[2], c1, "("); 
                print \$1, a[1], a[3], c0[1], c1[1], \$2 
            }
        ' \
    ) \
    > peak-seqs.tsv

    """
}

process get_peak_coverage {

    tag "${id}"

    publishDir( 
        "${params.outputs}/peaks", 
        mode: 'copy',
        saveAs: { "${id}.${it}" },
    )

    input:
    tuple val( id ), path( peaks ), path( bams )

    output:
    tuple val( id ), path( "peak-cov.bed" ), emit: bed
    tuple val( id ), path( "peak-cov.tsv" ), emit: tsv

    script:
    """
    for f in ${bams}
    do 
        samtools index -@ ${task.cpus} "\$f"
    done

    bedtools multicov -p -bams ${bams} -bed ${peaks} > peak-cov.bed
    printf 'chr\\tpeak_start\\tpeak_end\\tpeak_name\\tpeak_score\\tpeak_strand\\t${bams.join('\\t')}\\n' > peak-cov.tsv
    cat peak-cov.bed >> peak-cov.tsv

    """
}


process get_per_base_coverage_within_peaks {

    tag "${id}"

    publishDir( 
        "${params.outputs}/peaks", 
        mode: 'copy',
        saveAs: { "${id}.${it}" },
    )

    input:
    tuple val( id ), path( peaks ), path( bams )

    output:
    tuple val( id ), path( "coverage.tsv" ), emit: counts
    tuple val( id ), path( "coverage-norm.tsv" ), emit: normalized

    script:
    """
    for f in ${bams}
    do
        bedfile=\$(basename "\$f" .bam).fragments.bed
        covfile=\$(basename "\$f" .bam).coverage.tsv
        samtools sort -N -@ ${task.cpus} -m ${Math.round(task.memory.getGiga() * 0.8)}G "\$f" \
        | samtools view -bf 0x2 - \
        | bedtools bamtobed -bedpe -i stdin \
        | awk -v OFS='\\t' '
            { print \$1, \$2, \$6, \$7, \$8, \$9; next }
        ' \
        | sort -k1,1 -k2,2n \
        > "\$bedfile"
        printf 'chr\\tpeak_start\\tpeak_end\\tpeak_name\\tpeak_score\\tpeak_strand\\tpeak_coordinate\\t'"\$(basename "\$f" .sorted.bam)"'\\n' > "\$covfile"
        bedtools coverage -sorted -s -d -a ${peaks} -b "\$bedfile" >> "\$covfile"
    done

    #covfiles=(*.coverage.tsv)
    #head -n1 "\${covfiles[0]}" \
    #| cat - <(tail -n+2 -q "\${covfiles[@]}") \
    #> "coverage.tsv"
    #exit 1

    python -c '
    from glob import glob
    
    import numpy as np
    import pandas as pd
    from scipy.special import xlogy
    from scipy.stats import skew, kurtosis

    def norm_df(df):
        return df.values / df.values.sum(axis=1, keepdims=True)
        

    def entropy(df):
        normed = df.values / df.values.sum(axis=1, keepdims=True)
        return -np.nansum(xlogy(normed, normed), axis=1)

    
    def sparsity_index(df):
        normed = df.values / df.values.sum(axis=1, keepdims=True)
        n = df.shape[1]
        sqrt_n = np.sqrt(n)
        return (
            sqrt_n - np.linalg.norm(df.values, 1, axis=1, keepdims=True) 
            / np.linalg.norm(df.values, 2, axis=1, keepdims=True)
        ) / (sqrt_n - 1.)


    def _var(df):
        v = np.arange(df.shape[1], dtype=np.float64)[None]
        normed = df.values / df.values.sum(axis=1, keepdims=True)
        means = (normed * v).sum(axis=1, keepdims=True)
        means2 = (normed * np.square(v)).sum(axis=1, keepdims=True)
        return (means2 - np.square(means)).flatten()

    
    def _bc(df):
        normed = norm_df(df)
        g = skew(normed, axis=1, keepdims=True)
        k = kurtosis(normed, fisher=False, axis=1, keepdims=True) 
        n = normed.shape[1]
        return (np.square(g) + 1.) / (k + 3. * np.square(n - 1.)) / ((n - 2.) * (n - 3.))

    files = sorted(glob("*.coverage.tsv"))
    df = pd.read_csv(files[0], sep="\\t")
    for f in files[1:]:
        df = df.merge(pd.read_csv(f, sep="\\t"), how="outer")
    df.fillna(0.).to_csv("coverage.tsv", sep="\\t", index=False)

    df = df.set_index(df.columns[:7].tolist())
    peak_size_kb = np.abs(df.index.get_level_values("peak_start").values - df.index.get_level_values("peak_end").values) / 1000.
    df_rpkm = df.values / (peak_size_kb[:,None] * df.sum(axis=0).values[None,:] / 1_000_000.)
    df_rpkm = pd.DataFrame(df_rpkm, index=df.index, columns=df.columns)

    df2 = (
        df_rpkm
        .assign(
            base_sum=df_rpkm.sum(axis=1),
            base_var=_var(df_rpkm),
            base_entropy=entropy(df_rpkm),
            base_sparsity=sparsity_index(df_rpkm),
            base_bc=_bc(df_rpkm),
            base_max=df_rpkm.idxmax(axis="columns"),
            base_max_n=lambda x: x["base_max"].map(sorted(df_rpkm.columns.tolist()).index),
        )
    )
    df2.to_csv("coverage-norm.tsv", sep="\\t")
    
    '

    """
}


process make_peak_summary_table {

    tag "${id}"

    publishDir( 
        "${params.outputs}/peaks", 
        mode: 'copy',
        saveAs: { "${id}.${it}" },
    )

    input:
    tuple val( id ), path( summits ), path( peaks ), path( coverage )

    output:
    tuple val( id ), path( "peak-summary.tsv" )

    script:
    """
    #!/usr/bin/env python

    import pandas as pd

    (
        pd.read_csv("${summits}", sep="\\t")
        .merge(
            pd.read_csv("${coverage}", sep="\\t"),
            how="left",
        )
        .merge(
            pd.read_csv("${peaks}", sep="\\t"),
            how="left",
        )
        .to_csv("peak-summary.tsv", sep="\\t", index=False)
    )
    """
}