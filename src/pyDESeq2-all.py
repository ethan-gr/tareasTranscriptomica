# ==============================================================================
# IMPORTS

import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
from pydeseq2.dds import DeseqDataSet
from pydeseq2.default_inference import DefaultInference
from pydeseq2.ds import DeseqStats
import numpy as np
from numpy import log10
import os
import re

# ==============================================================================
for name in ["fc_RNAseq_bwt", "fc_RNAseq_ht2", "fc_microRNA"]:
    outdir = f"results/{name}-pyDESeq2"
    os.makedirs(outdir)

    # Load counts
    file = f"results/counts/{name}.txt"
    if name == "fc-microRNA": sampleNames = re.compile(r"((E12|E20)\.[AB])")
    else: sampleNames = re.compile(r"((5xFAD|WT)[-\.][AB])")

    counts = pd.read_csv(file, sep='\t').rename(columns=lambda x: sample.group() if (sample := re.search(sampleNames, x)) else x).T

    # Metadata
    if name == "fc-microRNA": 
        metadata = pd.DataFrame.from_dict({"Condition": {sample: "E12" if "E12" in sample else "E20" for sample in counts.index}})
    else: 
        metadata = pd.DataFrame.from_dict({"Condition": {sample: "WT" if "WT" in sample else "5xFAD" for sample in counts.index}})

    # ==============================================================================

    # Filter NaN in condition metadata
    counts = counts.loc[(~metadata.Condition.isna()).index]

    # Filter lowly expressed genes
    # Don't forget to restart when changing minExpression value
    minExpression = 3
    genesToKeep = counts.columns[counts.sum(axis=0) >= minExpression]
    counts = counts[genesToKeep]

    # ==============================================================================

    # Start DESeq
    dds = DeseqDataSet(counts=counts,  metadata=metadata, design_factors="Condition")
    dds.deseq2()

    # Statistics (DEA)
    if name == "fc-microRNA": stat_res = DeseqStats(dds, contrast = ("Condition", "E20", "E12"))
    else: stat_res = DeseqStats(dds, contrast = ("Condition", "5xFAD", "WT"))
    stat_res.summary()

    # Results
    res = stat_res.results_df

    # Remove low expressed genes
    res = res[res.baseMean >= 10]

    # Filter significants
    sigs = res[(res.padj < 0.05) & (abs(res.log2FoldChange) > 1)]
    sigs.to_csv(os.path.join(outdir, f"{name}-DEGs.tsv"), sep='\t')

    # ==============================================================================

    # Volcano plot
    dds.layers["log1p"] = np.log1p(dds.layers["normed_counts"])

    dds_sigs = dds[:, sigs.index]

    grapher = pd.DataFrame(dds_sigs.layers["log1p"].T,index=dds_sigs.var_names, columns=dds_sigs.obs_names)

    x = res.log2FoldChange
    y = -log10(res.padj.fillna(1))

    def mapColor(x, t):
        if x > t: return "UP"
        elif x < -t: return "Down"
        else: return "None"

    temp = pd.DataFrame.from_dict({"log2FoldChange":x, "-log10Padj":y})
    temp["Change"] = temp.log2FoldChange.apply(lambda x: mapColor(x, 1))

    fig,ax = plt.subplots()
    sns.scatterplot(temp, x="log2FoldChange", y="-log10Padj", hue="Change", palette=["gray", "blue", "red"])
    ax.axhline(y=1.3, linestyle='--', color="gray", alpha=0.5, linewidth=0.8)
    ax.axvline(x=1, linestyle='--', color="gray", alpha=0.5, linewidth=0.8)
    ax.axvline(x=-1, linestyle='--', color="gray", alpha=0.5, linewidth=0.8)
    xticks = [i for i in range(-10,11,5)]
    xtick_labels = [str(tick) for tick in xticks]
    plt.xticks(xticks, xtick_labels)
    plt.savefig(f"{os.path.join(outdir, name)}-volcano-plot.png")

    # ==============================================================================

    # heatmap
    sns.clustermap(grapher, z_score=0, cmap="vlag", figsize=(6,6))
    plt.savefig(f"{os.path.join(outdir, name)}-heatmap.png")
