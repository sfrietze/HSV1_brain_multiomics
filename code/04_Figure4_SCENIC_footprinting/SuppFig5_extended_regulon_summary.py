#!/usr/bin/env python3

import argparse
import os
import numpy as np
import pandas as pd
import mudata as md
import matplotlib.pyplot as plt
import seaborn as sns


def clean_celltype(x):
    x = str(x)
    x = x.replace("_Mock", "").replace("_HSV1", "")
    return x


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--scplus", required=True, help="SCENIC+ scplusmdata.h5mu")
    parser.add_argument("--direct", required=True, help="eRegulon_direct.tsv")
    parser.add_argument("--extended", required=True, help="eRegulons_extended.tsv")
    parser.add_argument("--outdir", default="outputs/SuppFig5")
    parser.add_argument("--top_n", type=int, default=40)
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    group_var = "scRNA_counts:celltype_condition"

    # -------------------------------
    # Panel A: extended regulon activity heatmap
    # -------------------------------
    mdata = md.read_h5mu(args.scplus)

    auc = mdata["extended_gene_based_AUC"].to_df()
    obs = mdata.obs.copy()

    obs["celltype"] = obs[group_var].map(clean_celltype)
    auc["celltype"] = obs["celltype"].values

    avg_auc = auc.groupby("celltype").mean(numeric_only=True)

    # keep extended regulons only
    reg_cols = [c for c in avg_auc.columns if "extended" in c]
    avg_auc = avg_auc[reg_cols]

    # z-score each regulon across celltypes
    z = avg_auc.copy()
    z = (z - z.mean(axis=0)) / z.std(axis=0)
    z = z.replace([np.inf, -np.inf], np.nan).dropna(axis=1)

    # choose most variable regulons for readable supplemental panel
    top_cols = z.var(axis=0).sort_values(ascending=False).head(args.top_n).index
    heatmap_data = z[top_cols]

    heatmap_data.to_csv(os.path.join(args.outdir, "SuppFig5A_extended_regulon_activity_zscores.csv"))

    g = sns.clustermap(
        heatmap_data,
        cmap="RdBu_r",
        center=0,
        metric="correlation",
        method="average",
        linewidths=0.2,
        figsize=(8, 8),
        cbar_kws={"label": "Z-scored AUCell"},
    )

    g.savefig(os.path.join(args.outdir, "SuppFig5A_extended_regulon_activity_heatmap.pdf"),
              dpi=300, bbox_inches="tight")
    g.savefig(os.path.join(args.outdir, "SuppFig5A_extended_regulon_activity_heatmap.png"),
              dpi=300, bbox_inches="tight")
    plt.close("all")

    # -------------------------------
    # Panel B: regulon size distribution
    # -------------------------------
    direct = pd.read_csv(args.direct, sep="\t")
    extended = pd.read_csv(args.extended, sep="\t")

    direct_sizes = (
        direct.groupby("TF")["Gene"]
        .nunique()
        .reset_index(name="n_targets")
    )
    direct_sizes["type"] = "Direct"

    extended_sizes = (
        extended.groupby("TF")["Gene"]
        .nunique()
        .reset_index(name="n_targets")
    )
    extended_sizes["type"] = "Extended"

    sizes = pd.concat([direct_sizes, extended_sizes], ignore_index=True)
    sizes.to_csv(os.path.join(args.outdir, "SuppFig5B_regulon_size_distribution_source.csv"),
                 index=False)

    plt.figure(figsize=(6, 5))
    for label, df in sizes.groupby("type"):
        plt.hist(df["n_targets"], bins=20, alpha=0.6, label=label)

    plt.xlabel("Number of Target Genes per Regulon")
    plt.ylabel("Frequency")
    plt.legend()
    plt.tight_layout()

    plt.savefig(os.path.join(args.outdir, "SuppFig5B_regulon_size_distribution.pdf"),
                dpi=300, bbox_inches="tight")
    plt.savefig(os.path.join(args.outdir, "SuppFig5B_regulon_size_distribution.png"),
                dpi=300, bbox_inches="tight")
    plt.close("all")

    print("Saved outputs to:", args.outdir)
    print("Direct median targets:", direct_sizes["n_targets"].median())
    print("Extended median targets:", extended_sizes["n_targets"].median())


if __name__ == "__main__":
    main()
