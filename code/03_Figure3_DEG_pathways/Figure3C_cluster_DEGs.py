#!/usr/bin/env python3

import argparse
import os

import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, help="Input annotated h5ad")
    parser.add_argument("--outdir", default="outputs/Figure3")
    parser.add_argument("--prefix", default="Figure3C_cluster_DEGs")
    parser.add_argument("--min_cells", type=int, default=50)
    parser.add_argument("--padj", type=float, default=0.05)
    parser.add_argument("--logfc", type=float, default=0.25)
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    adata = sc.read(args.input)

    if "annotated_clusters" not in adata.obs.columns:
        raise ValueError("Missing adata.obs['annotated_clusters']")
    if "condition" not in adata.obs.columns:
        raise ValueError("Missing adata.obs['condition']")

    adata.obs["annotated_clusters"] = (
        adata.obs["annotated_clusters"]
        .astype(str)
        .str.replace("‚Äì", "-", regex=False)
        .str.strip()
    )
    adata.obs["condition"] = adata.obs["condition"].astype(str).str.strip()

    myeloid_states = [
        "Homeostatic Microglia",
        "Transiently Activated Microglia",
        "IFN-Responsive Microglia",
        "Primed Microglia",
        "Mitochondrial-Activated Microglia",
        "IEG-High Microglia",
        "Infiltrating Macrophages",
    ]

    plot_order = [
        "Infiltrating Macrophages",
        "Transiently Activated Microglia",
        "Homeostatic Microglia",
        "Mitochondrial-Activated Microglia",
        "Primed Microglia",
        "IEG-High Microglia",
    ]

    adata_myeloid = adata[
        adata.obs["annotated_clusters"].isin(myeloid_states)
        & adata.obs["condition"].isin(["Mock", "HSV1"])
    ].copy()

    cell_counts = (
        adata_myeloid.obs
        .groupby(["annotated_clusters", "condition"])
        .size()
        .unstack(fill_value=0)
    )

    cell_counts.to_csv(
        os.path.join(args.outdir, f"{args.prefix}_cell_counts_by_cluster_condition.csv")
    )

    stable_clusters = cell_counts[
        (cell_counts.get("Mock", 0) >= args.min_cells)
        & (cell_counts.get("HSV1", 0) >= args.min_cells)
    ].index.tolist()

    all_deg_tables = []
    deg_rows = []

    for clust in stable_clusters:
        sub = adata_myeloid[
            adata_myeloid.obs["annotated_clusters"] == clust
        ].copy()

        sc.tl.rank_genes_groups(
            sub,
            groupby="condition",
            groups=["HSV1"],
            reference="Mock",
            method="wilcoxon",
            key_added="deg_tmp",
            pts=False,
        )

        df = sc.get.rank_genes_groups_df(
            sub,
            group="HSV1",
            key="deg_tmp"
        )

        df["cluster"] = clust
        df["comparison"] = "HSV1_vs_Mock"

        sig = df[
            (df["pvals_adj"] < args.padj)
            & (np.abs(df["logfoldchanges"]) >= args.logfc)
        ].copy()

        deg_rows.append({
            "cluster": clust,
            "Upregulated": int((sig["logfoldchanges"] > 0).sum()),
            "Downregulated": int((sig["logfoldchanges"] < 0).sum()),
        })

        all_deg_tables.append(df)

    all_deg = pd.concat(all_deg_tables, axis=0, ignore_index=True)
    all_deg.to_csv(
        os.path.join(args.outdir, f"{args.prefix}_all_wilcoxon_results.csv"),
        index=False,
    )

    sig_deg = all_deg[
        (all_deg["pvals_adj"] < args.padj)
        & (np.abs(all_deg["logfoldchanges"]) >= args.logfc)
    ].copy()

    sig_deg.to_csv(
        os.path.join(args.outdir, f"{args.prefix}_significant_DEGs.csv"),
        index=False,
    )

    deg_df = pd.DataFrame(deg_rows)

    deg_df["Total_DEGs"] = (
        deg_df["Upregulated"] + deg_df["Downregulated"]
    )

    deg_df = (
        deg_df
        .set_index("cluster")
        .reindex(plot_order)
        .dropna(how="all")
        .fillna(0)
        .reset_index()
    )

    deg_df[["Upregulated", "Downregulated", "Total_DEGs"]] = (
        deg_df[["Upregulated", "Downregulated", "Total_DEGs"]].astype(int)
    )

    deg_df.to_csv(
        os.path.join(args.outdir, f"{args.prefix}_summary_counts.csv"),
        index=False,
    )

    fig, ax = plt.subplots(figsize=(3.0, 2.2))

    y = np.arange(len(deg_df))

    ax.barh(
        y,
        deg_df["Downregulated"],
        color="#cfcfd1",
        edgecolor="white",
        linewidth=0.4,
        label="Downregulated",
    )

    ax.barh(
        y,
        deg_df["Upregulated"],
        left=deg_df["Downregulated"],
        color="#b2182b",
        edgecolor="white",
        linewidth=0.4,
        label="Upregulated",
    )

    for i, total in enumerate(deg_df["Total_DEGs"]):
        ax.text(
            total + 3,
            i,
            str(total),
            va="center",
            ha="left",
            fontsize=8,
        )

    ax.set_yticks(y)
    ax.set_yticklabels(deg_df["cluster"], fontsize=8)
    ax.invert_yaxis()

    ax.set_xlabel("Number of DEGs (HSV-1 vs Mock)", fontsize=9)
    ax.set_xlim(0, max(190, deg_df["Total_DEGs"].max() * 1.15))

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    ax.legend(
        frameon=False,
        loc="lower right",
        fontsize=7,
        handlelength=1.4,
        labelspacing=0.25,
    )

    fig.tight_layout()

    fig.savefig(
        os.path.join(args.outdir, f"{args.prefix}.pdf"),
        dpi=300,
        bbox_inches="tight",
    )

    fig.savefig(
        os.path.join(args.outdir, f"{args.prefix}.png"),
        dpi=300,
        bbox_inches="tight",
    )


if __name__ == "__main__":
    main()
