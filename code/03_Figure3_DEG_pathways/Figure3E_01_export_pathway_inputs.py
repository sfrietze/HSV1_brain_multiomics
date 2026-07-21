#!/usr/bin/env python3

import argparse
import os

import scanpy as sc
import pandas as pd
import numpy as np


def deg_table_for_cluster(adata, cluster, max_padj, min_abs_logfc):
    sub = adata[adata.obs["annotated_clusters"] == cluster].copy()

    sub.obs["condition"] = sub.obs["condition"].astype("category")
    sub.obs["condition"] = sub.obs["condition"].cat.set_categories(
        ["Mock", "HSV1"],
        ordered=True,
    )

    sc.tl.rank_genes_groups(
        sub,
        groupby="condition",
        groups=["HSV1"],
        reference="Mock",
        method="wilcoxon",
        key_added="deg_tmp",
        pts=False,
    )

    df = sc.get.rank_genes_groups_df(sub, group="HSV1", key="deg_tmp")
    df = df.rename(
        columns={
            "names": "gene",
            "logfoldchanges": "logFC",
            "pvals_adj": "padj",
        }
    )

    df = df[["gene", "logFC", "padj", "scores"]].copy()
    sig = df[
        (df["padj"] <= max_padj)
        & (df["logFC"].abs() >= min_abs_logfc)
    ].copy()

    up = sig[sig["logFC"] > 0].copy()
    down = sig[sig["logFC"] < 0].copy()

    return df, sig, up, down


def ensure_markers(adata, key="markers"):
    if key not in adata.uns:
        sc.tl.rank_genes_groups(
            adata,
            groupby="annotated_clusters",
            method="wilcoxon",
            key_added=key,
            pts=False,
        )


def marker_genes_for_cluster(
    adata,
    cluster,
    n=300,
    min_logfc=0.25,
    max_padj=0.05,
):
    ensure_markers(adata, key="markers")

    df = sc.get.rank_genes_groups_df(
        adata,
        group=cluster,
        key="markers",
    )

    df = df.rename(
        columns={
            "names": "gene",
            "logfoldchanges": "logFC",
            "pvals_adj": "padj",
        }
    )

    keep = df[
        (df["padj"] <= max_padj)
        & (df["logFC"] >= min_logfc)
    ].copy()

    return keep.sort_values("padj").head(n)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True)
    parser.add_argument("--outdir", default="data/03_Figure3_DEG_pathways/pathway_inputs")
    parser.add_argument("--max_padj", type=float, default=0.05)
    parser.add_argument("--min_abs_logfc", type=float, default=0.25)
    parser.add_argument("--ifn_marker_n", type=int, default=300)
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    adata = sc.read(args.input)

    adata.obs["annotated_clusters"] = (
        adata.obs["annotated_clusters"]
        .astype(str)
        .str.replace("‚Äì", "-", regex=False)
        .str.strip()
    )
    adata.obs["condition"] = adata.obs["condition"].astype(str).str.strip()

    stable_clusters = [
        "Infiltrating Macrophages",
        "Transiently Activated Microglia",
        "Homeostatic Microglia",
    ]

    ifn_cluster = "IFN-Responsive Microglia"

    all_metascape_lists = {}
    summary_rows = []

    for cluster in stable_clusters:
        full_df, sig, up, down = deg_table_for_cluster(
            adata,
            cluster,
            max_padj=args.max_padj,
            min_abs_logfc=args.min_abs_logfc,
        )

        tag = cluster.replace(" ", "_")

        full_df.to_csv(
            os.path.join(args.outdir, f"{tag}_DE_full.csv"),
            index=False,
        )
        sig.to_csv(
            os.path.join(args.outdir, f"{tag}_DE_sig_ALL.csv"),
            index=False,
        )
        up.to_csv(
            os.path.join(args.outdir, f"{tag}_DE_sig_UP.csv"),
            index=False,
        )
        down.to_csv(
            os.path.join(args.outdir, f"{tag}_DE_sig_DOWN.csv"),
            index=False,
        )

        up_genes = up["gene"].dropna().astype(str).unique().tolist()
        down_genes = down["gene"].dropna().astype(str).unique().tolist()
        all_genes = sig["gene"].dropna().astype(str).unique().tolist()

        pd.Series(up_genes).to_csv(
            os.path.join(args.outdir, f"{tag}_UP_genes.txt"),
            index=False,
            header=False,
        )
        pd.Series(down_genes).to_csv(
            os.path.join(args.outdir, f"{tag}_DOWN_genes.txt"),
            index=False,
            header=False,
        )
        pd.Series(all_genes).to_csv(
            os.path.join(args.outdir, f"{tag}_ALL_genes.txt"),
            index=False,
            header=False,
        )

        all_metascape_lists[f"{cluster} (HSV1_UP)"] = up_genes
        all_metascape_lists[f"{cluster} (HSV1_DOWN)"] = down_genes
        all_metascape_lists[f"{cluster} (ALL_DEGs)"] = all_genes

        summary_rows.append(
            {
                "cluster": cluster,
                "n_all_deg": len(all_genes),
                "n_up": len(up_genes),
                "n_down": len(down_genes),
            }
        )

    if ifn_cluster in adata.obs["annotated_clusters"].unique():
        ifn_markers = marker_genes_for_cluster(
            adata,
            ifn_cluster,
            n=args.ifn_marker_n,
            min_logfc=args.min_abs_logfc,
            max_padj=args.max_padj,
        )

        tag = ifn_cluster.replace(" ", "_")

        ifn_markers.to_csv(
            os.path.join(args.outdir, f"{tag}_MARKERS.csv"),
            index=False,
        )

        ifn_genes = ifn_markers["gene"].dropna().astype(str).unique().tolist()

        pd.Series(ifn_genes).to_csv(
            os.path.join(args.outdir, f"{tag}_MARKERS_genes.txt"),
            index=False,
            header=False,
        )

        all_metascape_lists[f"{ifn_cluster} (MARKERS_pos)"] = ifn_genes

        summary_rows.append(
            {
                "cluster": ifn_cluster,
                "n_all_deg": len(ifn_genes),
                "n_up": len(ifn_genes),
                "n_down": 0,
            }
        )

    summary_df = pd.DataFrame(summary_rows)
    summary_df.to_csv(
        os.path.join(args.outdir, "DEG_counts_summary.csv"),
        index=False,
    )

    max_len = max((len(v) for v in all_metascape_lists.values()), default=0)
    metascape_df = pd.DataFrame(
        {
            k: v + [""] * (max_len - len(v))
            for k, v in all_metascape_lists.items()
        }
    )

    metascape_df.to_csv(
        os.path.join(args.outdir, "Metascape_multilist.csv"),
        index=False,
    )


if __name__ == "__main__":
    main()
