#!/usr/bin/env python3

import argparse
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def read_table(path):
    if path.endswith((".tsv", ".txt")):
        return pd.read_csv(path, sep="\t")
    return pd.read_csv(path)


def plot_tf_landscape(
    stats_path,
    eregulon_path,
    outdir,
    prefix,
    title,
    padj_cutoff,
    min_delta,
    size_scale,
    label_top
):
    os.makedirs(outdir, exist_ok=True)

    stats = read_table(stats_path)
    ereg = read_table(eregulon_path)

    required_stats = {"TF", "padj", "delta_auc"}
    required_ereg = {"TF"}

    missing_stats = required_stats - set(stats.columns)
    missing_ereg = required_ereg - set(ereg.columns)

    if missing_stats:
        raise ValueError(f"Missing columns in {stats_path}: {missing_stats}")
    if missing_ereg:
        raise ValueError(f"Missing columns in {eregulon_path}: {missing_ereg}")

    sig = stats[
        (stats["padj"] < padj_cutoff) &
        (stats["delta_auc"] > min_delta)
    ].copy()

    print(f"{prefix}: significant positive regulons = {sig.shape[0]}")

    mean_delta = (
        sig.groupby("TF", as_index=False)["delta_auc"]
        .mean()
        .rename(columns={"delta_auc": "mean_delta_auc"})
    )

    edge_counts = (
        ereg[ereg["TF"].isin(mean_delta["TF"])]
        .groupby("TF")
        .size()
        .reset_index(name="target_edges")
    )

    tf_df = edge_counts.merge(mean_delta, on="TF", how="inner")

    if tf_df.empty:
        raise ValueError(f"{prefix}: no TFs remained after filtering and merging.")

    tf_df["normalized_activity_connectivity"] = (
        tf_df["mean_delta_auc"] *
        (tf_df["target_edges"] / tf_df["target_edges"].max())
    )

    tf_df = tf_df.sort_values(
        "normalized_activity_connectivity",
        ascending=False
    )

    source_csv = os.path.join(outdir, f"{prefix}_source_data.csv")
    tf_df.to_csv(source_csv, index=False)

    sizes = (
        np.sqrt(tf_df["normalized_activity_connectivity"].clip(lower=0)) *
        size_scale
    )
    sizes = np.clip(sizes, 120, None)

    fig, ax = plt.subplots(figsize=(6, 5))

    sc = ax.scatter(
        tf_df["target_edges"],
        tf_df["mean_delta_auc"],
        s=sizes,
        c=tf_df["normalized_activity_connectivity"],
        cmap="magma",
        edgecolor="black",
        linewidth=0.8,
        alpha=0.95
    )

    label_df = tf_df.head(label_top)

    for _, row in label_df.iterrows():
        ax.text(
            row["target_edges"] * 1.03,
            row["mean_delta_auc"] * 1.01,
            row["TF"],
            fontsize=8
        )

    ax.set_xscale("log")
    ax.set_xlabel("Δ TF target genes (HSV-1 vs Mock)")
    ax.set_ylabel("Mean Δ regulon activity (AUCell)")
    ax.set_title(title, fontsize=10)

    cbar = fig.colorbar(sc, ax=ax)
    cbar.set_label("Mean ΔAUCell")

    plt.tight_layout()

    pdf = os.path.join(outdir, f"{prefix}.pdf")
    png = os.path.join(outdir, f"{prefix}.png")

    fig.savefig(pdf, bbox_inches="tight")
    fig.savefig(png, dpi=600, bbox_inches="tight")
    plt.close(fig)

    print(f"Saved: {pdf}")
    print(f"Saved: {png}")
    print(f"Saved source data: {source_csv}")


def main():
    parser = argparse.ArgumentParser(
        description="Generate Figure 4C/D TF influence landscape panels."
    )

    parser.add_argument(
        "--stats-c",
        required=True,
        help="Stats table for Figure 4C; columns: TF, padj, delta_auc"
    )
    parser.add_argument(
        "--stats-d",
        required=True,
        help="Stats table for Figure 4D; columns: TF, padj, delta_auc"
    )
    parser.add_argument(
        "--eregulon",
        required=True,
        help="Extended eRegulon metadata table; must contain TF"
    )
    parser.add_argument(
        "--outdir",
        default="outputs/Figure4"
    )
    parser.add_argument(
        "--prefix-c",
        default="Figure4C_TF_influence_landscape"
    )
    parser.add_argument(
        "--prefix-d",
        default="Figure4D_TF_influence_landscape"
    )
    parser.add_argument(
        "--title-c",
        default="TF Influence Landscape"
    )
    parser.add_argument(
        "--title-d",
        default="TF Influence Landscape"
    )
    parser.add_argument(
        "--padj",
        type=float,
        default=1e-10
    )
    parser.add_argument(
        "--min-delta",
        type=float,
        default=0.0
    )
    parser.add_argument(
        "--size-scale",
        type=float,
        default=5000
    )
    parser.add_argument(
        "--label-top",
        type=int,
        default=25
    )

    args = parser.parse_args()

    plot_tf_landscape(
        stats_path=args.stats_c,
        eregulon_path=args.eregulon,
        outdir=args.outdir,
        prefix=args.prefix_c,
        title=args.title_c,
        padj_cutoff=args.padj,
        min_delta=args.min_delta,
        size_scale=args.size_scale,
        label_top=args.label_top
    )

    plot_tf_landscape(
        stats_path=args.stats_d,
        eregulon_path=args.eregulon,
        outdir=args.outdir,
        prefix=args.prefix_d,
        title=args.title_d,
        padj_cutoff=args.padj,
        min_delta=args.min_delta,
        size_scale=args.size_scale,
        label_top=args.label_top
    )


if __name__ == "__main__":
    main()
