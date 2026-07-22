#!/usr/bin/env python3

import argparse
import os
import pandas as pd
import matplotlib.pyplot as plt


def read_gene_list(path):
    if path.endswith(".csv"):
        df = pd.read_csv(path)
        if {"gene", "padj", "logFC"}.issubset(df.columns):
            return df[(df["padj"] < 0.05) & (df["logFC"] > 0.25)]["gene"].astype(str).tolist()
        return df.iloc[:, 0].astype(str).tolist()
    return pd.read_csv(path, header=None)[0].astype(str).tolist()


def make_panel(stats, ereg, pos_genes, ref_genes, pos_label, ref_label, padj):
    sig = stats[(stats["padj"] < padj) & (stats["delta_auc"] > 0)].copy()

    mean_delta = (
        sig.groupby("TF")["delta_auc"]
        .mean()
        .reset_index()
        .rename(columns={"delta_auc": "mean_delta"})
    )

    rows = []
    for tf in mean_delta["TF"]:
        targets = set(ereg.loc[ereg["TF"] == tf, "Gene"].astype(str))
        e_pos = len(targets & set(pos_genes))
        e_ref = len(targets & set(ref_genes))
        rows.append({
            "TF": tf,
            f"edges_{pos_label}": e_pos,
            f"edges_{ref_label}": e_ref,
            "delta_edges": e_pos - e_ref
        })

    edge_df = pd.DataFrame(rows)

    tf = mean_delta.merge(edge_df, on="TF")
    tf = tf[tf["delta_edges"] > 0].copy()
    tf = tf.sort_values(["delta_edges", "mean_delta"], ascending=False)

    return tf


def plot_panel(ax, tf, title, xlabel, highlight):
    size_scale = 1200
    sizes = 150 + size_scale * (tf["delta_edges"] / tf["delta_edges"].max())

    sc = ax.scatter(
        tf["delta_edges"],
        tf["mean_delta"],
        s=sizes,
        c=tf["mean_delta"],
        cmap="magma",
        vmin=0,
        vmax=0.12,
        edgecolor="black",
        linewidth=0.6,
        alpha=0.9
    )

    for _, row in tf.iterrows():
        if row["TF"] in highlight:
            ax.annotate(
                row["TF"],
                (row["delta_edges"], row["mean_delta"]),
                xytext=(6, 6),
                textcoords="offset points",
                fontsize=9,
                weight="bold"
            )

    ax.axvline(0, linestyle="--", color="gray", linewidth=1)
    ax.axhline(0, linestyle="--", color="gray", linewidth=1)
    ax.set_xlabel(xlabel)
    ax.set_ylabel("Mean Regulon Activity Shift (ΔAUCell)")
    ax.set_title(title)
    ax.set_ylim(0, 0.13)
    ax.set_xlim(-5, tf["delta_edges"].max() * 1.15)

    return sc


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--eregulon", required=True)
    parser.add_argument("--stats-ifn", required=True)
    parser.add_argument("--stats-mac", required=True)
    parser.add_argument("--ifn-up", required=True)
    parser.add_argument("--homeo-up", required=True)
    parser.add_argument("--mac-up", required=True)
    parser.add_argument("--mac-down", required=True)
    parser.add_argument("--outdir", default="outputs/Figure4")
    parser.add_argument("--padj", type=float, default=1e-10)
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    ereg = pd.read_csv(args.eregulon)
    ereg.columns = ereg.columns.str.strip()
    stats_ifn = pd.read_csv(args.stats_ifn)
    stats_mac = pd.read_csv(args.stats_mac)

    ifn_genes = read_gene_list(args.ifn_up)
    homeo_genes = read_gene_list(args.homeo_up)
    mac_up = read_gene_list(args.mac_up)
    mac_down = read_gene_list(args.mac_down)

    tf_ifn = make_panel(
        stats_ifn, ereg,
        ifn_genes, homeo_genes,
        "IFN", "Homeostatic",
        args.padj
    )

    tf_mac = make_panel(
        stats_mac, ereg,
        mac_up, mac_down,
        "HSV1", "Mock",
        args.padj
    )

    tf_ifn.to_csv(os.path.join(args.outdir, "Figure4C_IFN_Rewiring_source_data.csv"), index=False)
    tf_mac.to_csv(os.path.join(args.outdir, "Figure4D_Macrophage_Rewiring_source_data.csv"), index=False)

    fig, axes = plt.subplots(1, 2, figsize=(12, 5.5), constrained_layout=True)

    sc1 = plot_panel(
        axes[0],
        tf_ifn,
        "Figure 4C: IFN-responsive Microglia",
        "Δ Supported Target Edges (IFN − Homeostatic)",
        ["Stat1", "Stat2", "Irf1", "Etv6", "Jun", "Junb", "Jund", "Atf3"]
    )

    sc2 = plot_panel(
        axes[1],
        tf_mac,
        "Figure 4D: Infiltrating Macrophages",
        "Δ Supported Target Edges (HSV-1 − Mock)",
        ["Stat1", "Stat2", "Irf1", "Stat3", "Jun", "Junb", "Jund", "Atf3", "Etv6"]
    )

    cbar = fig.colorbar(sc2, ax=axes, fraction=0.03)
    cbar.set_label("Mean ΔAUCell")

    pdf = os.path.join(args.outdir, "Figure4CD_Regulon_Rewiring.pdf")
    png = os.path.join(args.outdir, "Figure4CD_Regulon_Rewiring.png")

    fig.savefig(pdf, dpi=300, bbox_inches="tight")
    fig.savefig(png, dpi=600, bbox_inches="tight")
    plt.close(fig)

    print(f"Saved: {pdf}")
    print(f"Saved: {png}")


if __name__ == "__main__":
    main()
