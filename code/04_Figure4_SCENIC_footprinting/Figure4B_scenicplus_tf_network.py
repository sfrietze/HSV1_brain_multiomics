#!/usr/bin/env python3

import argparse
import os
import re

import mudata as md
import numpy as np
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt


def clean_state(x):
    x = str(x)
    x = x.replace("_Mock", "").replace("_HSV1", "")
    x = x.replace("_", " ")
    return x


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, help="SCENIC+ scplusmdata.h5mu")
    parser.add_argument("--outdir", default="outputs/Figure4")
    parser.add_argument("--prefix", default="Figure4B_SCENICplus_TF_network")
    parser.add_argument("--min_shared_genes", type=int, default=8)
    parser.add_argument("--top_tfs", type=int, default=35)
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    group_var = "scRNA_counts:celltype_condition"

    states_of_interest = [
        "IFN-Responsive Microglia",
        "Infiltrating Macrophages",
        "Primed Microglia",
        "Transiently Activated Microglia",
    ]

    state_colors = {
        "IFN-Responsive Microglia": "#5e3c99",
        "Infiltrating Macrophages": "#8b1a1a",
        "Primed Microglia": "#ef3b2c",
        "Transiently Activated Microglia": "#fcbba1",
    }

    mdata = md.read_h5mu(args.input)

    meta = pd.DataFrame(mdata.uns["direct_e_regulon_metadata"]).copy()

    required = ["TF", "Gene", "Gene_signature_name"]
    for c in required:
        if c not in meta.columns:
            raise ValueError(f"Missing column in direct_e_regulon_metadata: {c}")

    auc = mdata["direct_gene_based_AUC"].to_df()
    obs = mdata.obs.copy()

    obs["state_clean"] = obs[group_var].map(clean_state)

    # Average AUC per state
    auc["state_clean"] = obs["state_clean"].values
    auc_state = auc.groupby("state_clean").mean(numeric_only=True)

    available_states = [s for s in states_of_interest if s in auc_state.index]
    if len(available_states) == 0:
        raise ValueError("None of the requested states found after cleaning group labels.")

    # TF to target genes
    tf_targets = (
        meta.groupby("TF")["Gene"]
        .apply(lambda x: set(x.dropna().astype(str)))
        .to_dict()
    )

    # TF to SCENIC signature columns
    tf_sigs = (
        meta.groupby("TF")["Gene_signature_name"]
        .apply(lambda x: sorted(set(x.dropna().astype(str))))
        .to_dict()
    )

    # Dominant state and max AUC per TF
    tf_rows = []
    for tf, sigs in tf_sigs.items():
        sigs_present = [s for s in sigs if s in auc_state.columns]
        if len(sigs_present) == 0:
            continue

        vals = auc_state.loc[available_states, sigs_present].mean(axis=1)
        dominant = vals.idxmax()
        max_auc = float(vals.max())

        tf_rows.append({
            "TF": tf,
            "dominant_state": dominant,
            "max_auc": max_auc,
            "n_targets": len(tf_targets.get(tf, set())),
        })

    tf_df = pd.DataFrame(tf_rows)

    # Keep TFs dominated by the requested states
    tf_df = tf_df[tf_df["dominant_state"].isin(states_of_interest)].copy()

    # Prioritize TFs with high AUC and target count
    tf_df["rank_score"] = tf_df["max_auc"] * np.log2(tf_df["n_targets"] + 1)
    tf_df = tf_df.sort_values("rank_score", ascending=False).head(args.top_tfs)

    keep_tfs = tf_df["TF"].tolist()

    # Build network based on shared target genes
    G = nx.Graph()

    for _, row in tf_df.iterrows():
        tf = row["TF"]
        G.add_node(
            tf,
            dominant_state=row["dominant_state"],
            max_auc=row["max_auc"],
            n_targets=row["n_targets"],
        )

    for i, tf1 in enumerate(keep_tfs):
        for tf2 in keep_tfs[i + 1:]:
            g1 = tf_targets.get(tf1, set())
            g2 = tf_targets.get(tf2, set())
            shared = len(g1 & g2)
            if shared >= args.min_shared_genes:
                union = len(g1 | g2)
                jaccard = shared / union if union > 0 else 0
                G.add_edge(tf1, tf2, shared_genes=shared, jaccard=jaccard)

    # Drop isolated nodes only if many exist
    isolates = list(nx.isolates(G))
    if len(isolates) > 5:
        G.remove_nodes_from(isolates)

    degree = dict(G.degree())

    # Layout
    pos = nx.spring_layout(G, seed=7, k=0.7, iterations=500, weight="shared_genes")

    node_colors = [
        state_colors.get(G.nodes[n]["dominant_state"], "lightgrey")
        for n in G.nodes()
    ]

    node_sizes = [
        250 + 90 * degree.get(n, 0)
        for n in G.nodes()
    ]

    edge_widths = [
        0.4 + 0.12 * G.edges[e]["shared_genes"]
        for e in G.edges()
    ]

    fig, ax = plt.subplots(figsize=(8, 7))

    nx.draw_networkx_edges(
        G,
        pos,
        ax=ax,
        edge_color="grey",
        width=edge_widths,
        alpha=0.35,
    )

    nx.draw_networkx_nodes(
        G,
        pos,
        ax=ax,
        node_color=node_colors,
        node_size=node_sizes,
        edgecolors="black",
        linewidths=1.0,
    )

    nx.draw_networkx_labels(
        G,
        pos,
        ax=ax,
        font_size=8,
        font_weight="bold",
        font_color="white",
    )

    # Manual legend
    handles = []
    for state, color in state_colors.items():
        handles.append(
            plt.Line2D(
                [0], [0],
                marker="o",
                color="w",
                markerfacecolor=color,
                markeredgecolor="black",
                markersize=10,
                label=state,
            )
        )

    leg1 = ax.legend(
        handles=handles,
        title="Dominant Regulon",
        loc="center left",
        bbox_to_anchor=(1.02, 0.72),
        frameon=False,
    )
    ax.add_artist(leg1)

    size_handles = [
        plt.scatter([], [], s=250, facecolors="white", edgecolors="black", label="Low (≤10 edges)"),
        plt.scatter([], [], s=800, facecolors="white", edgecolors="black", label="Medium (11–25 edges)"),
        plt.scatter([], [], s=1600, facecolors="white", edgecolors="black", label="High (>25 edges)"),
    ]

    ax.legend(
        handles=size_handles,
        title="TF Connectivity (Degree)",
        loc="center left",
        bbox_to_anchor=(1.02, 0.35),
        frameon=False,
    )

    ax.set_axis_off()
    plt.tight_layout()

    pdf = os.path.join(args.outdir, f"{args.prefix}.pdf")
    png = os.path.join(args.outdir, f"{args.prefix}.png")
    node_csv = os.path.join(args.outdir, f"{args.prefix}_nodes.csv")
    edge_csv = os.path.join(args.outdir, f"{args.prefix}_edges.csv")

    fig.savefig(pdf, dpi=300, bbox_inches="tight")
    fig.savefig(png, dpi=300, bbox_inches="tight")

    nodes_out = pd.DataFrame([
        {
            "TF": n,
            "dominant_state": G.nodes[n]["dominant_state"],
            "max_auc": G.nodes[n]["max_auc"],
            "n_targets": G.nodes[n]["n_targets"],
            "degree": degree.get(n, 0),
        }
        for n in G.nodes()
    ])
    nodes_out.to_csv(node_csv, index=False)

    edges_out = pd.DataFrame([
        {
            "TF1": u,
            "TF2": v,
            "shared_genes": G.edges[u, v]["shared_genes"],
            "jaccard": G.edges[u, v]["jaccard"],
        }
        for u, v in G.edges()
    ])
    edges_out.to_csv(edge_csv, index=False)

    print(f"Saved: {pdf}")
    print(f"Saved: {png}")
    print(f"Saved: {node_csv}")
    print(f"Saved: {edge_csv}")


if __name__ == "__main__":
    main()
