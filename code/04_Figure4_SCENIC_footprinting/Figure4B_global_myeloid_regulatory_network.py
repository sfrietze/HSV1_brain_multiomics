#!/usr/bin/env python3

import argparse
import os
import numpy as np
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
import mudata as md
from matplotlib.lines import Line2D


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--auc", required=True, help="AUCell_extended.h5mu")
    parser.add_argument("--eregulons", required=True, help="eRegulons_extended.tsv")
    parser.add_argument("--outdir", default="outputs/Figure4")
    parser.add_argument("--prefix", default="Figure4B_global_myeloid_regulatory_network")
    parser.add_argument("--threshold", type=float, default=1.0)
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    # ------------------------------------------------------------
    # Load AUCell extended object
    # ------------------------------------------------------------
    mdata = md.read_h5mu(args.auc)

    if "extended_gene_based_AUC" in mdata.mod:
        auc = mdata["extended_gene_based_AUC"].to_df()
    else:
        raise ValueError("Missing modality: extended_gene_based_AUC")

    meta = mdata.obs.copy()

    # ------------------------------------------------------------
    # Find celltype column
    # ------------------------------------------------------------
    possible_cols = [
        "celltype",
        "scRNA_counts:celltype",
        "scRNA_counts:annotated_clusters",
        "annotated_clusters",
    ]

    celltype_col = None
    for c in possible_cols:
        if c in meta.columns:
            celltype_col = c
            break

    if celltype_col is None:
        raise ValueError(f"No celltype column found. Available obs columns: {list(meta.columns)}")

    auc["celltype"] = meta[celltype_col].astype(str).values

    # Normalize labels to match old code
    auc["celltype"] = (
        auc["celltype"]
        .str.replace(" ", "_", regex=False)
        .str.replace("-", "-", regex=False)
        .str.strip()
    )

    # ------------------------------------------------------------
    # Restrict to selected myeloid clusters
    # ------------------------------------------------------------
    myeloid_clusters = [
        "Homeostatic_Microglia",
        "Transiently_Activated_Microglia",
        "IFN-Responsive_Microglia",
        "Primed_Microglia",
        "Mitochondrial-Activated_Microglia",
        "Infiltrating_Macrophages",
    ]

    # If labels include condition suffix, remove it
    auc["celltype"] = (
        auc["celltype"]
        .str.replace("_Mock", "", regex=False)
        .str.replace("_HSV1", "", regex=False)
    )

    df_myeloid = auc[auc["celltype"].isin(myeloid_clusters)].copy()

    if df_myeloid.shape[0] == 0:
        print("Observed celltype labels:")
        print(auc["celltype"].value_counts().head(30))
        raise ValueError("No selected myeloid cells found after label cleanup.")

    # ------------------------------------------------------------
    # Extract extended regulon columns
    # ------------------------------------------------------------
    reg_cols = [c for c in df_myeloid.columns if "extended" in c]

    if len(reg_cols) == 0:
        print("AUC columns:")
        print(df_myeloid.columns.tolist()[:50])
        raise ValueError("No extended regulon columns found.")

    reg_to_tf = {reg: reg.split("_")[0].upper() for reg in reg_cols}

    # ------------------------------------------------------------
    # Compute dominant cluster via Z-score
    # ------------------------------------------------------------
    tf_cluster_matrix = df_myeloid.groupby("celltype")[reg_cols].mean()

    dominant_state = {}
    dominant_score = {}

    for reg in reg_cols:
        tf = reg_to_tf[reg]
        cluster_means = tf_cluster_matrix[reg]

        if cluster_means.std() == 0:
            continue

        zscores = (cluster_means - cluster_means.mean()) / cluster_means.std()

        dominant_state[tf] = zscores.idxmax()
        dominant_score[tf] = zscores.max()

    # ------------------------------------------------------------
    # Build TF network from extended eRegulon target overlap
    # ------------------------------------------------------------
    extended = pd.read_csv(args.eregulons, sep="\t")
    extended["TF"] = extended["TF"].astype(str).str.upper()

    G = nx.Graph()

    for tf in extended["TF"].unique():
        G.add_node(tf)

    tf_targets = extended.groupby("TF")["Gene"].apply(set).to_dict()
    tfs = list(tf_targets.keys())

    for i, tf1 in enumerate(tfs):
        for tf2 in tfs[i + 1:]:
            overlap = len(tf_targets[tf1] & tf_targets[tf2])
            if overlap > 0:
                G.add_edge(tf1, tf2, weight=overlap)

    print("TF nodes:", len(G.nodes))
    print("TF edges:", len(G.edges))

    degree_dict = dict(G.degree())
    node_sizes = [max(degree_dict[n] * 80, 80) for n in G.nodes]

    # ------------------------------------------------------------
    # Palette
    # ------------------------------------------------------------
    cluster_palette = {
        "Homeostatic_Microglia": "darkorange",
        "Transiently_Activated_Microglia": "#fcae91",
        "IFN-Responsive_Microglia": "indigo",
        "Primed_Microglia": "#cb181d",
        "Mitochondrial-Activated_Microglia": "#99000d",
        "Infiltrating_Macrophages": "#67000d",
    }

    node_colors = []
    for n in G.nodes:
        if n in dominant_state and dominant_score.get(n, 0) >= args.threshold:
            cluster_label = dominant_state[n]
            color = cluster_palette.get(cluster_label, "#cccccc")
        else:
            color = "#cccccc"
        node_colors.append(color)

    # ------------------------------------------------------------
    # Fixed layout
    # ------------------------------------------------------------
    pos = nx.spring_layout(G, seed=42, k=0.5)

    # ------------------------------------------------------------
    # Plot
    # ------------------------------------------------------------
    plt.figure(figsize=(10, 8))

    nx.draw_networkx_edges(G, pos, alpha=0.25)

    nx.draw_networkx_nodes(
        G,
        pos,
        node_color=node_colors,
        node_size=node_sizes,
        edgecolors="black",
        linewidths=0.6,
    )

    label_colors = {}

    for n in G.nodes:
        if n in dominant_state and dominant_score.get(n, 0) >= args.threshold:
            cluster = dominant_state[n]
            if cluster in [
                "Infiltrating_Macrophages",
                "Primed_Microglia",
                "IFN-Responsive_Microglia",
                "Mitochondrial-Activated_Microglia",
            ]:
                label_colors[n] = "white"
            else:
                label_colors[n] = "black"
        else:
            label_colors[n] = "black"

    for node, (x, y) in pos.items():
        plt.text(
            x,
            y,
            s=node,
            fontsize=9,
            fontweight="bold",
            ha="center",
            va="center",
            color=label_colors[node],
        )

    plt.axis("off")

    # ------------------------------------------------------------
    # Dominant regulon legend
    # ------------------------------------------------------------
    legend_elements = []

    states_present = sorted(
        set([
            dominant_state[n]
            for n in G.nodes
            if n in dominant_state and dominant_score.get(n, 0) >= args.threshold
        ])
    )

    for state in states_present:
        legend_elements.append(
            Line2D(
                [0], [0],
                marker="o",
                color="w",
                label=state.replace("_", " "),
                markerfacecolor=cluster_palette.get(state, "lightgray"),
                markeredgecolor="black",
                markersize=9,
            )
        )

    legend1 = plt.legend(
        handles=legend_elements,
        loc="upper left",
        bbox_to_anchor=(1.02, 1),
        frameon=False,
        title="Dominant Regulon",
    )

    plt.gca().add_artist(legend1)

    # ------------------------------------------------------------
    # Connectivity legend
    # ------------------------------------------------------------
    min_deg = min(degree_dict.values())
    mid_deg = int(np.median(list(degree_dict.values())))
    max_deg = max(degree_dict.values())

    size_legend = [
        plt.scatter([], [], s=max(min_deg * 80, 80), edgecolors="black",
                    facecolors="white", label=f"{min_deg} connections"),
        plt.scatter([], [], s=mid_deg * 80, edgecolors="black",
                    facecolors="white", label=f"{mid_deg} connections"),
        plt.scatter([], [], s=max_deg * 80, edgecolors="black",
                    facecolors="white", label=f"{max_deg} connections"),
    ]

    plt.legend(
        handles=size_legend,
        title="TF connectivity",
        loc="lower left",
        bbox_to_anchor=(1.02, 0.3),
        frameon=False,
    )

    plt.tight_layout()

    pdf = os.path.join(args.outdir, f"{args.prefix}.pdf")
    png = os.path.join(args.outdir, f"{args.prefix}.png")

    plt.savefig(pdf, bbox_inches="tight")
    plt.savefig(png, dpi=300, bbox_inches="tight")

    node_table = pd.DataFrame({
        "TF": list(G.nodes),
        "degree": [degree_dict[n] for n in G.nodes],
        "dominant_state": [dominant_state.get(n, "none") for n in G.nodes],
        "dominant_score": [dominant_score.get(n, np.nan) for n in G.nodes],
    })

    edge_table = pd.DataFrame([
        {"TF1": u, "TF2": v, "shared_targets": d["weight"]}
        for u, v, d in G.edges(data=True)
    ])

    node_table.to_csv(os.path.join(args.outdir, f"{args.prefix}_nodes.csv"), index=False)
    edge_table.to_csv(os.path.join(args.outdir, f"{args.prefix}_edges.csv"), index=False)

    print(f"Saved: {pdf}")
    print(f"Saved: {png}")
    print("Saved node and edge tables.")


if __name__ == "__main__":
    main()
