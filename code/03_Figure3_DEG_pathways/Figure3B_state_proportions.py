#!/usr/bin/env python3

import argparse
import os

import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, help="Input annotated h5ad")
    parser.add_argument("--outdir", default="outputs/Figure3")
    parser.add_argument("--prefix", default="Figure3B_state_proportions")
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    adata = sc.read(args.input)

    cluster_col = "annotated_clusters"
    condition_col = "condition"

    for col in [cluster_col, condition_col]:
        if col not in adata.obs.columns:
            raise ValueError(f"Missing required adata.obs column: {col}")

    obs = adata.obs[[cluster_col, condition_col]].copy()
    obs[cluster_col] = (
        obs[cluster_col]
        .astype(str)
        .str.replace("‚Äì", "-", regex=False)
        .str.strip()
    )
    obs[condition_col] = obs[condition_col].astype(str).str.strip()

    myeloid_states = [
        "Homeostatic Microglia",
        "IFN-Responsive Microglia",
        "Transiently Activated Microglia",
        "IEG-High Microglia",
        "Primed Microglia",
        "Mitochondrial-Activated Microglia",
        "Infiltrating Macrophages",
    ]

    palette = {
        "Homeostatic Microglia": "darkorange",
        "IFN-Responsive Microglia": "#4B0082",
        "Transiently Activated Microglia": "#fcae91",
        "IEG-High Microglia": "#fb6a4a",
        "Primed Microglia": "#cb181d",
        "Mitochondrial-Activated Microglia": "#99000d",
        "Infiltrating Macrophages": "#67000d",
    }

    condition_order = ["Mock", "HSV1"]

    obs = obs[
        obs[cluster_col].isin(myeloid_states)
        & obs[condition_col].isin(condition_order)
    ].copy()

    counts = (
        obs.groupby([condition_col, cluster_col], observed=False)
        .size()
        .unstack(fill_value=0)
        .reindex(index=condition_order, columns=myeloid_states, fill_value=0)
    )

    fractions = counts.div(counts.sum(axis=1), axis=0)

    counts.to_csv(os.path.join(args.outdir, f"{args.prefix}_counts.csv"))
    fractions.to_csv(os.path.join(args.outdir, f"{args.prefix}_fractions.csv"))

    fig, ax = plt.subplots(figsize=(6.5, 1.8))

    y_positions = range(len(condition_order))
    left = [0.0] * len(condition_order)

    for state in myeloid_states:
        vals = fractions[state].values
        ax.barh(
            y_positions,
            vals,
            left=left,
            height=0.55,
            color=palette[state],
            edgecolor="white",
            linewidth=0.6,
            label=state,
        )
        left = [l + v for l, v in zip(left, vals)]

    ax.set_yticks(list(y_positions))
    ax.set_yticklabels(["Mock", "HSV-1"])
    ax.set_xlim(0, 1)
    ax.set_xlabel("Fraction of myeloid cells")
    ax.set_xticks([0, 0.5, 1.0])
    ax.invert_yaxis()

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    handles = [
        mpatches.Patch(color=palette[state], label=state)
        for state in myeloid_states
    ]

    ax.legend(
        handles=handles,
        frameon=False,
        loc="center left",
        bbox_to_anchor=(1.03, 0.5),
        fontsize=7,
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
