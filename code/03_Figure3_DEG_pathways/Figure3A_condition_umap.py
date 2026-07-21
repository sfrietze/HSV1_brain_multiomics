#!/usr/bin/env python3

import argparse
import os

import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, help="Input annotated h5ad")
    parser.add_argument("--outdir", default="outputs/Figure3")
    parser.add_argument("--prefix", default="Figure3A_condition_umap")
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    adata = sc.read(args.input)

    required_obs = ["annotated_clusters", "condition"]
    for col in required_obs:
        if col not in adata.obs.columns:
            raise ValueError(f"Missing required adata.obs column: {col}")

    if "X_umap" not in adata.obsm:
        raise ValueError("Missing adata.obsm['X_umap']")

    adata.obs["annotated_clusters"] = (
        adata.obs["annotated_clusters"]
        .astype(str)
        .str.replace("‚Äì", "-", regex=False)
        .str.strip()
    )

    adata.obs["condition"] = adata.obs["condition"].astype(str).str.strip()

    X = adata.obsm["X_umap"]

    is_mock = adata.obs["condition"].values == "Mock"
    is_hsv = adata.obs["condition"].values == "HSV1"

    is_ifn = adata.obs["annotated_clusters"].values == "IFN-Responsive Microglia"
    is_mac = adata.obs["annotated_clusters"].values == "Infiltrating Macrophages"

    dot_size = 8
    point_size = dot_size**2 / 10

    background = "#d9d9d9"
    ifn_color = "indigo"
    mac_color = "#67000d"

    fig, axes = plt.subplots(1, 2, figsize=(7, 3), sharex=True, sharey=True)

    axes[0].scatter(
        X[is_mock & ~is_ifn & ~is_mac, 0],
        X[is_mock & ~is_ifn & ~is_mac, 1],
        s=point_size,
        c=background,
        rasterized=True,
        linewidths=0,
    )

    axes[0].scatter(
        X[is_mock & is_mac, 0],
        X[is_mock & is_mac, 1],
        s=point_size * 1.2,
        c=mac_color,
        rasterized=True,
        linewidths=0,
    )

    axes[0].scatter(
        X[is_mock & is_ifn, 0],
        X[is_mock & is_ifn, 1],
        s=point_size * 1.4,
        c=ifn_color,
        rasterized=True,
        linewidths=0,
    )

    axes[0].set_title("Mock")

    axes[1].scatter(
        X[is_hsv & ~is_ifn & ~is_mac, 0],
        X[is_hsv & ~is_ifn & ~is_mac, 1],
        s=point_size,
        c=background,
        rasterized=True,
        linewidths=0,
    )

    axes[1].scatter(
        X[is_hsv & is_mac, 0],
        X[is_hsv & is_mac, 1],
        s=point_size * 1.2,
        c=mac_color,
        rasterized=True,
        linewidths=0,
    )

    axes[1].scatter(
        X[is_hsv & is_ifn, 0],
        X[is_hsv & is_ifn, 1],
        s=point_size * 1.4,
        c=ifn_color,
        rasterized=True,
        linewidths=0,
    )

    axes[1].set_title("HSV-1")

    for ax in axes:
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_frame_on(False)

    plt.tight_layout()

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

    summary = pd.DataFrame({
        "condition": ["Mock", "HSV1"],
        "total_cells": [int(is_mock.sum()), int(is_hsv.sum())],
        "IFN_responsive_microglia": [
            int((is_mock & is_ifn).sum()),
            int((is_hsv & is_ifn).sum()),
        ],
        "infiltrating_macrophages": [
            int((is_mock & is_mac).sum()),
            int((is_hsv & is_mac).sum()),
        ],
    })

    summary.to_csv(
        os.path.join(args.outdir, f"{args.prefix}_cell_counts.csv"),
        index=False,
    )


if __name__ == "__main__":
    main()
