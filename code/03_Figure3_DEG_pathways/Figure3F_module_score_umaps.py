#!/usr/bin/env python3

import argparse
import os
import textwrap

import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable


def load_gmt(gmt_path):
    gene_sets = {}
    with open(gmt_path) as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 3:
                gene_sets[parts[0]] = parts[2:]
    return gene_sets


def match_genes(adata, genes):
    adata_gene_map = {g.upper(): g for g in adata.var_names}
    matched = []
    for g in genes:
        gu = str(g).upper()
        if gu in adata_gene_map:
            matched.append(adata_gene_map[gu])
    return sorted(set(matched))


def score_pathway(adata, gene_sets, pathway, score_key):
    if pathway not in gene_sets:
        raise KeyError(f"{pathway} not found in GMT")

    genes = match_genes(adata, gene_sets[pathway])

    print(f"{pathway}: {len(genes)} genes matched in adata")
    if len(genes) < 5:
        raise ValueError(f"Too few genes matched for {pathway}: {len(genes)}")

    if score_key not in adata.obs.columns:
        sc.tl.score_genes(
            adata,
            gene_list=genes,
            score_name=score_key,
            use_raw=False,
        )

    return genes


def plot_panel(adata, score_keys, row_labels, outdir, prefix):
    X = adata.obsm["X_umap"]

    cond = adata.obs["condition"].astype(str).str.strip()
    is_mock = cond.isin(["Mock", "MOCK"]).values
    is_hsv = cond.isin(["HSV1", "HSV-1", "HSV"]).values

    if is_mock.sum() == 0 or is_hsv.sum() == 0:
        raise ValueError(f"Condition counts look wrong: Mock={is_mock.sum()}, HSV1={is_hsv.sum()}")

    fig, axes = plt.subplots(
        nrows=len(score_keys),
        ncols=2,
        figsize=(6.2, 7.4),
        sharex=True,
        sharey=True,
    )

    point_size = 4.5
    bg_size = 2.5

    for i, (score_key, row_label) in enumerate(zip(score_keys, row_labels)):
        score = adata.obs[score_key].astype(float)

        fixed_scales = {
            "GOBP_DEFENSE_RESPONSE_TO_VIRUS_score": (-0.2, 0.6),
            "GOBP_INTERLEUKIN_1_PRODUCTION_score": (0.0, 0.25),
            "KEGG_TOLL_LIKE_RECEPTOR_SIGNALING_PATHWAY_score": (-0.25, 1.0),
        }

        vmin, vmax = fixed_scales.get(
            score_key,
            (np.nanquantile(score, 0.01), np.nanquantile(score, 0.99))
        )

        for j, (mask, title) in enumerate([(is_mock, "Mock"), (is_hsv, "HSV-1")]):
            ax = axes[i, j]

            ax.scatter(
                X[~mask, 0],
                X[~mask, 1],
                c="#d9d9d9",
                s=bg_size,
                alpha=0.35,
                linewidths=0,
                rasterized=True,
            )

            im = ax.scatter(
                X[mask, 0],
                X[mask, 1],
                c=score.loc[mask],
                cmap="magma",
                s=point_size,
                linewidths=0,
                vmin=vmin,
                vmax=vmax,
                rasterized=True,
            )

            ax.set_xticks([])
            ax.set_yticks([])
            ax.set_frame_on(False)

            if i == 0:
                ax.set_title(title, fontsize=10, fontweight="bold")

            if j == 1:
                divider = make_axes_locatable(ax)
                cax = divider.append_axes("right", size="4%", pad=0.03)
                cb = fig.colorbar(im, cax=cax)
                cb.ax.tick_params(labelsize=6, length=2)
                cb.set_label(
                    "\n".join(textwrap.wrap(row_label, width=16)),
                    fontsize=10,
                    rotation=90,
                    labelpad=10,
                )

    plt.tight_layout(w_pad=0.4, h_pad=0.4)

    pdf = os.path.join(outdir, f"{prefix}.pdf")
    png = os.path.join(outdir, f"{prefix}.png")

    fig.savefig(pdf, dpi=600, bbox_inches="tight")
    fig.savefig(png, dpi=600, bbox_inches="tight")
    print(f"Saved: {pdf}")
    print(f"Saved: {png}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, help="Annotated h5ad")
    parser.add_argument("--go_gmt", required=True, help="GO BP GMT file")
    parser.add_argument("--kegg_gmt", required=True, help="KEGG/CP GMT file")
    parser.add_argument("--outdir", default="outputs/Figure3")
    parser.add_argument("--prefix", default="Figure3F_module_score_umaps")
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    adata = sc.read_h5ad(args.input)

    if "X_umap" not in adata.obsm:
        raise ValueError("Missing adata.obsm['X_umap']")
    if "condition" not in adata.obs.columns:
        raise ValueError("Missing adata.obs['condition']")

    go = load_gmt(args.go_gmt)
    kegg = load_gmt(args.kegg_gmt)

    pathways = [
        ("GOBP_DEFENSE_RESPONSE_TO_VIRUS", go, "Defense Response\nTo Virus"),
        ("GOBP_INTERLEUKIN_1_PRODUCTION", go, "IL1 Production"),
        ("KEGG_TOLL_LIKE_RECEPTOR_SIGNALING_PATHWAY", kegg, "TLR Signaling"),
    ]

    score_keys = []
    labels = []

    for pathway, gmt, label in pathways:
        score_key = f"{pathway}_score"
        score_pathway(adata, gmt, pathway, score_key)
        score_keys.append(score_key)
        labels.append(label)

    plot_panel(
        adata=adata,
        score_keys=score_keys,
        row_labels=labels,
        outdir=args.outdir,
        prefix=args.prefix,
    )

    adata.obs[score_keys].to_csv(
        os.path.join(args.outdir, f"{args.prefix}_scores.csv")
    )


if __name__ == "__main__":
    main()
