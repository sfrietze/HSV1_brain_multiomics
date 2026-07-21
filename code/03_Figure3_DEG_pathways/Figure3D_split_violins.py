#!/usr/bin/env python3

import argparse
import os

import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, help="Input annotated h5ad")
    parser.add_argument("--outdir", default="outputs/Figure3")
    parser.add_argument("--prefix", default="Figure3D_split_violins")
    parser.add_argument("--cluster", default="Infiltrating Macrophages")
    parser.add_argument("--n_genes", type=int, default=4)
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

    adata_sub = adata[
        (adata.obs["annotated_clusters"] == args.cluster)
        & (adata.obs["condition"].isin(["Mock", "HSV1"]))
    ].copy()

    adata_sub.obs["condition"] = adata_sub.obs["condition"].astype("category")
    adata_sub.obs["condition"] = adata_sub.obs["condition"].cat.set_categories(
        ["Mock", "HSV1"],
        ordered=True,
    )

    mock_blue = "#377eb8"
    hsv1_red = "#d7301f"
    adata_sub.uns["condition_colors"] = [mock_blue, hsv1_red]

    sc.tl.rank_genes_groups(
        adata_sub,
        groupby="condition",
        groups=["Mock"],
        reference="HSV1",
        method="wilcoxon",
        key_added="deg_mock_vs_hsv1",
    )

    de_df = sc.get.rank_genes_groups_df(
        adata_sub,
        group="Mock",
        key="deg_mock_vs_hsv1",
    )

    de_df.to_csv(
        os.path.join(args.outdir, f"{args.prefix}_InfiltratingMacrophages_DEG_table.csv"),
        index=False,
    )

    mock_up = (
        de_df.query("logfoldchanges > 0")
        .sort_values("pvals_adj")
        .head(args.n_genes)["names"]
        .tolist()
    )

    hsv1_up = (
        de_df.query("logfoldchanges < 0")
        .sort_values("pvals_adj")
        .head(args.n_genes)["names"]
        .tolist()
    )

    pd.DataFrame({
        "Mock_up_genes": pd.Series(mock_up),
        "HSV1_up_genes": pd.Series(hsv1_up),
    }).to_csv(
        os.path.join(args.outdir, f"{args.prefix}_selected_genes.csv"),
        index=False,
    )

    mpl.rcParams["axes.prop_cycle"] = mpl.cycler(color=[mock_blue, hsv1_red])

    def plot_split_violin(genes, label, outfile):
        with plt.rc_context({"figure.figsize": (5.2, 1.35)}):
            sc.pl.rank_genes_groups_violin(
                adata_sub,
                groups=["Mock"],
                gene_names=genes,
                key="deg_mock_vs_hsv1",
                jitter=False,
                split=True,
                scale="width",
                show=False,
            )

            fig = plt.gcf()
            axes = fig.axes

            for ax in axes:
                ax.set_xlabel("")
                ax.tick_params(axis="x", labelrotation=55, labelsize=8)
                ax.tick_params(axis="y", labelsize=8)
                ax.spines["top"].set_visible(True)
                ax.spines["right"].set_visible(True)
                ax.spines["left"].set_linewidth(0.8)
                ax.spines["bottom"].set_linewidth(0.8)
                ax.spines["top"].set_linewidth(0.8)
                ax.spines["right"].set_linewidth(0.8)

            axes[0].set_ylabel("Expression", fontsize=9)
            fig.suptitle("")
            fig.tight_layout(pad=0.25)

            fig.savefig(
                os.path.join(args.outdir, outfile),
                dpi=300,
                bbox_inches="tight",
            )
            plt.close(fig)

    plot_split_violin(
        mock_up,
        "Mock-up",
        f"{args.prefix}_Mock_up_top{args.n_genes}.pdf",
    )

    plot_split_violin(
        hsv1_up,
        "HSV1-up",
        f"{args.prefix}_HSV1_up_top{args.n_genes}.pdf",
    )

    mpl.rcParams.update(mpl.rcParamsDefault)


if __name__ == "__main__":
    main()
