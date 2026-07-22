#!/usr/bin/env python3

import argparse
import os
import pandas as pd
import mudata as md

from scenicplus.plotting.dotplot import heatmap_dotplot
from plotnine import theme, element_text


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, help="SCENIC+ MuData (.h5mu)")
    parser.add_argument("--outdir", default="outputs/Figure4")
    parser.add_argument("--prefix", default="Figure4A_SCENICplus_heatmap")
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    group_var = "scRNA_counts:celltype_condition"

    keep_groups = [
        "Homeostatic_Microglia_Mock",
        "Homeostatic_Microglia_HSV1",
        "Transiently_Activated_Microglia_Mock",
        "Transiently_Activated_Microglia_HSV1",
        "Primed_Microglia_Mock",
        "Primed_Microglia_HSV1",
        "IEG-High_Microglia_Mock",
        "IEG-High_Microglia_HSV1",
        "Infiltrating_Macrophages_Mock",
        "Infiltrating_Macrophages_HSV1",
        "IFN-Responsive_Microglia_HSV1",
    ]

    scplus_mdata = md.read_h5mu(args.input)

    mask = scplus_mdata.obs[group_var].isin(keep_groups)
    scplus_mdata_filt = scplus_mdata[mask].copy()

    # restore metadata needed after MuData subsetting
    scplus_mdata_filt.uns["direct_e_regulon_metadata"] = \
        scplus_mdata.uns["direct_e_regulon_metadata"]

    group_order = [
        g for g in keep_groups
        if g in scplus_mdata_filt.obs[group_var].astype(str).unique()
    ]

    scplus_mdata_filt.obs[group_var] = pd.Categorical(
        scplus_mdata_filt.obs[group_var],
        categories=group_order,
        ordered=True,
    )

    p = heatmap_dotplot(
        scplus_mudata=scplus_mdata_filt,
        color_modality="direct_gene_based_AUC",
        size_modality="direct_region_based_AUC",
        group_variable=group_var,
        group_variable_order=group_order,
        eRegulon_metadata_key="direct_e_regulon_metadata",
        color_feature_key="Gene_signature_name",
        size_feature_key="Region_signature_name",
        feature_name_key="eRegulon_name",
        sort_data_by="direct_gene_based_AUC",
        figsize=(9, 8),
        orientation="vertical",
        split_repressor_activator=False,
    )

    p += theme(axis_text_x=element_text(rotation=90, ha="center"))

    pdf = os.path.join(args.outdir, f"{args.prefix}.pdf")
    png = os.path.join(args.outdir, f"{args.prefix}.png")

    p.save(pdf, width=9, height=8, dpi=300)
    p.save(png, width=9, height=8, dpi=300)

    pd.Series(group_order, name="group_order").to_csv(
        os.path.join(args.outdir, f"{args.prefix}_group_order.csv"),
        index=False,
    )

    print(f"Saved: {pdf}")
    print(f"Saved: {png}")


if __name__ == "__main__":
    main()
