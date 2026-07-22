#!/usr/bin/env python3

import argparse
import os
import re
import pandas as pd
import mudata as md


def clean_tf(x):
    x = str(x)
    x = re.sub(r"\(\+\)|\(-\)", "", x)
    x = x.replace("_extended", "")
    x = x.replace("_", "")
    return x.lower()


def main():

    parser = argparse.ArgumentParser(
        description="Extract SCENIC+ motif metadata for manuscript figures."
    )

    parser.add_argument(
        "--input",
        required=True,
        help="SCENIC+ MuData (.h5mu)"
    )

    parser.add_argument(
        "--outdir",
        default="outputs/Figure4"
    )

    parser.add_argument(
        "--prefix",
        default="Figure4B_motif_metadata"
    )

    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    print("\nLoading MuData...")
    mdata = md.read_h5mu(args.input)

    print("\nuns keys:")
    for k in sorted(mdata.uns.keys()):
        print("  ", k)

    if "direct_e_regulon_metadata" not in mdata.uns:
        raise KeyError(
            "direct_e_regulon_metadata not found."
        )

    meta = mdata.uns["direct_e_regulon_metadata"].copy()

    print("\nMetadata dimensions:")
    print(meta.shape)

    print("\nMetadata columns:")
    for c in meta.columns:
        print("  ", c)

    tf_column = None

    for candidate in [
        "TF",
        "TF_name",
        "Transcription_factor",
        "eRegulon_name",
    ]:
        if candidate in meta.columns:
            tf_column = candidate
            break

    if tf_column is None:
        raise ValueError("Could not determine TF column.")

    motif_columns = [
        c for c in meta.columns
        if any(
            x in c.lower()
            for x in [
                "motif",
                "logo",
                "matrix",
                "jaspar",
                "cb",
                "annotation",
                "pwm",
                "cluster"
            ]
        )
    ]

    print("\nMotif-related columns:")
    print(motif_columns)

    selected = [
        "Irf9",
        "Mafb",
        "Jun",
        "Ctcf",
        "Atf3",
        "Ikzf1",
        "Fli1"
    ]

    meta["TF_clean"] = meta[tf_column].map(clean_tf)

    rows = []

    print("\n=====================================")
    print("Selected TF matches")
    print("=====================================\n")

    for tf in selected:

        hit = meta[
            meta["TF_clean"].str.contains(
                clean_tf(tf),
                regex=False,
                na=False
            )
        ]

        if hit.empty:
            print(f"{tf}: NOT FOUND\n")
            continue

        row = hit.iloc[0]

        rows.append(row)

        print(f"{tf}")
        print("-" * len(tf))

        print(row)

        print()

    out_csv = os.path.join(
        args.outdir,
        f"{args.prefix}_selected.csv"
    )

    pd.DataFrame(rows).to_csv(out_csv, index=False)

    full_csv = os.path.join(
        args.outdir,
        f"{args.prefix}_full.csv"
    )

    meta.to_csv(full_csv, index=False)

    print("Saved:")
    print(out_csv)
    print(full_csv)


if __name__ == "__main__":
    main()

