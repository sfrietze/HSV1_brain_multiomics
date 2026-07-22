#!/usr/bin/env python3

import argparse
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import logomaker
from pyjaspar import jaspardb


TF_TO_JASPAR = {
    "Irf1": "MA0050.2",
    "Irf8": "MA0652.1",
    "Jun": "MA0099.3",
    "Junb": "MA0491.1",
    "Jund": "MA0492.1",
    "Atf3": "MA0605.1",
    "Nfe2l2": "MA0150.2",
    "Cebpb": "MA0466.2",
    "Egr1": "MA0162.3",
    "Maf": "MA0496.1",
    "Mafb": "MA0497.1",
    "Klf2": "MA0036.3",
    "Klf12": "MA0596.1",
    "Pbx1": "MA0070.1",
    "Pbx3": "MA0794.1",
    "Ets1": "MA0098.3",
    "Etv6": "MA0769.1",
    "Fli1": "MA0470.1",
    "Ikzf1": "MA0164.1",
    "Ikzf3": "MA0165.1",
    "Irf2": "MA0051.1",
    "Irf9": "MA0653.1",
    "Sp100": "MA0746.1",
    "Ctcf": "MA0139.1",
    "Erg": "MA0471.1",
    "Zfp148": "MA0759.1",
}


def pwm_object_to_dataframe(pwm_obj):
    """
    Convert Biopython PositionWeightMatrix from pyjaspar to
    position x base pandas DataFrame.
    """
    pwm = pd.DataFrame({
        "A": list(pwm_obj["A"]),
        "C": list(pwm_obj["C"]),
        "G": list(pwm_obj["G"]),
        "T": list(pwm_obj["T"]),
    }).astype(float)

    pwm = pwm.div(pwm.sum(axis=1), axis=0)
    return pwm


def pwm_to_information_content(pwm):
    """
    Convert PWM probabilities to information-content matrix in bits.
    """
    ic_total = 2 + (pwm * np.log2(pwm + 1e-9)).sum(axis=1)
    ic = pwm.mul(ic_total, axis=0)
    return ic


def trim_low_information_positions(ic, threshold=0.2):
    keep = ic.sum(axis=1) > threshold
    return ic.loc[keep].reset_index(drop=True)


def plot_logo(tf, motif_id, jdb, outdir, fmt):
    motif = jdb.fetch_motif_by_id(motif_id)

    pwm = pwm_object_to_dataframe(motif.pwm)
    ic = pwm_to_information_content(pwm)
    ic = trim_low_information_positions(ic)

    if ic.shape[0] == 0:
        raise ValueError(f"{tf} / {motif_id}: empty motif after trimming")

    fig, ax = plt.subplots(figsize=(3, 1.4))
    logomaker.Logo(ic, ax=ax)

    ax.set_title(tf, fontsize=9)
    ax.set_ylim(0, 2)
    ax.set_xticks(range(len(ic)))
    ax.set_xticklabels(range(1, len(ic) + 1), fontsize=6)
    ax.set_ylabel("bits", fontsize=7)
    ax.set_xlabel("Position", fontsize=7)
    ax.spines[["top", "right"]].set_visible(False)

    plt.tight_layout()

    outfile = os.path.join(outdir, f"{tf}_motif.{fmt}")
    fig.savefig(outfile, format=fmt, bbox_inches="tight")
    plt.close(fig)

    return outfile


def main():
    parser = argparse.ArgumentParser(
        description="Generate JASPAR reference TF motif logos."
    )

    parser.add_argument(
        "--outdir",
        default="outputs/Figure4/motif_logos",
        help="Output directory"
    )

    parser.add_argument(
        "--format",
        default="svg",
        choices=["svg", "pdf", "png"],
        help="Output file format"
    )

    parser.add_argument(
        "--tfs",
        default="Irf9,Mafb,Jun,Ctcf,Atf3,Ikzf1,Fli1",
        help="Comma-separated TF list"
    )

    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    selected_tfs = [x.strip() for x in args.tfs.split(",") if x.strip()]
    jdb = jaspardb()

    generated = []
    failed = []

    for tf in selected_tfs:
        if tf not in TF_TO_JASPAR:
            failed.append({
                "TF": tf,
                "JASPAR_ID": "",
                "reason": "TF not present in TF_TO_JASPAR"
            })
            continue

        motif_id = TF_TO_JASPAR[tf]

        try:
            outfile = plot_logo(tf, motif_id, jdb, args.outdir, args.format)
            generated.append({
                "TF": tf,
                "JASPAR_ID": motif_id,
                "output_file": outfile
            })
            print(f"Saved: {outfile}")

        except Exception as e:
            failed.append({
                "TF": tf,
                "JASPAR_ID": motif_id,
                "reason": str(e)
            })

    generated_df = pd.DataFrame(generated)
    failed_df = pd.DataFrame(failed)

    generated_df.to_csv(
        os.path.join(args.outdir, "motif_logos_generated.csv"),
        index=False
    )

    failed_df.to_csv(
        os.path.join(args.outdir, "motif_logos_failed.csv"),
        index=False
    )

    print("\nMotifs generated:")
    print(generated_df)

    if len(failed_df) > 0:
        print("\nMissing / failed:")
        print(failed_df)


if __name__ == "__main__":
    main()
