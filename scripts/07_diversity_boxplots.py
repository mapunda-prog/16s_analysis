#!/usr/bin/env python3
"""
Step 07 -- alpha diversity boxplots per region, grouped by specimen type and
by collection site (the same two groupings as 05_diversity_analysis.py's
summary tables).

Reads alpha_diversity_per_sample.tsv (written by 05_diversity_analysis.py)
and produces, per region, two 4-panel figures -- one boxplot per metric
(observed richness, Shannon, Simpson, Pielou evenness) -- one figure grouped
by type_ofsample, one by hf_name. Each group value gets a fixed color across
every region's chart (Blood is always the same color everywhere, etc.), and
each box is annotated with its sample count since group sizes vary a lot
(single digits to 20+) in this dataset.

Requires matplotlib -- run locally, like the other 05/06 scripts.

Usage:
    python3 scripts/07_diversity_boxplots.py
    python3 scripts/07_diversity_boxplots.py --alpha-table alpha_beta_diversity/alpha_diversity_per_sample.tsv
"""

import argparse
import csv
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
PROJECT = os.path.dirname(HERE)

METRICS = [
    ("observed_richness", "Observed richness"),
    ("shannon", "Shannon"),
    ("simpson", "Simpson"),
    ("pielou_evenness", "Pielou evenness"),
]
GROUP_COLUMNS = [("type_ofsample", "specimen type", ["Blood", "Swab", "Blood+Swab", "Unknown"]),
                 ("hf_name", "collection site", None)]

# dataviz skill categorical palette, light mode, fixed adjacent-pairlist order
SERIES_COLORS = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100",
                 "#e87ba4", "#008300", "#4a3aa7", "#e34948"]
SURFACE = "#fcfcfb"
INK_PRIMARY = "#0b0b0b"
INK_SECONDARY = "#52514e"
INK_MUTED = "#898781"
GRIDLINE = "#e1e0d9"
BASELINE = "#c3c2b7"


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--alpha-table",
                    default=os.path.join(PROJECT, "alpha_beta_diversity", "alpha_diversity_per_sample.tsv"))
    ap.add_argument("--outdir", default=os.path.join(PROJECT, "alpha_beta_diversity", "boxplots"))
    ap.add_argument("--regions", nargs="*", help="default: every region in --alpha-table")
    args = ap.parse_args()

    if not os.path.isfile(args.alpha_table):
        sys.exit(f"not found: {args.alpha_table}\nRun scripts/05_diversity_analysis.py first.")
    rows = read_alpha_table(args.alpha_table)
    regions = args.regions or sorted({r["region"] for r in rows})
    os.makedirs(args.outdir, exist_ok=True)

    for group_col, group_label, fixed_order in GROUP_COLUMNS:
        order = fixed_order or sorted({r[group_col] for r in rows} - {"Unknown"}) + ["Unknown"]
        order = [g for g in order if any(r[group_col] == g for r in rows)]
        color_of = {g: SERIES_COLORS[i % len(SERIES_COLORS)] for i, g in enumerate(order)}

        for region in regions:
            region_rows = [r for r in rows if r["region"] == region]
            if not region_rows:
                continue
            present = [g for g in order if any(r[group_col] == g for r in region_rows)]
            out_path = os.path.join(args.outdir, f"{region}_by_{group_col}.png")
            plot_region(region, group_col, group_label, present, color_of, region_rows, out_path)
            print(f"{region} x {group_col}: {len(present)} group(s) -> {out_path}")


def read_alpha_table(path):
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        rows = []
        for r in reader:
            for metric, _ in METRICS:
                r[metric] = float(r[metric]) if r.get(metric) not in (None, "") else None
            rows.append(r)
        return rows


def plot_region(region, group_col, group_label, groups, color_of, rows, out_path):
    fig, axes = plt.subplots(2, 2, figsize=(10, 8), dpi=150)
    fig.patch.set_facecolor(SURFACE)
    fig.suptitle(f"{region} -- alpha diversity by {group_label}", color=INK_PRIMARY, fontsize=13)

    for ax, (metric, metric_label) in zip(axes.flat, METRICS):
        ax.set_facecolor(SURFACE)
        data, counts = [], []
        for g in groups:
            vals = [r[metric] for r in rows if r[group_col] == g and r[metric] is not None]
            data.append(vals)
            counts.append(len(vals))

        bp = ax.boxplot(data, patch_artist=True, widths=0.6, showfliers=True,
                        medianprops=dict(color=INK_PRIMARY, linewidth=1.5),
                        whiskerprops=dict(color=INK_MUTED, linewidth=1),
                        capprops=dict(color=INK_MUTED, linewidth=1),
                        flierprops=dict(marker="o", markersize=3, markeredgecolor=INK_MUTED,
                                        markerfacecolor="none"))
        for patch, g in zip(bp["boxes"], groups):
            patch.set_facecolor(color_of[g])
            patch.set_edgecolor(SURFACE)
            patch.set_alpha(0.85)

        ymax = max((v for vals in data for v in vals), default=1)
        for i, (g, n) in enumerate(zip(groups, counts), start=1):
            ax.text(i, ymax * 1.06, f"n={n}", ha="center", va="bottom",
                    fontsize=7, color=INK_SECONDARY)

        ax.set_xticks(range(1, len(groups) + 1))
        ax.set_xticklabels(groups, fontsize=8, color=INK_SECONDARY, rotation=15, ha="right")
        ax.set_title(metric_label, color=INK_PRIMARY, fontsize=10)
        ax.yaxis.grid(True, color=GRIDLINE, linewidth=0.6)
        ax.set_axisbelow(True)
        for spine in ("top", "right"):
            ax.spines[spine].set_visible(False)
        for spine in ("left", "bottom"):
            ax.spines[spine].set_color(BASELINE)
        ax.tick_params(colors=INK_MUTED)
        ax.set_ylim(top=ymax * 1.15 if ymax > 0 else 1)

    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(out_path, facecolor=SURFACE)
    plt.close(fig)


if __name__ == "__main__":
    main()
