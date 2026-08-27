#!/usr/bin/env python3
"""
Step 12 -- alpha rarefaction curves per region: expected ASV richness as a
function of sequencing depth, one line per sample.

Uses Hurlbert's (1971) exact rarefaction formula rather than repeated random
subsampling -- deterministic, and cheap even at depths in the millions:

    E[S(n)] = S_total - sum_i C(N - N_i, n) / C(N, n)

for a sample with N total reads, S_total observed ASVs, and per-ASV counts
N_i, computed in log space via math.lgamma for numerical stability. This is
the same formula behind R's vegan::rarefy and QIIME's alpha rarefaction.

Reads each region's dada2/ASV_table.tsv + taxonomy via lib/taxonomy.py
(already excludes Mitochondria/Chloroplast, consistent with every other
compiled output). Lines are colored by specimen type (config/metadata.tsv),
using the same fixed color-per-type mapping as 07_diversity_boxplots.py. A
dashed reference line marks the rarefaction depth 05_diversity_analysis.py
would use for that region (the lowest depth among samples with >=1,000
reads) -- where a sample's curve is still climbing steeply at that line, its
diversity numbers from step 5 are on a less flat part of the curve than
samples further right.

Requires matplotlib.

Usage:
    python3 scripts/12_rarefaction_plots.py --regions V3V4 V4V5 ITS
"""

import argparse
import csv
import math
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
PROJECT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(HERE, "lib"))
from taxonomy import load_region_asvs  # noqa: E402

MIN_READS = 1000  # same floor used by build_metadata.py / 05_diversity_analysis.py
N_DEPTH_POINTS = 20
TYPE_ORDER = ["Blood", "Swab", "Blood+Swab", "Unknown"]
SERIES_COLORS = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100"]
SURFACE = "#fcfcfb"
INK_PRIMARY = "#0b0b0b"
INK_SECONDARY = "#52514e"
INK_MUTED = "#898781"
GRIDLINE = "#e1e0d9"
BASELINE = "#c3c2b7"


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--results-dir", default=os.path.join(PROJECT, "results_by_region"))
    ap.add_argument("--metadata", default=os.path.join(PROJECT, "config", "metadata.tsv"))
    ap.add_argument("--regions", nargs="+", required=True, help="e.g. --regions V3V4 V4V5 ITS")
    ap.add_argument("--outdir", default=os.path.join(PROJECT, "alpha_beta_diversity", "rarefaction_curves"))
    args = ap.parse_args()

    metadata = read_metadata(args.metadata)
    os.makedirs(args.outdir, exist_ok=True)
    color_of = {g: SERIES_COLORS[i] for i, g in enumerate(TYPE_ORDER)}

    for region in args.regions:
        region_dir = os.path.join(args.results_dir, region)
        try:
            samples, asvs = load_region_asvs(region_dir)
        except FileNotFoundError as e:
            print(f"{region}: WARNING {e}, skipped")
            continue

        per_sample_counts = {s: [] for s in samples}
        for asv in asvs:
            for s, c in asv["counts"].items():
                if c > 0:
                    per_sample_counts[s].append(c)
        depths_by_sample = {s: sum(c) for s, c in per_sample_counts.items()}

        max_depth = max(depths_by_sample.values(), default=0)
        if max_depth == 0:
            print(f"{region}: no reads after contaminant filtering, skipped")
            continue
        depth_grid = sorted(set(int(x) for x in np.linspace(1, max_depth, N_DEPTH_POINTS)))

        qualifying = [n for n in depths_by_sample.values() if n >= MIN_READS]
        ref_depth = min(qualifying) if qualifying else None

        out_path = os.path.join(args.outdir, f"{region}_rarefaction.png")
        plot_region(region, samples, per_sample_counts, depths_by_sample, depth_grid,
                   ref_depth, metadata, color_of, out_path)
        note = f"ref depth (>=1000-read min, n={len(qualifying)} samples qualify) = {ref_depth}" \
            if ref_depth else "no samples reach the 1,000-read floor -- no reference line"
        print(f"{region}: {len(samples)} samples, max depth {max_depth} -- {note} -> {out_path}")


def expected_richness(counts, N, n):
    """Hurlbert's exact rarefaction: expected # distinct ASVs in a random
    subsample of size n (without replacement) from a sample with per-ASV
    counts `counts` summing to N. None if n > N (undefined beyond total depth)."""
    if n > N:
        return None
    if n <= 0:
        return 0.0
    log_denom = math.lgamma(N - n + 1) - math.lgamma(N + 1)
    s_total = len(counts)
    p_absent_sum = 0.0
    for Ni in counts:
        if N - Ni < n:
            continue  # species guaranteed present -> contributes 0
        log_num = math.lgamma(N - Ni + 1) - math.lgamma(N - Ni - n + 1)
        p_absent_sum += math.exp(log_num + log_denom)
    return s_total - p_absent_sum


def plot_region(region, samples, per_sample_counts, depths_by_sample, depth_grid,
                ref_depth, metadata, color_of, out_path):
    fig, ax = plt.subplots(figsize=(9, 6), dpi=150)
    fig.patch.set_facecolor(SURFACE)
    ax.set_facecolor(SURFACE)

    group_counts = {}
    for sid in samples:
        counts = per_sample_counts[sid]
        N = depths_by_sample[sid]
        if N == 0:
            continue
        t = metadata.get(sid, {}).get("type_ofsample", "") or "Unknown"
        if t not in color_of:
            t = "Unknown"
        group_counts[t] = group_counts.get(t, 0) + 1

        xs, ys = [], []
        for n in depth_grid:
            if n > N:
                break
            xs.append(n)
            ys.append(expected_richness(counts, N, n))
        ax.plot(xs, ys, color=color_of[t], alpha=0.55, linewidth=0.9)

    if ref_depth:
        ax.axvline(ref_depth, color=BASELINE, linestyle="--", linewidth=1)
        ax.text(ref_depth, ax.get_ylim()[1] if ax.get_ylim()[1] > 0 else 1, "  step 5 rarefaction depth",
                rotation=90, va="top", ha="left", fontsize=7, color=INK_SECONDARY)

    ax.set_xlabel("Sequencing depth (reads)", color=INK_PRIMARY, fontsize=10)
    ax.set_ylabel("Expected ASV richness", color=INK_PRIMARY, fontsize=10)
    ax.set_title(f"{region} -- rarefaction curves ({len(samples)} samples)",
                color=INK_PRIMARY, fontsize=12)
    ax.xaxis.grid(True, color=GRIDLINE, linewidth=0.6)
    ax.yaxis.grid(True, color=GRIDLINE, linewidth=0.6)
    ax.set_axisbelow(True)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    for spine in ("left", "bottom"):
        ax.spines[spine].set_color(BASELINE)
    ax.tick_params(colors=INK_MUTED)

    from matplotlib.lines import Line2D
    handles = [Line2D([0], [0], color=color_of[g], lw=2, label=f"{g} (n={group_counts.get(g, 0)})")
              for g in TYPE_ORDER if group_counts.get(g, 0) > 0]
    ax.legend(handles=handles, loc="lower right", frameon=False, fontsize=8, labelcolor=INK_SECONDARY)

    fig.tight_layout()
    fig.savefig(out_path, facecolor=SURFACE)
    plt.close(fig)


def read_metadata(path):
    with open(path, newline="") as fh:
        return {row["sample-id"]: row for row in csv.DictReader(fh, delimiter="\t")}


if __name__ == "__main__":
    main()
