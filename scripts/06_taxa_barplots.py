#!/usr/bin/env python3
"""
Step 06 -- per-region stacked taxonomic-composition barplots (genus level),
as a local substitute for QIIME2's own barplot when that process didn't run
or failed (see 05_diversity_analysis.py's docstring re: ITS/V5V7).

Reads compiled/genus_counts_by_region.tsv (already free of host mitochondria/
chloroplast -- see 04_compile_results.py) and config/metadata.tsv. Sequencing
controls (NTC, PC1*, PC2*) are excluded from the plot -- mixing spike-in
organisms into a composition chart of patient specimens is misleading, and
they're already covered in the main QC review.

The top 8 genera by total abundance across all requested regions get a fixed
color each (same genus = same color on every region's chart); everything
else folds into "Other". Samples are ordered by specimen type, then
collection site, so the Blood-vs-Swab compositional shift is visible at a
glance. Palette and mark choices follow the dataviz skill's categorical
formula (8-hue adjacent-pairlist order, validated CVD-safe for stacked bars).

Requires matplotlib (not part of the server pipeline's dependencies) -- run
this locally, like build_metadata.py / 05_diversity_analysis.py.

Usage:
    python3 scripts/06_taxa_barplots.py --regions V1V2 V2V3 V3V4 V4V5 V7V9
"""

import argparse
import csv
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

HERE = os.path.dirname(os.path.abspath(__file__))
PROJECT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(HERE, "lib"))
from samplesheets import is_control  # noqa: E402

N_TOP_GENERA = 8
TYPE_ORDER = ["Blood", "Swab", "Blood+Swab", "Unknown"]

# dataviz skill categorical palette, light mode, fixed adjacent-pairlist order
SERIES_COLORS = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100",
                 "#e87ba4", "#008300", "#4a3aa7", "#e34948"]
OTHER_COLOR = "#898781"      # skill's "muted" ink -- deliberately recessive
SURFACE = "#fcfcfb"
INK_PRIMARY = "#0b0b0b"
INK_SECONDARY = "#52514e"
INK_MUTED = "#898781"
GRIDLINE = "#e1e0d9"
BASELINE = "#c3c2b7"


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--genus-counts", default=os.path.join(PROJECT, "compiled", "genus_counts_by_region.tsv"))
    ap.add_argument("--metadata", default=os.path.join(PROJECT, "config", "metadata.tsv"))
    ap.add_argument("--outdir", default=os.path.join(PROJECT, "taxa_barplots"))
    ap.add_argument("--regions", nargs="+", required=True,
                    help="regions to plot, e.g. V1V2 V2V3 V3V4 V4V5 V7V9")
    args = ap.parse_args()

    metadata = read_metadata(args.metadata)
    counts = read_genus_counts(args.genus_counts, args.regions)  # region -> sample -> lineage -> count
    os.makedirs(args.outdir, exist_ok=True)

    top_genera = pick_top_genera(counts, metadata)
    color_of = {lineage: SERIES_COLORS[i] for i, lineage in enumerate(top_genera)}
    labels_of = {lineage: display_label(lineage) for lineage in top_genera}

    print(f"Top {N_TOP_GENERA} genera (fixed color across all regions):")
    for lineage in top_genera:
        print(f"  {color_of[lineage]}  {labels_of[lineage]}")
    print()

    for region in args.regions:
        if region not in counts:
            print(f"{region}: no data in {args.genus_counts}, skipped")
            continue
        samples = ordered_samples(counts[region], metadata)
        matrix = relative_abundance(counts[region], samples, top_genera)
        png_path = os.path.join(args.outdir, f"{region}_genus_barplot.png")
        tsv_path = os.path.join(args.outdir, f"{region}_genus_relative_abundance.tsv")
        plot_region(region, samples, matrix, top_genera, color_of, labels_of, metadata, png_path)
        write_matrix(samples, matrix, top_genera, labels_of, tsv_path)
        print(f"{region}: {len(samples)} samples -> {png_path}")


# --------------------------------------------------------------------------- #
# data
# --------------------------------------------------------------------------- #

def read_metadata(path):
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        return {row["sample-id"]: row for row in reader}


def read_genus_counts(path, wanted_regions):
    wanted = set(wanted_regions)
    counts = {}
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            region = row["region"]
            if region not in wanted or is_control(row["sampleID"]):
                continue
            counts.setdefault(region, {}).setdefault(row["sampleID"], {})[row["genus_lineage"]] = \
                int(row["count"])
    return counts


def pick_top_genera(counts, metadata):
    totals = {}
    for region_counts in counts.values():
        for sample, lineages in region_counts.items():
            for lineage, c in lineages.items():
                totals[lineage] = totals.get(lineage, 0) + c
    ranked = sorted(totals, key=lambda k: -totals[k])
    return ranked[:N_TOP_GENERA]


def ordered_samples(region_counts, metadata):
    def key(sid):
        meta = metadata.get(sid, {})
        t = meta.get("type_ofsample", "") or "Unknown"
        site = meta.get("hf_name", "") or "Unknown"
        type_rank = TYPE_ORDER.index(t) if t in TYPE_ORDER else len(TYPE_ORDER)
        return (type_rank, site, sid)
    return sorted(region_counts, key=key)


def relative_abundance(region_counts, samples, top_genera):
    """-> dict[lineage or 'Other'] = [proportion per sample, in `samples` order]."""
    matrix = {lineage: [] for lineage in top_genera}
    matrix["Other"] = []
    for sid in samples:
        lineages = region_counts.get(sid, {})
        total = sum(lineages.values())
        for lineage in top_genera:
            c = lineages.get(lineage, 0)
            matrix[lineage].append(c / total if total else 0.0)
        other = sum(c for lin, c in lineages.items() if lin not in top_genera)
        matrix["Other"].append(other / total if total else 0.0)
    return matrix


def display_label(lineage):
    parts = [p for p in lineage.split(";") if p]
    return parts[-1] if parts else "Unassigned"


# --------------------------------------------------------------------------- #
# plotting
# --------------------------------------------------------------------------- #

def plot_region(region, samples, matrix, top_genera, color_of, labels_of, metadata, out_path):
    n = len(samples)
    fig_w = max(10, n * 0.22)
    fig, ax = plt.subplots(figsize=(fig_w, 6), dpi=150)
    fig.patch.set_facecolor(SURFACE)
    ax.set_facecolor(SURFACE)

    x = range(n)
    bottom = [0.0] * n
    for lineage in top_genera + ["Other"]:
        color = color_of.get(lineage, OTHER_COLOR)
        ax.bar(x, matrix[lineage], bottom=bottom, width=0.85, color=color,
               edgecolor=SURFACE, linewidth=0.6)
        bottom = [b + v for b, v in zip(bottom, matrix[lineage])]

    # group dividers + labels (specimen type, in the sample ordering already applied)
    groups = []
    for i, sid in enumerate(samples):
        meta = metadata.get(sid, {})
        t = meta.get("type_ofsample", "") or "Unknown"
        if not groups or groups[-1][0] != t:
            groups.append([t, i, i])
        else:
            groups[-1][2] = i
    for t, start, end in groups:
        mid = (start + end) / 2
        ax.text(mid, 1.03, t, ha="center", va="bottom", fontsize=9, color=INK_SECONDARY)
        if start > 0:
            ax.axvline(start - 0.5, color=BASELINE, linewidth=0.8, linestyle="--")

    ax.set_xlim(-0.6, n - 0.4)
    ax.set_ylim(0, 1)
    ax.set_xticks(list(x))
    ax.set_xticklabels(samples, rotation=90, fontsize=5, color=INK_SECONDARY)
    ax.set_ylabel("Relative abundance", color=INK_PRIMARY, fontsize=10)
    ax.set_title(f"{region} -- genus-level composition (top {N_TOP_GENERA} + Other)",
                color=INK_PRIMARY, fontsize=12, pad=20)

    ax.yaxis.grid(True, color=GRIDLINE, linewidth=0.6)
    ax.set_axisbelow(True)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    for spine in ("left", "bottom"):
        ax.spines[spine].set_color(BASELINE)
    ax.tick_params(colors=INK_MUTED)

    handles = [Patch(facecolor=color_of[l], label=labels_of[l]) for l in top_genera]
    handles.append(Patch(facecolor=OTHER_COLOR, label="Other"))
    ax.legend(handles=handles, bbox_to_anchor=(1.01, 1), loc="upper left",
             frameon=False, fontsize=8, labelcolor=INK_SECONDARY)

    fig.tight_layout()
    fig.savefig(out_path, facecolor=SURFACE)
    plt.close(fig)


def write_matrix(samples, matrix, top_genera, labels_of, path):
    cols = top_genera + ["Other"]
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["sample-id"] + [labels_of.get(c, c) for c in cols])
        for i, sid in enumerate(samples):
            w.writerow([sid] + [f"{matrix[c][i]:.6f}" for c in cols])


if __name__ == "__main__":
    main()
