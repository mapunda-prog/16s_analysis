#!/usr/bin/env python3
"""
Step 11 -- per-sample organism detections collapsed to a fixed taxonomic
rank (phylum, family, genus, or species), across multiple regions
including ITS.

Unlike 08_export_for_biostatistician.py (ASV-level, V3V4/V4V5 only, built
from ampliseq's merged QIIME2 rel-abundance-with-taxonomy table), this reads
each region's own dada2/ASV_table.tsv + dada2/ASV_tax_species.<ref>.tsv
directly -- the same layout regardless of reference database, so it works
uniformly for SILVA regions (V3V4, V4V5) and UNITE-fungi (ITS), which has no
merged QIIME2 table of its own. Mitochondria/Chloroplast ASVs are excluded
before relative abundance is computed, same as 04_compile_results.py.

For each requested rank, ASVs sharing the same lineage *up to that rank*
(deeper ranks ignored, same semantics as QIIME2's own `taxa collapse`) are
merged into one row per sample x region:
  - relative_abundance: summed across the merged ASVs
  - confidence: abundance-weighted mean of the merged ASVs' own classifier
    confidence (there is no such thing as "the confidence of a genus" --
    only per-ASV classification calls have one, so this is a derived
    summary, not a native score -- documented here and in the output README)
  - n_asvs_collapsed: how many ASVs went into this row, for transparency

Writes one file per rank into --outdir: phylum_level.tsv, family_level.tsv,
genus_level.tsv, species_level.tsv. Controls (NTC, PC1*, PC2*) are INCLUDED (flagged via
is_control), matching 08_export_for_biostatistician.py's convention for a
complete handoff.

Pure standard library.

Usage:
    python3 scripts/11_export_by_rank.py
    python3 scripts/11_export_by_rank.py --regions V3V4 V4V5 ITS --ranks phylum family genus species
"""

import argparse
import csv
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
PROJECT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(HERE, "lib"))
from samplesheets import is_control  # noqa: E402
from taxonomy import load_region_asvs  # noqa: E402

RANK_COLUMNS = {
    "phylum": ["Kingdom", "Phylum"],
    "family": ["Kingdom", "Phylum", "Class", "Order", "Family"],
    "genus": ["Kingdom", "Phylum", "Class", "Order", "Family", "Genus"],
    "species": ["Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Species_exact"],
}


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--results-dir", default=os.path.join(PROJECT, "results_by_region"))
    ap.add_argument("--regions", nargs="+", default=["V3V4", "V4V5", "ITS"])
    ap.add_argument("--ranks", nargs="+", default=["phylum", "family", "genus", "species"],
                    choices=list(RANK_COLUMNS))
    ap.add_argument("--outdir", default=os.path.join(PROJECT, "organisms_by_sample"))
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    # load once per region, reuse across ranks
    region_data = {}
    for region in args.regions:
        region_dir = os.path.join(args.results_dir, region)
        try:
            samples, asvs = load_region_asvs(region_dir)
        except FileNotFoundError as e:
            print(f"{region}: WARNING {e}, skipped")
            continue
        totals = per_sample_totals(samples, asvs)
        region_data[region] = (samples, asvs, totals)
        print(f"{region}: {len(samples)} samples, {len(asvs)} non-contaminant ASVs")
    print()

    for rank in args.ranks:
        out_path = os.path.join(args.outdir, f"{rank}_level.tsv")
        n_rows = write_rank_level(rank, args.regions, region_data, out_path)
        print(f"{rank}: {n_rows} rows -> {out_path}")


def per_sample_totals(samples, asvs):
    totals = {s: 0 for s in samples}
    for asv in asvs:
        for s, c in asv["counts"].items():
            totals[s] += c
    return totals


def write_rank_level(rank, regions, region_data, out_path):
    cols = RANK_COLUMNS[rank]
    out_cols = (["sample_id", "region", "is_control"] + cols
                + ["organism", "relative_abundance", "confidence", "n_asvs_collapsed"])
    n_rows = 0
    with open(out_path, "w", newline="") as out_fh:
        w = csv.writer(out_fh, delimiter="\t", lineterminator="\n")
        w.writerow(out_cols)

        for region in regions:
            if region not in region_data:
                continue
            samples, asvs, totals = region_data[region]
            # (sample, lineage_key) -> {abund_sum, weighted_conf_sum, n_asvs}
            groups = {}
            for asv in asvs:
                key = tuple(asv["tax"].get(c, "") for c in cols)
                conf = _float(asv["confidence"])
                for s, c in asv["counts"].items():
                    if c <= 0 or totals[s] == 0:
                        continue
                    rel = c / totals[s]
                    g = groups.setdefault((s, key), {"abund": 0.0, "wconf": 0.0, "n": 0})
                    g["abund"] += rel
                    g["wconf"] += rel * conf
                    g["n"] += 1

            for (sample, key), g in sorted(groups.items()):
                organism = organism_label(dict(zip(cols, key)), cols)
                confidence = g["wconf"] / g["abund"] if g["abund"] else 0.0
                w.writerow([sample, region, is_control(sample)] + list(key)
                          + [organism, f"{g['abund']:.8f}", f"{confidence:.6f}", g["n"]])
                n_rows += 1
    return n_rows


def organism_label(tax, cols):
    """Most specific available name given this rank's column set: 'Genus
    species_exact' > 'Genus species' > Genus > 'Family (family)' ... else
    'Unclassified'. Mirrors 08_export_for_biostatistician.py's convention."""
    genus = tax.get("Genus", "")
    if "Species_exact" in cols and tax.get("Species_exact"):
        return f"{genus} {tax['Species_exact']}".strip()
    if "Species" in cols and tax.get("Species"):
        return f"{genus} {tax['Species']}".strip()
    if genus:
        return genus
    for c in reversed(cols):
        if c in ("Genus", "Species", "Species_exact"):
            continue
        if tax.get(c):
            return f"{tax[c]} ({c.lower()})"
    return "Unclassified"


def _float(v):
    try:
        return float(v)
    except (TypeError, ValueError):
        return 0.0


if __name__ == "__main__":
    main()
