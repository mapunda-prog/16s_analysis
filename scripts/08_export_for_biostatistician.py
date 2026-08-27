#!/usr/bin/env python3
"""
Step 08 -- per-sample organism-detection export for downstream biostatistics.

Reads each requested region's
results_by_region/<REGION>/qiime2/rel_abundance_tables/rel-table-ASV_with-DADA2-tax.tsv
(ampliseq's own QIIME2 output: per-ASV taxonomy + classifier confidence +
per-sample relative abundance, already free of Mitochondria/Chloroplast --
that filter is applied upstream of this table) and reshapes it from one
row per ASV / one column per sample into one row per sample x detected
organism, tidy-format:

    sample_id  region  is_control  Kingdom..Species_exact  organism
    ASV_ID  relative_abundance  confidence

"Detected" means relative_abundance > 0 in that sample x region. Sequencing
controls (NTC, PC1*, PC2*) are INCLUDED here (flagged via is_control), unlike
06_taxa_barplots.py which excludes them from the composition chart -- this
export is meant as a complete, unfiltered handoff for someone else's
analysis, where dropping data silently would be the wrong call.

Pure standard library.

Usage:
    python3 scripts/08_export_for_biostatistician.py --regions V3V4 V4V5
    python3 scripts/08_export_for_biostatistician.py --results-dir /path/to/results_by_region --regions V3V4 V4V5
"""

import argparse
import csv
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
PROJECT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(HERE, "lib"))
from samplesheets import is_control  # noqa: E402

TAX_RANKS = ["Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Species_exact"]


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--results-dir", default=os.path.join(PROJECT, "results_by_region"))
    ap.add_argument("--regions", nargs="+", required=True,
                    help="e.g. --regions V3V4 V4V5")
    ap.add_argument("--out", default=os.path.join(PROJECT, "organisms_by_sample", "organisms_by_sample.tsv"))
    args = ap.parse_args()

    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    out_cols = (["sample_id", "region", "is_control"] + TAX_RANKS
                + ["organism", "ASV_ID", "relative_abundance", "confidence"])

    n_rows = 0
    n_samples_seen = set()
    with open(args.out, "w", newline="") as out_fh:
        w = csv.writer(out_fh, delimiter="\t", lineterminator="\n")
        w.writerow(out_cols)

        for region in args.regions:
            path = os.path.join(args.results_dir, region, "qiime2", "rel_abundance_tables",
                                "rel-table-ASV_with-DADA2-tax.tsv")
            if not os.path.isfile(path):
                print(f"{region}: WARNING no rel-table-ASV_with-DADA2-tax.tsv at {path}, skipped")
                continue

            with open(path, newline="") as fh:
                reader = csv.reader(fh, delimiter="\t")
                header = [h.strip('"') for h in next(reader)]
                idx = {c: i for i, c in enumerate(header)}
                sample_cols = [c for c in header
                              if c not in ("ID", "confidence", "sequence") and c not in TAX_RANKS]

                region_rows = 0
                for row in reader:
                    row = [v.strip('"') for v in row]
                    asv_id = row[idx["ID"]]
                    tax = {r: row[idx[r]] for r in TAX_RANKS}
                    confidence = row[idx["confidence"]]
                    organism = organism_label(tax)
                    for sample in sample_cols:
                        abund = _float(row[idx[sample]])
                        if not abund:
                            continue
                        w.writerow([sample, region, is_control(sample)]
                                  + [tax[r] for r in TAX_RANKS]
                                  + [organism, asv_id, f"{abund:.8f}", confidence])
                        region_rows += 1
                        n_samples_seen.add(sample)
                n_rows += region_rows
                print(f"{region}: {len(sample_cols)} samples, {region_rows} detection rows")

    print(f"\nWrote {n_rows} rows, {len(n_samples_seen)} distinct samples -> {args.out}")


def organism_label(tax):
    """Most specific available name: 'Genus species_exact' > 'Genus species' >
    Genus > 'Family (family)' ... down to Kingdom, else 'Unclassified'."""
    genus = tax["Genus"]
    if genus and tax["Species_exact"]:
        return f"{genus} {tax['Species_exact']}"
    if genus and tax["Species"]:
        return f"{genus} {tax['Species']}"
    if genus:
        return genus
    for rank in ["Family", "Order", "Class", "Phylum", "Kingdom"]:
        if tax[rank]:
            return f"{tax[rank]} ({rank.lower()})"
    return "Unclassified"


def _float(v):
    try:
        return float(v)
    except (TypeError, ValueError):
        return 0.0


if __name__ == "__main__":
    main()
