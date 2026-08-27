#!/usr/bin/env python3
"""
Step 10 -- the reverse of check_metadata.py: samples that have real clinical
metadata (a genuine AFI or ARI match, panel_source non-empty in
config/metadata.tsv -- so not a blank control row, not one of the "no
clinical metadata found" gaps) but no detection rows at all in the organism
export (step 8) for the requested regions.

For each such sample and region, explains why using the same samplesheets
step 8 was built from:
  - "excluded (low reads): N read pairs, threshold M"   -- never ran ampliseq
  - "kept, 0 ASVs survived (DADA2_input=X, merged=Y, nonchim=Z)"  -- ran, but
     nothing came out (usually a merge/truncation failure for that sample)
  - "not sequenced in this region"                       -- not in the
     samplesheet or excluded_samples.tsv either; worth checking demux directly

Pure standard library.

Usage:
    python3 scripts/10_samples_missing_results.py --regions V3V4 V4V5
"""

import argparse
import csv
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
PROJECT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(HERE, "lib"))
from samplesheets import collect_sample_ids, read_excluded  # noqa: E402


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--metadata", default=os.path.join(PROJECT, "config", "metadata.tsv"))
    ap.add_argument("--organisms", default=os.path.join(PROJECT, "organisms_by_sample", "organisms_by_sample.tsv"))
    ap.add_argument("--results-dir", default=os.path.join(PROJECT, "results_by_region"))
    ap.add_argument("--samplesheet-dir", default=os.path.join(PROJECT, "samplesheets"))
    ap.add_argument("--regions", nargs="+", required=True, help="e.g. --regions V3V4 V4V5")
    args = ap.parse_args()

    meta_ids = read_metadata_ids_with_real_data(args.metadata)
    result_ids = read_result_sample_ids(args.organisms)
    missing = sorted(meta_ids - result_ids)

    print(f"{len(meta_ids)} samples have real clinical metadata; "
          f"{len(result_ids)} have >=1 detection row in {args.regions}")
    print(f"{len(missing)} sample(s) have metadata but no results:\n")

    if not missing:
        return

    sequenced = collect_sample_ids(args.samplesheet_dir)   # sample -> set(regions kept)
    excluded = read_excluded_detail(args.samplesheet_dir)  # (region, sample) -> (n, threshold)
    summaries = {r: read_overall_summary(os.path.join(args.results_dir, r, "overall_summary.tsv"))
                for r in args.regions}

    for sid in missing:
        print(sid)
        for region in args.regions:
            print(f"  {region}: {explain(sid, region, sequenced, excluded, summaries[region])}")


def explain(sid, region, sequenced, excluded, summary):
    if (region, sid) in excluded:
        n, threshold = excluded[(region, sid)]
        return f"excluded (low reads): {n} read pairs, threshold {threshold}"
    if region in sequenced.get(sid, set()):
        row = summary.get(sid)
        if row is None:
            return "kept in samplesheet, but not in overall_summary.tsv -- check the ampliseq run"
        return (f"kept, 0 ASVs survived (DADA2_input={row.get('DADA2_input','?')}, "
                f"merged={row.get('merged','?')}, nonchim={row.get('nonchim','?')})")
    return "not sequenced in this region (not in samplesheet or excluded_samples.tsv)"


def read_metadata_ids_with_real_data(path):
    with open(path, newline="") as fh:
        return {row["sample-id"] for row in csv.DictReader(fh, delimiter="\t")
                if row.get("panel_source")}


def read_result_sample_ids(path):
    with open(path, newline="") as fh:
        return {row["sample_id"] for row in csv.DictReader(fh, delimiter="\t")}


def read_excluded_detail(samplesheet_dir):
    path = os.path.join(samplesheet_dir, "excluded_samples.tsv")
    out = {}
    if not os.path.isfile(path):
        return out
    with open(path, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            out[(row["region"], row["sampleID"])] = (row["read_pairs"], row["min_reads_threshold"])
    return out


def read_overall_summary(path):
    if not os.path.isfile(path):
        return {}
    with open(path, newline="") as fh:
        return {row["sample"]: row for row in csv.DictReader(fh, delimiter="\t")}


if __name__ == "__main__":
    main()
