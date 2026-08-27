#!/usr/bin/env python3
"""
Step 09 -- left-join the per-sample organism export (step 8) with the AFI/ARI
clinical metadata (config/metadata.tsv, from build_metadata.py), so each
detected-organism row also carries its sample's demographics and pathogen
panel results.

Every column from config/metadata.tsv is carried over except sample-id
itself (organisms_by_sample.tsv already has sample_id) -- this stays in sync
automatically if build_metadata.py's schema changes, rather than hardcoding
the AFI/ARI column list here too.

A sample in organisms_by_sample.tsv with no row in metadata.tsv gets blank
metadata columns and is reported at the end; per check_metadata.py, this
should not happen for a metadata.tsv that's been kept up to date, but the
join does not silently drop the organism rows either way.

Pure standard library.

Usage:
    python3 scripts/09_merge_organisms_with_metadata.py
    python3 scripts/09_merge_organisms_with_metadata.py --organisms organisms_by_sample/organisms_by_sample.tsv --metadata config/metadata.tsv
"""

import argparse
import csv
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
PROJECT = os.path.dirname(HERE)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--organisms", default=os.path.join(PROJECT, "organisms_by_sample", "organisms_by_sample.tsv"))
    ap.add_argument("--metadata", default=os.path.join(PROJECT, "config", "metadata.tsv"))
    ap.add_argument("--out", default=os.path.join(PROJECT, "organisms_by_sample", "organisms_by_sample_with_metadata.tsv"))
    args = ap.parse_args()

    if not os.path.isfile(args.organisms):
        sys.exit(f"not found: {args.organisms}\nRun scripts/08_export_for_biostatistician.py first.")
    if not os.path.isfile(args.metadata):
        sys.exit(f"not found: {args.metadata}\nRun scripts/build_metadata.py first.")

    meta_by_id, meta_cols = read_metadata(args.metadata)

    unmatched = set()
    matched_samples = set()
    n_rows = 0
    with open(args.organisms, newline="") as in_fh, open(args.out, "w", newline="") as out_fh:
        reader = csv.DictReader(in_fh, delimiter="\t")
        out_cols = reader.fieldnames + meta_cols
        w = csv.DictWriter(out_fh, fieldnames=out_cols, delimiter="\t", lineterminator="\n")
        w.writeheader()

        for row in reader:
            sid = row["sample_id"]
            meta = meta_by_id.get(sid)
            if meta is None:
                unmatched.add(sid)
                meta = {c: "" for c in meta_cols}
            else:
                matched_samples.add(sid)
            w.writerow({**row, **meta})
            n_rows += 1

    print(f"Wrote {n_rows} rows, {len(matched_samples)} samples matched to metadata -> {args.out}")
    if unmatched:
        print(f"\nWARNING: {len(unmatched)} sample(s) in {args.organisms} have no row in "
              f"{args.metadata} -- their metadata columns are blank:")
        print("  " + ", ".join(sorted(unmatched)))


def read_metadata(path):
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        cols = [c for c in reader.fieldnames if c != "sample-id"]
        by_id = {}
        for row in reader:
            by_id[row["sample-id"]] = {c: row.get(c, "") for c in cols}
        return by_id, cols


if __name__ == "__main__":
    main()
