#!/usr/bin/env python3
"""
Sanity-check config/metadata.tsv against the actual per-region samplesheets
before spending hours on scripts/03_run_ampliseq.sh --extra "--metadata ...".

QIIME2's diversity plugins silently drop any sample missing from the metadata
file (or error, depending on the command) -- better to catch that here than
after a multi-hour ampliseq run.

Compares the union of sampleID across every samplesheets/<REGION>.tsv against
the sample-id column of config/metadata.tsv:

  - in samplesheets but missing from metadata  -> real problem, fix before running
  - in samplesheets/excluded_samples.tsv        -> expected absence, not flagged
  - in metadata but never sequenced/no sheet    -> stale row, warning only

Exit status is nonzero only for the first kind (an unexplained gap).

Usage:
    python3 scripts/check_metadata.py
    python3 scripts/check_metadata.py --samplesheet-dir samplesheets --metadata config/metadata.tsv
"""

import argparse
import csv
import glob
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
PROJECT = os.path.dirname(HERE)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--samplesheet-dir", default=os.path.join(PROJECT, "samplesheets"))
    ap.add_argument("--metadata", default=os.path.join(PROJECT, "config", "metadata.tsv"))
    args = ap.parse_args()

    if not os.path.isdir(args.samplesheet_dir):
        sys.exit(f"samplesheet dir not found: {args.samplesheet_dir}\n"
                  f"Run scripts/02_make_samplesheets.py first.")
    if not os.path.isfile(args.metadata):
        sys.exit(f"metadata file not found: {args.metadata}\n"
                  f"Run scripts/build_metadata.py first.")

    sequenced = {}  # sample_id -> set of regions it appears in
    sheets = sorted(glob.glob(os.path.join(args.samplesheet_dir, "*.tsv")))
    sheets = [s for s in sheets if os.path.basename(s) != "excluded_samples.tsv"]
    if not sheets:
        sys.exit(f"no region samplesheets in {args.samplesheet_dir}")
    for sheet in sheets:
        region = os.path.splitext(os.path.basename(sheet))[0]
        with open(sheet) as fh:
            reader = csv.reader(fh, delimiter="\t")
            header = next(reader)
            col = header.index("sampleID") if "sampleID" in header else header.index("sample")
            for row in reader:
                if row:
                    sequenced.setdefault(row[col], set()).add(region)

    excluded = read_excluded(os.path.join(args.samplesheet_dir, "excluded_samples.tsv"))

    with open(args.metadata) as fh:
        reader = csv.reader(fh, delimiter="\t")
        header = next(reader)
        col = header.index("sample-id")
        metadata_ids = {row[col] for row in reader if row}

    missing = sorted(set(sequenced) - metadata_ids - excluded)
    stale = sorted(metadata_ids - set(sequenced))
    ok = sorted(set(sequenced) & metadata_ids)

    print(f"samplesheets: {len(sheets)} region(s), {len(sequenced)} distinct sample ID(s)")
    print(f"metadata.tsv: {len(metadata_ids)} sample ID(s)")
    print(f"  matched:                    {len(ok)}")
    print(f"  excluded (expected, low-read, not flagged): {len(excluded & set(sequenced))}")
    print(f"  stale in metadata (no samplesheet, not fatal): {len(stale)}")
    if stale:
        print("    " + ", ".join(stale))
    print(f"  MISSING from metadata (real problem):  {len(missing)}")
    if missing:
        print("    " + ", ".join(f"{sid} ({','.join(sorted(sequenced[sid]))})" for sid in missing))
        sys.exit(1)


def read_excluded(path):
    if not os.path.isfile(path):
        return set()
    with open(path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        return {row["sampleID"] for row in reader}


if __name__ == "__main__":
    main()
