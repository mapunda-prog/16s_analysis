#!/usr/bin/env python3
"""
Sanity-check config/metadata.tsv against the actual per-region samplesheets
before spending hours on scripts/03_run_ampliseq.sh --extra "--metadata ...".

QIIME2's diversity plugins silently drop any sample missing from the metadata
file (or error, depending on the command) -- better to catch that here than
after a multi-hour ampliseq run.

Compares the union of sampleID across every samplesheets/<REGION>.tsv against
the sample-id column of config/metadata.tsv:

  - kept/runnable in >=1 region but missing from metadata -> real problem,
    fix before running. A sample excluded (low reads) in one region does NOT
    excuse it here if it's kept and running in another region.
  - in metadata but never sequenced/no sheet                -> stale row, warning only
  - excluded (low reads) everywhere it appears               -> informational only

Exit status is nonzero only for the first kind (an unexplained gap).

Usage:
    python3 scripts/check_metadata.py
    python3 scripts/check_metadata.py --samplesheet-dir samplesheets --metadata config/metadata.tsv
"""

import argparse
import csv
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
PROJECT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(HERE, "lib"))
from samplesheets import collect_sample_ids, is_control, read_excluded  # noqa: E402


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

    sequenced = collect_sample_ids(args.samplesheet_dir)  # kept/runnable only
    if not sequenced:
        sys.exit(f"no region samplesheets in {args.samplesheet_dir}")
    excluded = read_excluded(args.samplesheet_dir)

    with open(args.metadata) as fh:
        reader = csv.reader(fh, delimiter="\t")
        header = next(reader)
        col = header.index("sample-id")
        metadata_ids = {row[col] for row in reader if row}

    missing = sorted(set(sequenced) - metadata_ids)
    stale = sorted(metadata_ids - set(sequenced))
    ok = sorted(set(sequenced) & metadata_ids)
    partially_excluded = sorted(set(sequenced) & excluded)

    print(f"samplesheets: {len(sequenced)} distinct sample ID(s) kept in >=1 region")
    print(f"metadata.tsv: {len(metadata_ids)} sample ID(s)")
    print(f"  matched:                        {len(ok)}")
    print(f"  stale in metadata (no samplesheet, not fatal): {len(stale)}")
    if stale:
        print("    " + ", ".join(stale))
    print(f"  kept in one region, excluded (low reads) in another (informational): {len(partially_excluded)}")
    if partially_excluded:
        print("    " + ", ".join(partially_excluded))
    print(f"  MISSING from metadata (real problem):  {len(missing)}")
    if missing:
        print("    " + ", ".join(
            f"{sid} ({','.join(sorted(sequenced[sid]))})"
            + (" [control, but no blank row in metadata.tsv]" if is_control(sid) else "")
            for sid in missing
        ))
        sys.exit(1)


if __name__ == "__main__":
    main()
