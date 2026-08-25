#!/usr/bin/env python3
"""Shared helper: collect sample IDs from the per-region ampliseq samplesheets
written by scripts/02_make_samplesheets.py, used by build_metadata.py and
check_metadata.py so both agree on what "the sequenced samples" means."""

import csv
import glob
import os


def collect_sample_ids(samplesheet_dir):
    """dict[sample_id] -> set(regions) it's kept/runnable in. A sample can be
    kept in one region and excluded (low reads) in another -- both states
    coexist, so callers must not treat "excluded somewhere" as "excluded
    everywhere"."""
    sequenced = {}
    for sheet in sorted(glob.glob(os.path.join(samplesheet_dir, "*.tsv"))):
        if os.path.basename(sheet) == "excluded_samples.tsv":
            continue
        region = os.path.splitext(os.path.basename(sheet))[0]
        with open(sheet) as fh:
            reader = csv.reader(fh, delimiter="\t")
            header = next(reader)
            col = header.index("sampleID") if "sampleID" in header else header.index("sample")
            for row in reader:
                if row:
                    sequenced.setdefault(row[col], set()).add(region)
    return sequenced


def read_excluded(samplesheet_dir):
    """set of sample IDs excluded (low reads) in at least one region."""
    path = os.path.join(samplesheet_dir, "excluded_samples.tsv")
    if not os.path.isfile(path):
        return set()
    with open(path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        return {row["sampleID"] for row in reader}


def is_control(sample_id):
    """Sequencing controls (NTC, positive controls, undetermined reads) are
    expected to have no clinical metadata -- naming varies per run (e.g.
    PC1NPHL vs PC1Ecoli), so match by prefix rather than a fixed ID list."""
    return sample_id in ("NTC", "Undetermined") or sample_id.startswith(("PC1", "PC2"))
