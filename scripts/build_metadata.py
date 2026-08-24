#!/usr/bin/env python3
"""
Build config/metadata.tsv (QIIME2 format) from the clinical metadata Stata
files in metadata/.

Two files land in metadata/:
    "Meta data for AFI_submit.dta"  -- blood-draw febrile-illness panel
    "Meta data for ARI_submit.dta"  -- swab respiratory-infection panel

Both are keyed by enrolment_id_or, e.g. "MLTP-AIC-125", with no BL/SW suffix --
that only appears once a patient's fastq is demultiplexed and the sample ID
gets a specimen-type suffix (AIC125BL, AIC125SW) if both specimens were drawn.

This script maps each *sequenced* sample ID (as it appears in
results*/overall_summary.tsv, or a samplesheet you point it at) back to its
enrolment code, picks the matching panel by specimen-type suffix (BL -> AFI,
SW -> ARI), and writes one combined row per sample to config/metadata.tsv.

Ambiguous cases (an enrolment has both panels but the sequenced sample ID has
no BL/SW suffix to tell them apart -- duplicate/re-run IDs like AIC057a/057b
or AIC115_S5/_S12) are NOT guessed at: both panels are merged into that row
and the sample ID is listed under "AMBIGUOUS" in the run report so a person
can confirm which specimen it actually was.

Requires pandas (not part of the server pipeline's dependencies) -- run this
locally to produce config/metadata.tsv, then commit that file; nothing else
in the pipeline needs pandas.

Usage:
    python3 scripts/build_metadata.py
    python3 scripts/build_metadata.py --sample-ids-file samplesheets/V3V4.tsv
    python3 scripts/build_metadata.py --afi metadata/AFI.dta --ari metadata/ARI.dta
"""

import argparse
import csv
import os
import re
import sys

import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
PROJECT = os.path.dirname(HERE)

# sample IDs that are sequencing controls, not patient specimens -- expected
# to have no clinical metadata, so they're not reported as unmatched.
CONTROL_IDS = {"NTC", "PC1NPHL", "PC2EXT", "Undetermined"}

SHARED_COLUMNS = [
    "recruitment_month", "sex", "hf_name", "location", "agegroup", "seasonal",
]

AFI_PATHOGEN_COLUMNS = [
    "brucella_pathogen", "coxiella_pathogen", "denv_pathogen", "chikv_pathogen",
    "zikv_pathogen", "wnv_pathogen", "uniplex_pcr_brucella",
]

ARI_PATHOGEN_COLUMNS = [
    "rv_general_result", "rsv_pathogen", "sars_cov_2_pathogen", "flu_a_pathogen",
    "flu_b_pathogen", "mpv_pathogen", "adv_pathogen", "hrv_pathogen",
    "piv_pathogen", "seegene_general_result", "hbov_pathogen",
    "coronavirus_229e_pathogen", "coronavirus_nl63_pathogen",
    "coronavirus_oc43_pathogen", "panels3_hrv_pathogen", "bpp_pathogen",
    "bp_pathogen", "cp_pathogen", "hi_pathogen", "lp_pathogen", "mp_pathogen",
    "sp_pathogen", "pathogen", "virus", "bacteria", "type_patho",
]


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--afi", default=os.path.join(PROJECT, "metadata", "Meta data for AFI_submit.dta"))
    ap.add_argument("--ari", default=os.path.join(PROJECT, "metadata", "Meta data for ARI_submit.dta"))
    ap.add_argument("--sample-ids-file", default=None,
                    help="TSV with a sampleID/sample column to convert "
                         "[default: results*/overall_summary.tsv, sample column]")
    ap.add_argument("--out", default=os.path.join(PROJECT, "config", "metadata.tsv"))
    args = ap.parse_args()

    sample_ids = read_sample_ids(args.sample_ids_file)
    afi = load_panel(args.afi)
    ari = load_panel(args.ari)
    site_codes = {site_of(e) for e in list(afi) + list(ari)}

    rows = []
    ambiguous, no_metadata, controls = [], [], []

    for sid in sample_ids:
        parsed = parse_sample_id(sid, site_codes)
        if parsed is None:
            (controls if sid in CONTROL_IDS else no_metadata).append(sid)
            rows.append(blank_row(sid))
            continue

        enrolment_id_or, specimen = parsed
        afi_row = afi.get(enrolment_id_or)
        ari_row = ari.get(enrolment_id_or)

        if specimen == "BL":
            chosen = [("AFI", afi_row)] if afi_row is not None else []
        elif specimen == "SW":
            chosen = [("ARI", ari_row)] if ari_row is not None else []
        else:
            chosen = [p for p in [("AFI", afi_row), ("ARI", ari_row)] if p[1] is not None]

        if not chosen:
            no_metadata.append(sid)
            rows.append(blank_row(sid))
            continue
        if len(chosen) > 1:
            ambiguous.append(sid)

        rows.append(build_row(sid, enrolment_id_or, chosen))

    write_metadata(rows, args.out)

    print(f"wrote {len(rows)} rows -> {args.out}")
    print(f"  matched (unambiguous): {len(rows) - len(ambiguous) - len(no_metadata) - len(controls)}")
    print(f"  ambiguous (both AFI+ARI, no BL/SW suffix -- merged, confirm specimen type): {len(ambiguous)}")
    if ambiguous:
        print("    " + ", ".join(ambiguous))
    print(f"  no clinical metadata found: {len(no_metadata)}")
    if no_metadata:
        print("    " + ", ".join(no_metadata))
    print(f"  controls (expected, no metadata): {len(controls)}")


def load_panel(path):
    """dta -> dict[enrolment_id_or] = row dict of strings, NaN -> ''."""
    df = pd.read_stata(path)
    out = {}
    for _, row in df.iterrows():
        d = {}
        for col, val in row.items():
            d[col] = "" if pd.isna(val) else str(val).strip()
        out[d["enrolment_id_or"]] = d
    return out


def site_of(enrolment_id_or):
    return enrolment_id_or.split("-")[1]


_SID_RE = None


def parse_sample_id(sid, site_codes):
    """'AIC125BL' -> ('MLTP-AIC-125', 'BL'); 'AIC057a' -> (..., None); None if
    sid isn't a <site><digits><suffix> pattern with a known site code."""
    global _SID_RE
    if _SID_RE is None:
        _SID_RE = re.compile(r"^([A-Za-z]+)(\d+)(.*)$")
    m = _SID_RE.match(sid)
    if not m:
        return None
    site, num, suffix = m.group(1), m.group(2), m.group(3)
    if site not in site_codes:
        return None
    specimen = suffix if suffix in ("BL", "SW") else None
    return f"MLTP-{site}-{num}", specimen


def blank_row(sid):
    return {"sample-id": sid}


def build_row(sid, enrolment_id_or, chosen):
    row = {"sample-id": sid, "enrolment_id_or": enrolment_id_or}
    sources = [name for name, _ in chosen]
    row["panel_source"] = "+".join(sources)

    first = chosen[0][1]
    row["screening_id"] = first.get("screening_id", "")
    for col in SHARED_COLUMNS:
        row[col] = first.get(col, "")

    type_vals = sorted({r.get("type_ofsample", "") for _, r in chosen} - {""})
    row["type_ofsample"] = "+".join(type_vals)

    for name, r in chosen:
        cols = AFI_PATHOGEN_COLUMNS if name == "AFI" else ARI_PATHOGEN_COLUMNS
        prefix = "afi_" if name == "AFI" else "ari_"
        for col in cols:
            row[prefix + col] = r.get(col, "")
    return row


def write_metadata(rows, out_path):
    header = ["sample-id", "enrolment_id_or", "screening_id", "panel_source", "type_ofsample"]
    header += SHARED_COLUMNS
    header += ["afi_" + c for c in AFI_PATHOGEN_COLUMNS]
    header += ["ari_" + c for c in ARI_PATHOGEN_COLUMNS]

    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    with open(out_path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=header, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for row in rows:
            w.writerow({k: row.get(k, "") for k in header})


def read_sample_ids(path):
    if path is None:
        path = find_default_summary()
    col = None
    with open(path) as fh:
        reader = csv.reader(fh, delimiter="\t")
        header = next(reader)
        for candidate in ("sample", "sampleID", "sample-id"):
            if candidate in header:
                col = header.index(candidate)
                break
        if col is None:
            sys.exit(f"{path}: no sample/sampleID/sample-id column in header {header}")
        return [r[col] for r in reader if r]


def find_default_summary():
    for candidate in (
        os.path.join(PROJECT, "results", "overall_summary.tsv"),
        os.path.join(os.path.dirname(PROJECT), "results", "overall_summary.tsv"),
    ):
        if os.path.isfile(candidate):
            return candidate
    sys.exit("no results*/overall_summary.tsv found; pass --sample-ids-file")


if __name__ == "__main__":
    main()
