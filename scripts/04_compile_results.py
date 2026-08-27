#!/usr/bin/env python3
"""
Step 04 -- compile the per-region ampliseq results into cross-region tables.

Reads results_by_region/<REGION>/ and writes to compiled/:

  read_tracking_by_region.tsv   per sample per region, every column of that
                                region's overall_summary.tsv, region prepended
  region_qc_summary.tsv         one row per region: samples, reads in/out,
                                median % retained, ASV and genus counts
  asv_long.tsv.gz               region, sample, ASV_ID, count (nonzero only),
                                full taxonomy -- the tidy table to load in R
  genus_counts_by_region.tsv    region, sample, genus lineage, summed count
  genus_matrix_across_regions.tsv  genus lineage x region totals, plus how many
                                regions detected it -- the cross-region view
  demux_summary.tsv             copy of the step 01 assignment counts, if present

ASVs classified as Mitochondria or Chloroplast (host/plant DNA, not target
bacteria) are dropped from every output here -- ASV_table.tsv is upstream of
ampliseq's own --exclude_taxa mitochondria,chloroplast filter, which only
applies to the QIIME2 abundance tables, so without this these compiled files
disagree with the pipeline's own filtered results. This matters most for
blood-derived specimens, where host mitochondrial DNA can dominate a sample's
reads entirely.

Pure standard library: no pandas needed on the analysis server.

Usage:
    python3 scripts/04_compile_results.py
    python3 scripts/04_compile_results.py --results-dir results_by_region --outdir compiled
"""

import argparse
import csv
import glob
import gzip
import os
import statistics
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "lib"))
from taxonomy import is_host_contaminant  # noqa: E402

TAX_RANKS = ["Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"]


def read_tsv(path):
    """-> (header list, list of row lists). Skips biom-style '# Constructed' lines."""
    with open(path, newline="") as fh:
        rows = [r for r in csv.reader(fh, delimiter="\t")
                if r and not (len(r) == 1 and r[0].startswith("# "))]
    if not rows:
        return [], []
    return rows[0], rows[1:]


def find_one(pattern):
    hits = sorted(glob.glob(pattern))
    return hits[0] if hits else None


def main():
    here = os.path.dirname(os.path.abspath(__file__))
    project = os.path.dirname(here)

    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--results-dir", default=os.path.join(project, "results_by_region"))
    ap.add_argument("--outdir", default=os.path.join(project, "compiled"))
    ap.add_argument("--demux-counts", default=os.path.join(project, "demux", "demux_counts.tsv"))
    ap.add_argument("--regions", nargs="*",
                    help="only compile these regions (default: every subdirectory "
                         "of --results-dir that looks like an ampliseq run)")
    args = ap.parse_args()

    rroot = os.path.abspath(args.results_dir)
    if not os.path.isdir(rroot):
        sys.exit(f"results dir not found: {rroot}\nRun scripts/03_run_ampliseq.sh first.")
    os.makedirs(args.outdir, exist_ok=True)

    regions = args.regions or sorted(
        d for d in os.listdir(rroot)
        if os.path.isdir(os.path.join(rroot, d))
        and (os.path.isdir(os.path.join(rroot, d, "dada2"))
             or os.path.isfile(os.path.join(rroot, d, "overall_summary.tsv")))
    )
    if not regions:
        sys.exit(f"no completed region results found under {rroot}")
    print(f"Compiling {len(regions)} region(s): {', '.join(regions)}\n")

    tracking_rows, tracking_cols = [], []
    qc_rows = []
    asv_long_path = os.path.join(args.outdir, "asv_long.tsv.gz")
    genus_rows = []          # (region, sample, lineage, count)
    skipped = []

    with gzip.open(asv_long_path, "wt", newline="") as asv_fh:
        asv_w = csv.writer(asv_fh, delimiter="\t", lineterminator="\n")
        asv_w.writerow(["region", "sampleID", "ASV_ID", "count"] + TAX_RANKS)

        for region in regions:
            rdir = os.path.join(rroot, region)

            # ---------------- read tracking ----------------
            summ = os.path.join(rdir, "overall_summary.tsv")
            n_samples = 0
            retained = []
            total_in = total_out = 0
            if os.path.isfile(summ):
                hdr, rows = read_tsv(summ)
                for c in hdr:
                    if c not in tracking_cols:
                        tracking_cols.append(c)
                idx = {c: i for i, c in enumerate(hdr)}
                for row in rows:
                    d = {c: (row[i] if i < len(row) else "") for c, i in idx.items()}
                    tracking_rows.append((region, d))
                    n_samples += 1
                    total_in += _int(d.get("DADA2_input"))
                    total_out += _int(d.get("nonchim"))
                    pct = _float(d.get("retained_percent"))
                    if pct is not None:
                        retained.append(pct)
            else:
                print(f"  {region}: WARNING no overall_summary.tsv")

            # ---------------- ASV table + taxonomy ----------------
            asv_table = os.path.join(rdir, "dada2", "ASV_table.tsv")
            tax_file = (find_one(os.path.join(rdir, "dada2", "ASV_tax_species.*.tsv"))
                        or find_one(os.path.join(rdir, "dada2", "ASV_tax.*.tsv")))
            n_asv = 0
            genera = set()

            if not os.path.isfile(asv_table):
                print(f"  {region}: WARNING no dada2/ASV_table.tsv -- ASV tables skipped")
                skipped.append(region)
            else:
                tax = {}
                if tax_file:
                    thdr, trows = read_tsv(tax_file)
                    tidx = {c: i for i, c in enumerate(thdr)}
                    for row in trows:
                        asv = row[tidx["ASV_ID"]]
                        tax[asv] = [
                            (row[tidx[r]] if r in tidx and tidx[r] < len(row) else "")
                            for r in TAX_RANKS
                        ]
                else:
                    print(f"  {region}: WARNING no taxonomy file -- taxonomy left blank")

                hdr, rows = read_tsv(asv_table)
                samples = hdr[1:]
                genus_acc = {}   # (sample, lineage) -> count
                n_contaminant_asv = 0
                n_contaminant_reads = 0
                for row in rows:
                    asv = row[0]
                    lineage_full = tax.get(asv, [""] * len(TAX_RANKS))
                    if is_host_contaminant(lineage_full):
                        n_contaminant_asv += 1
                        n_contaminant_reads += sum(_int(v) for v in row[1:])
                        continue
                    n_asv += 1
                    # genus-level lineage, truncated at the first unassigned rank
                    parts = []
                    for val in lineage_full[:6]:
                        if not val:
                            break
                        parts.append(val)
                    lineage = ";".join(parts) if parts else "Unassigned"
                    if len(parts) == 6:
                        genera.add(lineage)
                    for j, sample in enumerate(samples, start=1):
                        c = _int(row[j] if j < len(row) else 0)
                        if c <= 0:
                            continue
                        asv_w.writerow([region, sample, asv, c] + lineage_full)
                        key = (sample, lineage)
                        genus_acc[key] = genus_acc.get(key, 0) + c
                for (sample, lineage), c in genus_acc.items():
                    genus_rows.append((region, sample, lineage, c))
                if n_contaminant_asv:
                    print(f"  {region}: excluded {n_contaminant_asv} mitochondria/chloroplast "
                          f"ASV(s), {n_contaminant_reads} reads")

            qc_rows.append({
                "region": region,
                "n_samples": n_samples,
                "reads_into_dada2": total_in,
                "reads_nonchimeric": total_out,
                "median_retained_percent": (f"{statistics.median(retained):.2f}"
                                            if retained else ""),
                "min_retained_percent": (f"{min(retained):.2f}" if retained else ""),
                "n_ASVs": n_asv,
                "n_genera": len(genera),
            })
            print(f"  {region}: {n_samples} samples, {n_asv} ASVs, {len(genera)} genera")

    # ---------------- write read tracking ----------------
    track_path = os.path.join(args.outdir, "read_tracking_by_region.tsv")
    with open(track_path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["region"] + tracking_cols)
        for region, d in tracking_rows:
            w.writerow([region] + [d.get(c, "") for c in tracking_cols])

    # ---------------- write QC summary ----------------
    qc_path = os.path.join(args.outdir, "region_qc_summary.tsv")
    qc_cols = ["region", "n_samples", "reads_into_dada2", "reads_nonchimeric",
               "median_retained_percent", "min_retained_percent", "n_ASVs", "n_genera"]
    with open(qc_path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(qc_cols)
        for r in qc_rows:
            w.writerow([r[c] for c in qc_cols])

    # ---------------- write genus tables ----------------
    genus_path = os.path.join(args.outdir, "genus_counts_by_region.tsv")
    with open(genus_path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["region", "sampleID", "genus_lineage", "count"])
        for row in sorted(genus_rows):
            w.writerow(row)

    matrix_path = os.path.join(args.outdir, "genus_matrix_across_regions.tsv")
    region_order = [r["region"] for r in qc_rows]
    totals = {}
    for region, _sample, lineage, count in genus_rows:
        totals.setdefault(lineage, {}).update(
            {region: totals.get(lineage, {}).get(region, 0) + count})
    with open(matrix_path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["genus_lineage"] + region_order + ["n_regions_detected", "total_count"])
        for lineage in sorted(totals, key=lambda k: -sum(totals[k].values())):
            vals = [totals[lineage].get(r, 0) for r in region_order]
            w.writerow([lineage] + vals + [sum(1 for v in vals if v > 0), sum(vals)])

    # ---------------- carry the demux summary forward ----------------
    demux_out = os.path.join(args.outdir, "demux_summary.tsv")
    if os.path.isfile(args.demux_counts):
        with open(args.demux_counts) as src, open(demux_out, "w") as dst:
            dst.write(src.read())
    else:
        demux_out = None

    print()
    print("Wrote:")
    for p in [track_path, qc_path, asv_long_path, genus_path, matrix_path, demux_out]:
        if p:
            print(f"  {p}")
    if skipped:
        print(f"\nWARNING: no ASV table found for: {', '.join(skipped)} -- these "
              f"regions contribute read tracking only.")

    print()
    print("=== region QC summary ===")
    widths = [max(len(c), max((len(str(r[c])) for r in qc_rows), default=0))
              for c in qc_cols]
    print("  ".join(c.ljust(w) for c, w in zip(qc_cols, widths)))
    for r in qc_rows:
        print("  ".join(str(r[c]).ljust(w) for c, w in zip(qc_cols, widths)))


def _int(v):
    try:
        return int(float(v))
    except (TypeError, ValueError):
        return 0


def _float(v):
    try:
        return float(v)
    except (TypeError, ValueError):
        return None


if __name__ == "__main__":
    main()
