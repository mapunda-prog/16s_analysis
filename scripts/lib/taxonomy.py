#!/usr/bin/env python3
"""Shared helpers: taxonomy ranks, host-contaminant filtering, and reading a
region's per-ASV table (counts + taxonomy + classifier confidence) in a form
that works the same way whether that region used SILVA (bacteria, e.g.
V3V4/V4V5) or UNITE (fungi, ITS) as its reference -- both write the same
dada2/ASV_table.tsv and dada2/ASV_tax_species.<ref>.tsv column layout."""

import csv
import glob
import os

TAX_RANKS = ["Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Species_exact"]

# ampliseq's own default (--exclude_taxa mitochondria,chloroplast) applies only
# to the QIIME2 abundance tables; ASV_table.tsv is upstream of that filter, so
# callers working from ASV_table.tsv directly must re-apply the same rule to
# stay consistent with the pipeline's own filtered results.
HOST_CONTAMINANT_TERMS = ("mitochondria", "chloroplast")


def is_host_contaminant(lineage_full):
    return any(term in val.lower() for val in lineage_full for term in HOST_CONTAMINANT_TERMS)


def read_tsv(path):
    """-> (header list, list of row lists). Skips biom-style '# Constructed' lines."""
    with open(path, newline="") as fh:
        rows = [r for r in csv.reader(fh, delimiter="\t")
                if r and not (len(r) == 1 and r[0].startswith("# "))]
    if not rows:
        return [], []
    return rows[0], rows[1:]


def load_region_asvs(region_dir):
    """Read dada2/ASV_table.tsv + dada2/ASV_tax_species.<ref>.tsv for one
    region. Returns (samples, asvs) where asvs is a list of dicts:
        {asv_id, tax: {rank: value}, confidence: str, counts: {sample: int}}
    Host-contaminant ASVs (Mitochondria/Chloroplast, any rank) are excluded.
    Raises FileNotFoundError with a clear message if either input is missing.
    """
    table_path = os.path.join(region_dir, "dada2", "ASV_table.tsv")
    tax_path = _find_one(os.path.join(region_dir, "dada2", "ASV_tax_species.*.tsv"))
    if not os.path.isfile(table_path):
        raise FileNotFoundError(f"no dada2/ASV_table.tsv under {region_dir}")
    if not tax_path:
        raise FileNotFoundError(f"no dada2/ASV_tax_species.*.tsv under {region_dir}")

    thdr, trows = read_tsv(tax_path)
    tidx = {c: i for i, c in enumerate(thdr)}
    tax_by_asv = {}
    conf_by_asv = {}
    for row in trows:
        asv_id = row[tidx["ASV_ID"]]
        tax_by_asv[asv_id] = {r: (row[tidx[r]] if r in tidx and tidx[r] < len(row) else "")
                              for r in TAX_RANKS}
        conf_by_asv[asv_id] = row[tidx["confidence"]] if "confidence" in tidx else ""

    hdr, rows = read_tsv(table_path)
    samples = hdr[1:]
    asvs = []
    for row in rows:
        asv_id = row[0]
        tax = tax_by_asv.get(asv_id, {r: "" for r in TAX_RANKS})
        if is_host_contaminant(list(tax.values())):
            continue
        counts = {s: _int(row[j]) for j, s in enumerate(samples, start=1) if j < len(row)}
        asvs.append({"asv_id": asv_id, "tax": tax, "confidence": conf_by_asv.get(asv_id, ""),
                    "counts": counts})
    return samples, asvs


def _find_one(pattern):
    hits = sorted(glob.glob(pattern))
    return hits[0] if hits else None


def _int(v):
    try:
        return int(float(v))
    except (TypeError, ValueError):
        return 0
