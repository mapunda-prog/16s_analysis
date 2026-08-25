#!/usr/bin/env python3
"""
Step 05 -- alpha and beta diversity per region, grouped by specimen type
(type_ofsample: Blood/Swab) and collection site (hf_name).

Reads compiled/asv_long.tsv.gz (written by 04_compile_results.py -- already
free of the mitochondria/chloroplast ASVs that dominate some blood specimens,
see that script's docstring) and config/metadata.tsv, then per region:

Deliberately NOT sourced from each region's qiime2/abundance_tables/
feature-table.tsv: that file only exists if QIIME2_BARPLOT succeeded for that
region, which it did not for every region in at least one real run (ITS and
V5V7 both had QIIME2_BARPLOT fail while DADA2 itself completed fine) --
asv_long.tsv.gz has no such dependency.

  1. Drops samples below --min-depth reads (default 1000, matching the
     --min-reads convention already used in step 02) -- comparing alpha
     diversity across wildly different depths is not meaningful. Every
     dropped sample is recorded, not silently discarded.
  2. Rarefies the remaining samples to the lowest surviving depth in that
     region (single draw, fixed seed -- reproducible, not "the" answer, since
     rarefaction has draw-to-draw noise, but noise a fair reference point).
  3. Computes per-sample alpha diversity: observed richness, Shannon (natural
     log), Simpson, Pielou's evenness.
  4. Computes a Bray-Curtis distance matrix on the rarefied table, and a
     PERMANOVA-style permutation test (Anderson 2001 pseudo-F, 999
     permutations) for each grouping variable that has >=2 groups of >=2
     samples.

Requires numpy (not part of the server pipeline's dependencies) -- run this
locally against a downloaded compiled/, same as build_metadata.py.

Usage:
    python3 scripts/05_diversity_analysis.py \\
        --asv-long /path/to/compiled/asv_long.tsv.gz \\
        --metadata config/metadata.tsv \\
        --outdir /path/to/alpha_beta_diversity
"""

import argparse
import csv
import gzip
import os
import statistics
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
PROJECT = os.path.dirname(HERE)

GROUP_COLUMNS = ["type_ofsample", "hf_name"]
SEED = 100
N_PERM = 999


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--asv-long", default=os.path.join(PROJECT, "compiled", "asv_long.tsv.gz"))
    ap.add_argument("--metadata", default=os.path.join(PROJECT, "config", "metadata.tsv"))
    ap.add_argument("--outdir", default=os.path.join(PROJECT, "alpha_beta_diversity"))
    ap.add_argument("--min-depth", type=int, default=1000,
                    help="samples below this many reads are excluded [1000]")
    args = ap.parse_args()

    if not os.path.isfile(args.asv_long):
        sys.exit(f"not found: {args.asv_long}\nRun scripts/04_compile_results.py first.")
    metadata = read_metadata(args.metadata)
    by_region = read_asv_long(args.asv_long)
    if not by_region:
        sys.exit(f"no rows in {args.asv_long}")
    regions = sorted(by_region)

    os.makedirs(args.outdir, exist_ok=True)
    os.makedirs(os.path.join(args.outdir, "beta_distance_matrices"), exist_ok=True)

    alpha_rows = []
    excluded_rows = []
    group_test_rows = []
    rng = np.random.default_rng(SEED)

    print(f"{'region':8} {'samples':>8} {'excluded':>9} {'rarefy_depth':>13}")
    for region in regions:
        samples, asv_ids, counts = pivot_region(by_region[region])  # counts: samples x ASVs

        depths = counts.sum(axis=1)
        keep_mask = depths >= args.min_depth
        for sid, d in zip(samples, depths):
            if d < args.min_depth:
                excluded_rows.append({"region": region, "sample-id": sid,
                                      "reads": int(d), "min_depth": args.min_depth})

        kept_samples = [s for s, k in zip(samples, keep_mask) if k]
        kept_counts = counts[keep_mask]
        n_excluded = len(samples) - len(kept_samples)

        if len(kept_samples) < 2:
            print(f"{region:8} {len(samples):8} {n_excluded:9} {'--':>13}   "
                  f"fewer than 2 samples above --min-depth, skipped")
            continue

        rarefy_depth = int(kept_counts.sum(axis=1).min())
        rarefied = np.array([rarefy(row, rarefy_depth, rng) for row in kept_counts])

        print(f"{region:8} {len(samples):8} {n_excluded:9} {rarefy_depth:13}")

        richness, shannon, simpson, pielou = alpha_metrics(rarefied)
        for i, sid in enumerate(kept_samples):
            meta = metadata.get(sid, {})
            alpha_rows.append({
                "region": region, "sample-id": sid,
                "type_ofsample": meta.get("type_ofsample", "") or "Unknown",
                "hf_name": meta.get("hf_name", "") or "Unknown",
                "reads_used_raw": int(depths[samples.index(sid)]),
                "rarefied_depth": rarefy_depth,
                "observed_richness": richness[i],
                "shannon": round(shannon[i], 4),
                "simpson": round(simpson[i], 4),
                "pielou_evenness": "" if pielou[i] is None else round(pielou[i], 4),
            })

        dist = bray_curtis_matrix(rarefied)
        write_distance_matrix(dist, kept_samples,
                              os.path.join(args.outdir, "beta_distance_matrices",
                                           f"{region}_bray_curtis.tsv"))

        for group_col in GROUP_COLUMNS:
            labels = [metadata.get(s, {}).get(group_col, "") for s in kept_samples]
            group_test_rows.append(permanova_row(region, group_col, labels, dist))

    write_alpha_tables(alpha_rows, args.outdir)
    write_excluded(excluded_rows, args.outdir)
    write_group_tests(group_test_rows, args.outdir)

    print()
    print("Wrote:")
    for name in ("alpha_diversity_per_sample.tsv", "alpha_diversity_summary_by_type.tsv",
                 "alpha_diversity_summary_by_site.tsv", "beta_diversity_group_tests.tsv",
                 "excluded_samples.tsv"):
        print(f"  {os.path.join(args.outdir, name)}")
    print(f"  {os.path.join(args.outdir, 'beta_distance_matrices', '<REGION>_bray_curtis.tsv')}")


# --------------------------------------------------------------------------- #
# I/O
# --------------------------------------------------------------------------- #

def read_asv_long(path):
    """compiled/asv_long.tsv.gz -> dict[region][(sample, ASV_ID)] = count.
    Only nonzero cells are stored in the file; missing = 0."""
    by_region = {}
    with gzip.open(path, "rt", newline="") as fh:
        reader = csv.reader(fh, delimiter="\t")
        header = next(reader)
        idx = {c: i for i, c in enumerate(header)}
        for row in reader:
            region = row[idx["region"]]
            sample = row[idx["sampleID"]]
            asv = row[idx["ASV_ID"]]
            count = float(row[idx["count"]])
            by_region.setdefault(region, {})[(sample, asv)] = count
    return by_region


def pivot_region(cell_dict):
    """{(sample, asv): count} -> (samples, asv_ids, samples x ASVs dense matrix)."""
    samples = sorted({s for s, _ in cell_dict})
    asv_ids = sorted({a for _, a in cell_dict})
    s_idx = {s: i for i, s in enumerate(samples)}
    a_idx = {a: i for i, a in enumerate(asv_ids)}
    matrix = np.zeros((len(samples), len(asv_ids)))
    for (s, a), c in cell_dict.items():
        matrix[s_idx[s], a_idx[a]] = c
    return samples, asv_ids, matrix


def read_metadata(path):
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        return {row["sample-id"]: row for row in reader}


def write_distance_matrix(dist, samples, path):
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow([""] + samples)
        for sid, row in zip(samples, dist):
            w.writerow([sid] + [f"{v:.6f}" for v in row])


def write_alpha_tables(rows, outdir):
    cols = ["region", "sample-id", "type_ofsample", "hf_name", "reads_used_raw",
            "rarefied_depth", "observed_richness", "shannon", "simpson", "pielou_evenness"]
    with open(os.path.join(outdir, "alpha_diversity_per_sample.tsv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t")
        w.writeheader()
        for row in rows:
            w.writerow(row)

    for group_col, out_name in (("type_ofsample", "alpha_diversity_summary_by_type.tsv"),
                                 ("hf_name", "alpha_diversity_summary_by_site.tsv")):
        summary = summarize_alpha(rows, group_col)
        with open(os.path.join(outdir, out_name), "w", newline="") as fh:
            w = csv.writer(fh, delimiter="\t", lineterminator="\n")
            w.writerow(["region", group_col, "n_samples",
                        "richness_mean", "richness_sd",
                        "shannon_mean", "shannon_sd",
                        "simpson_mean", "simpson_sd",
                        "pielou_mean", "pielou_sd"])
            for row in summary:
                w.writerow(row)


def summarize_alpha(rows, group_col):
    groups = {}
    for r in rows:
        key = (r["region"], r[group_col])
        groups.setdefault(key, []).append(r)
    out = []
    for (region, group_val), items in sorted(groups.items()):
        def stats(field):
            vals = [it[field] for it in items if it[field] != ""]
            if not vals:
                return "", ""
            if len(vals) < 2:
                return round(statistics.mean(vals), 4), ""
            return round(statistics.mean(vals), 4), round(statistics.stdev(vals), 4)
        r_mean, r_sd = stats("observed_richness")
        sh_mean, sh_sd = stats("shannon")
        si_mean, si_sd = stats("simpson")
        pi_mean, pi_sd = stats("pielou_evenness")
        out.append([region, group_val, len(items), r_mean, r_sd, sh_mean, sh_sd,
                    si_mean, si_sd, pi_mean, pi_sd])
    return out


def write_excluded(rows, outdir):
    with open(os.path.join(outdir, "excluded_samples.tsv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["region", "sample-id", "reads", "min_depth"],
                           delimiter="\t")
        w.writeheader()
        for row in rows:
            w.writerow(row)


def write_group_tests(rows, outdir):
    cols = ["region", "grouping_variable", "n_groups", "n_samples",
            "pseudo_F", "p_value", "n_permutations", "note"]
    with open(os.path.join(outdir, "beta_diversity_group_tests.tsv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t")
        w.writeheader()
        for row in rows:
            w.writerow(row)


# --------------------------------------------------------------------------- #
# math
# --------------------------------------------------------------------------- #

def rarefy(counts, depth, rng):
    """Single-draw subsample of an integer count vector down to `depth` reads."""
    total = int(counts.sum())
    if total == depth:
        return counts.copy()
    pool = np.repeat(np.arange(len(counts)), counts.astype(int))
    draw = rng.choice(pool, size=depth, replace=False)
    out = np.zeros(len(counts))
    idx, freq = np.unique(draw, return_counts=True)
    out[idx] = freq
    return out


def alpha_metrics(rarefied):
    richness = (rarefied > 0).sum(axis=1)
    totals = rarefied.sum(axis=1)
    p = np.divide(rarefied, totals[:, None], out=np.zeros_like(rarefied), where=totals[:, None] > 0)
    with np.errstate(divide="ignore", invalid="ignore"):
        shannon = -np.nansum(np.where(p > 0, p * np.log(p), 0), axis=1)
    simpson = 1 - np.nansum(p ** 2, axis=1)
    pielou = [
        (shannon[i] / np.log(richness[i])) if richness[i] > 1 else None
        for i in range(len(richness))
    ]
    return richness.tolist(), shannon.tolist(), simpson.tolist(), pielou


def bray_curtis_matrix(counts):
    n = counts.shape[0]
    dist = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            num = np.abs(counts[i] - counts[j]).sum()
            den = (counts[i] + counts[j]).sum()
            d = num / den if den > 0 else 0.0
            dist[i, j] = dist[j, i] = d
    return dist


def permanova_row(region, group_col, labels, dist):
    groups = {}
    for i, lab in enumerate(labels):
        if lab == "":
            continue
        groups.setdefault(lab, []).append(i)
    valid_groups = {g: idx for g, idx in groups.items() if len(idx) >= 2}
    n_used = sum(len(v) for v in valid_groups.values())

    if len(valid_groups) < 2:
        return {"region": region, "grouping_variable": group_col,
                "n_groups": len(valid_groups), "n_samples": n_used,
                "pseudo_F": "", "p_value": "",
                "n_permutations": 0, "note": "fewer than 2 groups with >=2 samples"}

    idx_all = [i for idxs in valid_groups.values() for i in idxs]
    labels_all = np.array([g for g, idxs in valid_groups.items() for _ in idxs])
    sub = dist[np.ix_(idx_all, idx_all)]

    f_obs = pseudo_f(sub, labels_all)
    rng = np.random.default_rng(SEED)
    count_ge = 0
    for _ in range(N_PERM):
        perm = rng.permutation(labels_all)
        if pseudo_f(sub, perm) >= f_obs:
            count_ge += 1
    p_value = (count_ge + 1) / (N_PERM + 1)

    return {"region": region, "grouping_variable": group_col,
            "n_groups": len(valid_groups), "n_samples": n_used,
            "pseudo_F": round(f_obs, 4), "p_value": round(p_value, 4),
            "n_permutations": N_PERM, "note": ""}


def pseudo_f(dist_sq_source, labels):
    """Anderson (2001) PERMANOVA pseudo-F from a Bray-Curtis distance matrix."""
    d2 = dist_sq_source ** 2
    n = len(labels)
    ss_total = d2.sum() / (2 * n)
    ss_within = 0.0
    for g in np.unique(labels):
        idx = np.where(labels == g)[0]
        if len(idx) < 2:
            continue
        sub = d2[np.ix_(idx, idx)]
        ss_within += sub.sum() / (2 * len(idx))
    ss_between = ss_total - ss_within
    a = len(np.unique(labels))
    df_between = a - 1
    df_within = n - a
    if df_within <= 0 or ss_within <= 0:
        return 0.0
    return (ss_between / df_between) / (ss_within / df_within)


if __name__ == "__main__":
    main()
