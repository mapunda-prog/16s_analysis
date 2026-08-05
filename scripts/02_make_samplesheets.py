#!/usr/bin/env python3
"""
Step 02 -- build one nf-core/ampliseq samplesheet per primer region.

Reads the per-region FASTQs written by step 01 and emits
samplesheets/<REGION>.tsv with the columns ampliseq expects
(sampleID / forwardReads / reverseReads, absolute paths).

Samples with very few read pairs are excluded per region rather than left in:
DADA2 either errors or produces meaningless error models on near-empty samples,
and a single bad sample fails the whole region's run. Everything excluded is
listed in samplesheets/excluded_samples.tsv, so nothing disappears silently.

Usage:
    python3 scripts/02_make_samplesheets.py --demux-dir demux
    python3 scripts/02_make_samplesheets.py --demux-dir demux --min-reads 500
"""

import argparse
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "lib"))
import regions as R  # noqa: E402


def main():
    here = os.path.dirname(os.path.abspath(__file__))
    project = os.path.dirname(here)

    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--demux-dir", default=os.path.join(project, "demux"),
                    help="output root from step 01 [demux]")
    ap.add_argument("--config", default=os.path.join(project, "config", "regions.tsv"),
                    help="region table [config/regions.tsv]")
    ap.add_argument("--outdir", default=os.path.join(project, "samplesheets"),
                    help="where to write samplesheets [samplesheets]")
    ap.add_argument("--min-reads", type=int, default=1000,
                    help="exclude a sample from a region below this many read "
                         "pairs [1000]; use 0 to keep everything")
    ap.add_argument("--layout", choices=["legacy", "standardized"], default="legacy",
                    help="legacy = sampleID/forwardReads/reverseReads (ampliseq "
                         "native), standardized = sample/fastq_1/fastq_2")
    args = ap.parse_args()

    cfg = R.read_regions(args.config)
    demux = os.path.abspath(args.demux_dir)
    if not os.path.isdir(demux):
        sys.exit(f"demux dir not found: {demux}\nRun scripts/01_demux_regions.sh first.")
    os.makedirs(args.outdir, exist_ok=True)

    header = (["sampleID", "forwardReads", "reverseReads"]
              if args.layout == "legacy" else ["sample", "fastq_1", "fastq_2"])

    excluded = []
    written = []

    print(f"{'region':10} {'samples':>8} {'excluded':>9}   samplesheet")
    for region in cfg:
        rdir = os.path.join(demux, region)
        if not os.path.isdir(rdir):
            print(f"{region:10} {'-':>8} {'-':>9}   WARNING: {rdir} missing, skipped")
            continue

        pairs = R.find_pairs(rdir, clean=False)
        keep, drop = [], []
        for sid, r1, r2 in pairs:
            n = count_read_pairs(r1)
            if n < args.min_reads:
                drop.append((sid, n))
            else:
                keep.append((sid, r1, r2, n))

        sheet = os.path.join(args.outdir, f"{region}.tsv")
        with open(sheet, "w") as fh:
            fh.write("\t".join(header) + "\n")
            for sid, r1, r2, _n in keep:
                fh.write("\t".join([sid, r1, r2]) + "\n")

        for sid, n in drop:
            excluded.append((region, sid, n))
        print(f"{region:10} {len(keep):8} {len(drop):9}   {sheet}")
        if keep:
            written.append(region)
        else:
            print(f"{'':10} {'':8} {'':9}   WARNING: no samples passed "
                  f"--min-reads {args.min_reads}; this region cannot be run")

    exc_path = os.path.join(args.outdir, "excluded_samples.tsv")
    with open(exc_path, "w") as fh:
        fh.write("region\tsampleID\tread_pairs\tmin_reads_threshold\n")
        for region, sid, n in excluded:
            fh.write(f"{region}\t{sid}\t{n}\t{args.min_reads}\n")

    print()
    print(f"Excluded {len(excluded)} sample x region combination(s) -> {exc_path}")
    print(f"Runnable regions: {', '.join(written) if written else '(none)'}")
    print()
    print("Next: scripts/03_run_ampliseq.sh --all")


def count_read_pairs(fastq_gz):
    """Line count / 4. Uses the counts file from step 01 when available, since
    re-decompressing every FASTQ is the slow part of this script."""
    cached = _counts_cache(fastq_gz)
    if cached is not None:
        return cached
    import gzip
    n = 0
    with gzip.open(fastq_gz, "rb") as fh:
        for _ in fh:
            n += 1
    return n // 4


_CACHE = None


def _counts_cache(fastq_gz):
    """Look up demux/demux_counts.tsv written by step 01."""
    global _CACHE
    if _CACHE is None:
        _CACHE = {}
        # demux/<REGION>/<SID>_R1.fastq.gz -> demux/demux_counts.tsv
        region_dir = os.path.dirname(os.path.abspath(fastq_gz))
        counts = os.path.join(os.path.dirname(region_dir), "demux_counts.tsv")
        if os.path.isfile(counts):
            with open(counts) as fh:
                next(fh, None)
                for line in fh:
                    parts = line.rstrip("\n").split("\t")
                    if len(parts) == 3:
                        _CACHE[(parts[1], parts[0])] = int(parts[2])
    if not _CACHE:
        return None
    base = os.path.basename(fastq_gz)
    sid = base[:-len("_R1.fastq.gz")] if base.endswith("_R1.fastq.gz") else base
    region = os.path.basename(os.path.dirname(os.path.abspath(fastq_gz)))
    return _CACHE.get((region, sid))


if __name__ == "__main__":
    main()
