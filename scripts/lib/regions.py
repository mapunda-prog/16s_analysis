#!/usr/bin/env python3
"""
Shared helpers: region config parsing, sample-ID cleaning, FASTQ pair discovery.

Also usable as a CLI:

    # region names, one per line
    python3 scripts/lib/regions.py regions config/regions.tsv

    # write cutadapt primer FASTAs + readthrough adapter table
    python3 scripts/lib/regions.py primers config/regions.tsv --outdir config

    # tab-separated: sample_id, R1, R2
    python3 scripts/lib/regions.py pairs /path/to/fastq

    # one field of one region: V3V4 trunclenf
    python3 scripts/lib/regions.py field config/regions.tsv V3V4 trunclenf
"""

import argparse
import os
import re
import sys
from collections import OrderedDict, defaultdict

# --------------------------------------------------------------------------- #
# region config
# --------------------------------------------------------------------------- #

REQUIRED_COLUMNS = [
    "region", "fwd_name", "fwd_seq", "rev_name", "rev_seq", "amplicon_bp",
    "trunclenf", "trunclenr", "ref_taxonomy", "trim_readthrough", "extra_args",
]

IUPAC = set("ACGTRYSWKMBDHVN")


def read_regions(path):
    """Parse config/regions.tsv -> OrderedDict[region] = {column: value}."""
    rows = OrderedDict()
    header = None
    with open(path) as fh:
        for lineno, raw in enumerate(fh, 1):
            line = raw.rstrip("\n\r")
            if not line.strip() or line.lstrip().startswith("#"):
                continue
            fields = line.split("\t")
            if header is None:
                header = fields
                missing = [c for c in REQUIRED_COLUMNS if c not in header]
                if missing:
                    sys.exit(f"{path}: missing required column(s): {', '.join(missing)}")
                continue
            if len(fields) != len(header):
                sys.exit(f"{path}:{lineno}: expected {len(header)} tab-separated "
                         f"fields, found {len(fields)}")
            row = dict(zip(header, fields))
            region = row["region"].strip()
            if not re.fullmatch(r"[A-Za-z0-9_]+", region):
                sys.exit(f"{path}:{lineno}: region name '{region}' must be "
                         f"letters/digits/underscore only (it becomes a directory name)")
            if region in rows:
                sys.exit(f"{path}:{lineno}: duplicate region '{region}'")
            for key in ("fwd_seq", "rev_seq"):
                seq = row[key].strip().upper()
                if not seq:
                    sys.exit(f"{path}:{lineno}: {key} is empty for {region}")
                bad = set(seq) - IUPAC
                if bad:
                    sys.exit(f"{path}:{lineno}: {key} for {region} has non-IUPAC "
                             f"character(s): {''.join(sorted(bad))}")
                row[key] = seq
            rows[region] = row
    if header is None:
        sys.exit(f"{path}: no header row found")
    if not rows:
        sys.exit(f"{path}: no region rows found")
    return rows


REVCOMP = str.maketrans("ACGTRYSWKMBDHVNacgtryswkmbdhvn",
                        "TGCAYRSWMKVHDBNtgcayrswmkvhdbn")


def revcomp(seq):
    return seq.translate(REVCOMP)[::-1]


def write_primer_files(regions, outdir):
    """Write cutadapt inputs. Order is identical in both FASTAs, which is what
    cutadapt --pair-adapters requires: i-th forward pairs with i-th reverse."""
    os.makedirs(outdir, exist_ok=True)
    fwd_path = os.path.join(outdir, "primers_fwd.fasta")
    rev_path = os.path.join(outdir, "primers_rev.fasta")
    rt_path = os.path.join(outdir, "readthrough.tsv")

    with open(fwd_path, "w") as fw, open(rev_path, "w") as rv:
        for region, row in regions.items():
            fw.write(f">{region}\n{row['fwd_seq']}\n")
            rv.write(f">{region}\n{row['rev_seq']}\n")

    with open(rt_path, "w") as fh:
        fh.write("region\ttrim_readthrough\tadapter_r1\tadapter_r2\n")
        for region, row in regions.items():
            fh.write("\t".join([
                region,
                row["trim_readthrough"].strip().lower(),
                revcomp(row["rev_seq"]),   # -a : 3' of R1 is revcomp of reverse primer
                revcomp(row["fwd_seq"]),   # -A : 3' of R2 is revcomp of forward primer
            ]) + "\n")

    shortest = min(len(r["fwd_seq"]) for r in regions.values())
    shortest = min(shortest, min(len(r["rev_seq"]) for r in regions.values()))
    return fwd_path, rev_path, rt_path, shortest


# --------------------------------------------------------------------------- #
# sample IDs and FASTQ pairing
# --------------------------------------------------------------------------- #

# Illumina: <sample>_S<num>_L<lane>_R<1|2>_001.fastq.gz
ILLUMINA_RE = re.compile(
    r"^(?P<sample>.+)_S(?P<snum>\d+)_L(?P<lane>\d+)_R(?P<read>[12])_001\.f(?:ast)?q\.gz$"
)
# Already-clean: <sample>_R<1|2>.fastq.gz  (what step 01 emits)
SIMPLE_RE = re.compile(
    r"^(?P<sample>.+)_R(?P<read>[12])\.f(?:ast)?q\.gz$"
)


def clean_id(raw):
    """MLTP-AIC-126-SW -> AIC126SW.

    ampliseq requires IDs that are unique, start with a letter, and contain only
    letters/digits/underscores. Kept identical to make_ampliseq_samplesheet.py so
    IDs stay comparable with the earlier pooled run.
    """
    s = raw
    if s.startswith("MLTP-"):
        s = s[len("MLTP-"):]
    s = s.replace("-", "")
    s = re.sub(r"[^A-Za-z0-9_]", "_", s)
    if not s or not s[0].isalpha():
        s = "S_" + s
    return s


def find_pairs(indir, clean=True):
    """Discover paired FASTQs in indir.

    Returns [(sample_id, r1_path, r2_path)] sorted by sample_id. Lanes belonging
    to the same (sample, S-index) are an error rather than a silent pick, because
    picking one would quietly drop data.
    """
    indir = os.path.abspath(indir)
    if not os.path.isdir(indir):
        sys.exit(f"Not a directory: {indir}")

    groups = defaultdict(dict)      # (sample, snum) -> {"1": path, "2": path}
    lanes = defaultdict(set)
    unparsed = []

    for fname in sorted(os.listdir(indir)):
        if not (fname.endswith(".fastq.gz") or fname.endswith(".fq.gz")):
            continue
        m = ILLUMINA_RE.match(fname)
        if m:
            key = (m["sample"], m["snum"])
            lanes[key].add(m["lane"])
        else:
            m = SIMPLE_RE.match(fname)
            if not m:
                unparsed.append(fname)
                continue
            key = (m["sample"], "")
        groups[key][m["read"]] = os.path.join(indir, fname)

    if unparsed:
        sys.stderr.write("WARNING: skipped files not matching a known FASTQ "
                         "naming pattern:\n")
        for f in unparsed:
            sys.stderr.write(f"    {f}\n")

    multilane = {k: v for k, v in lanes.items() if len(v) > 1}
    if multilane:
        msg = ", ".join(f"{k[0]}_S{k[1]} (lanes {sorted(v)})"
                        for k, v in sorted(multilane.items()))
        sys.exit("ERROR: these samples are split across multiple lanes: " + msg +
                 "\nConcatenate the lanes per read before demultiplexing, e.g.\n"
                 "  cat X_S1_L00*_R1_001.fastq.gz > merged/X_S1_L001_R1_001.fastq.gz")

    base_counts = defaultdict(int)
    for sample, _snum in groups:
        base_counts[clean_id(sample) if clean else sample] += 1

    out = []
    for (sample, snum), reads in sorted(groups.items()):
        base = clean_id(sample) if clean else sample
        # disambiguate same-name samples by Illumina S-index, as before
        sid = base if base_counts[base] == 1 or not snum else f"{base}_S{snum}"
        r1, r2 = reads.get("1"), reads.get("2")
        if not r1 or not r2:
            sys.stderr.write(f"WARNING: {sample} (S{snum or '-'}) is not paired "
                             f"(R1={bool(r1)} R2={bool(r2)}) -- skipping. "
                             f"This pipeline requires paired-end reads.\n")
            continue
        out.append((sid, r1, r2))

    seen = set()
    for sid, _, _ in out:
        if sid in seen:
            sys.exit(f"ERROR: duplicate sample ID after cleaning: {sid}")
        seen.add(sid)
    return out


# --------------------------------------------------------------------------- #
# CLI
# --------------------------------------------------------------------------- #

def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)

    p = sub.add_parser("regions", help="print region names")
    p.add_argument("config")

    p = sub.add_parser("primers", help="write cutadapt primer FASTAs")
    p.add_argument("config")
    p.add_argument("--outdir", required=True)

    p = sub.add_parser("pairs", help="print sample_id<TAB>R1<TAB>R2")
    p.add_argument("indir")
    p.add_argument("--no-clean", action="store_true",
                   help="keep original sample names instead of cleaning them")

    p = sub.add_parser("field", help="print one config field")
    p.add_argument("config")
    p.add_argument("region")
    p.add_argument("column")

    args = ap.parse_args()

    if args.cmd == "regions":
        for r in read_regions(args.config):
            print(r)

    elif args.cmd == "primers":
        regions = read_regions(args.config)
        fwd, rev, rt, shortest = write_primer_files(regions, args.outdir)
        print(f"regions:          {len(regions)}", file=sys.stderr)
        print(f"forward primers:  {fwd}", file=sys.stderr)
        print(f"reverse primers:  {rev}", file=sys.stderr)
        print(f"readthrough:      {rt}", file=sys.stderr)
        print(shortest)   # stdout: shortest primer length, used as cutadapt -O

    elif args.cmd == "pairs":
        for sid, r1, r2 in find_pairs(args.indir, clean=not args.no_clean):
            print(f"{sid}\t{r1}\t{r2}")

    elif args.cmd == "field":
        regions = read_regions(args.config)
        if args.region not in regions:
            sys.exit(f"unknown region '{args.region}'; have: "
                     f"{', '.join(regions)}")
        row = regions[args.region]
        if args.column not in row:
            sys.exit(f"unknown column '{args.column}'")
        print(row[args.column])


if __name__ == "__main__":
    main()
