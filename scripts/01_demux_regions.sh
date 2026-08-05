#!/usr/bin/env bash
#
# Step 01 -- split each sample's reads into the 7 QIAseq primer regions.
#
# One cutadapt pass per sample assigns each read PAIR to a region by matching the
# forward primer at the 5' end of R1 and the matching reverse primer at the 5' end
# of R2 (cutadapt --pair-adapters), then trims both. The 5' match is UNANCHORED,
# which is what makes this work with QIAseq phased primers: the 0-11 phasing bases
# sit upstream of the primer and a 5' trim removes everything up to and including
# the primer match. Pairs where both primers do not agree on a region land in
# demux/unknown/ for QC rather than being thrown away.
#
# Usage:
#   scripts/01_demux_regions.sh --indir <fastq_dir> [options]
#
#   --indir DIR        directory of paired FASTQs (required)
#   --outdir DIR       output root                      [demux]
#   --config FILE      region table                     [config/regions.tsv]
#   --cores N          cutadapt cores                   [8]
#   --error-rate F     cutadapt -e                      [0.15]
#   --min-length N     drop reads shorter than this after trimming [100]
#   --min-overlap N    cutadapt -O; default = shortest primer length
#   --swap-reads       treat the R2 file as R1 (see --probe)
#   --probe [N]        dry test on N read pairs (default 100000) of one sample,
#                      in both orientations, then exit
#   --force            overwrite existing per-sample output
#
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"
LIB="$SCRIPT_DIR/lib/regions.py"

INDIR=""
OUTDIR="$PROJECT_DIR/demux"
CONFIG="$PROJECT_DIR/config/regions.tsv"
CORES=8
ERROR_RATE=0.15
MIN_LENGTH=100
MIN_OVERLAP=""
SWAP_READS=0
PROBE=0
PROBE_N=100000
FORCE=0

while [[ $# -gt 0 ]]; do
    case "$1" in
        --indir)       INDIR="$2"; shift 2 ;;
        --outdir)      OUTDIR="$2"; shift 2 ;;
        --config)      CONFIG="$2"; shift 2 ;;
        --cores)       CORES="$2"; shift 2 ;;
        --error-rate)  ERROR_RATE="$2"; shift 2 ;;
        --min-length)  MIN_LENGTH="$2"; shift 2 ;;
        --min-overlap) MIN_OVERLAP="$2"; shift 2 ;;
        --swap-reads)  SWAP_READS=1; shift ;;
        --force)       FORCE=1; shift ;;
        --probe)
            PROBE=1; shift
            if [[ $# -gt 0 && "$1" =~ ^[0-9]+$ ]]; then PROBE_N="$1"; shift; fi
            ;;
        -h|--help)     sed -n '2,30p' "${BASH_SOURCE[0]}"; exit 0 ;;
        *) echo "Unknown option: $1" >&2; exit 2 ;;
    esac
done

[[ -n "$INDIR" ]] || { echo "ERROR: --indir is required (see --help)" >&2; exit 2; }
[[ -d "$INDIR" ]] || { echo "ERROR: not a directory: $INDIR" >&2; exit 2; }
[[ -f "$CONFIG" ]] || { echo "ERROR: config not found: $CONFIG" >&2; exit 2; }

for tool in cutadapt python3; do
    command -v "$tool" >/dev/null || { echo "ERROR: '$tool' not on PATH" >&2; exit 127; }
done
echo "cutadapt: $(cutadapt --version)"

# ---------------------------------------------------------------------------- #
# primer FASTAs (regenerated every run so they can never drift from the config)
# ---------------------------------------------------------------------------- #
mkdir -p "$OUTDIR"
SHORTEST_PRIMER="$(python3 "$LIB" primers "$CONFIG" --outdir "$PROJECT_DIR/config")"
FWD_FASTA="$PROJECT_DIR/config/primers_fwd.fasta"
REV_FASTA="$PROJECT_DIR/config/primers_rev.fasta"
: "${MIN_OVERLAP:=$SHORTEST_PRIMER}"

# read_lines_into VAR < input  -- portable stand-in for bash 4's mapfile, so these
# scripts also run under the bash 3.2 that ships with macOS
REGIONS=()
while IFS= read -r line; do [[ -n "$line" ]] && REGIONS+=("$line"); done \
    < <(python3 "$LIB" regions "$CONFIG")
echo "regions:  ${REGIONS[*]}"
echo "-O ${MIN_OVERLAP}  -e ${ERROR_RATE}  --minimum-length ${MIN_LENGTH}"

# cutadapt substitutes {name} but will not create directories, so pre-create them
for r in "${REGIONS[@]}" unknown; do mkdir -p "$OUTDIR/$r"; done

PAIRS=()
while IFS= read -r line; do [[ -n "$line" ]] && PAIRS+=("$line"); done \
    < <(python3 "$LIB" pairs "$INDIR")
[[ ${#PAIRS[@]} -gt 0 ]] || { echo "ERROR: no paired FASTQs found in $INDIR" >&2; exit 1; }
echo "samples:  ${#PAIRS[@]}"
echo

# ---------------------------------------------------------------------------- #
# --probe: how many pairs get assigned, in each orientation?
# ---------------------------------------------------------------------------- #
count_pairs() {
    # read pairs in a gzipped FASTQ = lines / 4
    local f="$1"
    if [[ -s "$f" ]]; then
        echo $(( $(gzip -cd "$f" | wc -l) / 4 ))
    else
        echo 0
    fi
}

if [[ $PROBE -eq 1 ]]; then
    IFS=$'\t' read -r SID R1 R2 <<< "${PAIRS[0]}"
    echo "PROBE: $SID, first $PROBE_N read pairs, both orientations"
    TMP="$(mktemp -d)"
    trap 'rm -rf "$TMP"' EXIT
    HEAD_LINES=$((PROBE_N * 4))
    for orient in as_given swapped; do
        if [[ "$orient" == "as_given" ]]; then A="$R1"; B="$R2"; else A="$R2"; B="$R1"; fi
        # `head` closing the pipe early makes gzip exit non-zero; that is expected
        gzip -cd "$A" 2>/dev/null | head -n "$HEAD_LINES" | gzip > "$TMP/a.fastq.gz" || true
        gzip -cd "$B" 2>/dev/null | head -n "$HEAD_LINES" | gzip > "$TMP/b.fastq.gz" || true

        # Demultiplex exactly as the real run does, into a throwaway tree, then
        # count per region. Counting the output is the only honest measure of the
        # assignment rate: cutadapt's "Pairs written" includes unassigned pairs.
        rm -rf "$TMP/out"
        for r in "${REGIONS[@]}" unknown; do mkdir -p "$TMP/out/$r"; done
        cutadapt \
            -g "file:$FWD_FASTA" -G "file:$REV_FASTA" \
            --pair-adapters -e "$ERROR_RATE" -O "$MIN_OVERLAP" \
            --minimum-length "$MIN_LENGTH" --cores "$CORES" \
            -o "$TMP/out/{name}/p_R1.fastq.gz" -p "$TMP/out/{name}/p_R2.fastq.gz" \
            "$TMP/a.fastq.gz" "$TMP/b.fastq.gz" > "$TMP/probe.log" 2>&1 || {
                echo "  cutadapt failed:"; sed 's/^/    /' "$TMP/probe.log"; exit 1; }

        echo
        echo "--- orientation: $orient"
        total=0
        for r in "${REGIONS[@]}" unknown; do
            n="$(count_pairs "$TMP/out/$r/p_R1.fastq.gz")"
            total=$((total + n))
            printf '%s\t%s\n' "$r" "$n"
        done > "$TMP/probe_counts.tsv"
        awk -F'\t' -v tot="$total" '
            { printf "    %-10s %10d  %6.2f%%\n", $1, $2, (tot ? 100*$2/tot : 0)
              if ($1 != "unknown") assigned += $2 }
            END { printf "    %-10s %10d  %6.2f%%\n", "ASSIGNED", assigned,
                         (tot ? 100*assigned/tot : 0) }' "$TMP/probe_counts.tsv"
    done
    cat <<'EOF'

Interpretation:
  * "as_given" should show a high ASSIGNED rate. If "swapped" is much higher,
    rerun step 01 with --swap-reads.
  * If BOTH orientations assign very little, the primers in config/regions.tsv do
    not match this library -- check your kit/panel before going further.
  * A region sitting at 0% in an otherwise healthy probe was probably not part of
    your panel; remove its row from config/regions.tsv.
EOF
    exit 0
fi

# ---------------------------------------------------------------------------- #
# demultiplex
# ---------------------------------------------------------------------------- #
LOGDIR="$OUTDIR/logs"
mkdir -p "$LOGDIR"
COUNTS="$OUTDIR/demux_counts.tsv"

i=0
for entry in "${PAIRS[@]}"; do
    i=$((i + 1))
    IFS=$'\t' read -r SID R1 R2 <<< "$entry"
    if [[ $SWAP_READS -eq 1 ]]; then TMP_R="$R1"; R1="$R2"; R2="$TMP_R"; fi

    DONE_MARK="$LOGDIR/$SID.done"
    if [[ -f "$DONE_MARK" && $FORCE -eq 0 ]]; then
        echo "[$i/${#PAIRS[@]}] $SID -- already done, skipping (--force to redo)"
    else
        echo "[$i/${#PAIRS[@]}] $SID"
        cutadapt \
            -g "file:$FWD_FASTA" \
            -G "file:$REV_FASTA" \
            --pair-adapters \
            -e "$ERROR_RATE" \
            -O "$MIN_OVERLAP" \
            --minimum-length "$MIN_LENGTH" \
            --cores "$CORES" \
            --json "$LOGDIR/$SID.cutadapt.json" \
            -o "$OUTDIR/{name}/${SID}_R1.fastq.gz" \
            -p "$OUTDIR/{name}/${SID}_R2.fastq.gz" \
            "$R1" "$R2" \
            > "$LOGDIR/$SID.cutadapt.log" 2>&1 \
          || { echo "ERROR: cutadapt failed for $SID; see $LOGDIR/$SID.cutadapt.log" >&2; exit 1; }
        touch "$DONE_MARK"
    fi
done

# ---------------------------------------------------------------------------- #
# 3' readthrough trim, for regions where the amplicon can be shorter than a read
# ---------------------------------------------------------------------------- #
RT="$PROJECT_DIR/config/readthrough.tsv"
while IFS=$'\t' read -r region flag ad_r1 ad_r2; do
    [[ "$region" == "region" ]] && continue
    [[ "$flag" == "yes" ]] || continue
    echo
    echo "3' readthrough trim: $region"
    for f1 in "$OUTDIR/$region"/*_R1.fastq.gz; do
        [[ -e "$f1" ]] || continue
        f2="${f1%_R1.fastq.gz}_R2.fastq.gz"
        sid="$(basename "${f1%_R1.fastq.gz}")"
        [[ -f "$OUTDIR/$region/.rt_$sid.done" && $FORCE -eq 0 ]] && continue
        # temp names must keep the .fastq.gz suffix: cutadapt chooses its output
        # compression from the extension, so a bare .tmp writes plain text and the
        # following mv would leave an uncompressed file under a .gz name
        t1="${f1%.fastq.gz}.rt_tmp.fastq.gz"
        t2="${f2%.fastq.gz}.rt_tmp.fastq.gz"
        cutadapt \
            -a "$ad_r1" -A "$ad_r2" \
            -e "$ERROR_RATE" --minimum-length "$MIN_LENGTH" --cores "$CORES" \
            -o "$t1" -p "$t2" "$f1" "$f2" \
            >> "$LOGDIR/${region}_readthrough.log" 2>&1
        mv "$t1" "$f1"
        mv "$t2" "$f2"
        touch "$OUTDIR/$region/.rt_$sid.done"
    done
done < "$RT"

# ---------------------------------------------------------------------------- #
# count read pairs -- after readthrough trimming, so the numbers match the files
# that step 02 and step 03 will actually consume
# ---------------------------------------------------------------------------- #
echo
echo "Counting read pairs per sample per region..."
printf 'sample\tregion\tread_pairs\n' > "$COUNTS"
for entry in "${PAIRS[@]}"; do
    IFS=$'\t' read -r SID _ _ <<< "$entry"
    for r in "${REGIONS[@]}" unknown; do
        n="$(count_pairs "$OUTDIR/$r/${SID}_R1.fastq.gz")"
        printf '%s\t%s\t%s\n' "$SID" "$r" "$n" >> "$COUNTS"
    done
done

# ---------------------------------------------------------------------------- #
# summary
# ---------------------------------------------------------------------------- #
echo
echo "=== assignment summary (all samples pooled) ==="
python3 - "$COUNTS" <<'PY'
import sys, collections
path = sys.argv[1]
tot = collections.OrderedDict()
per_sample = collections.defaultdict(int)
with open(path) as fh:
    next(fh)
    for line in fh:
        s, r, n = line.rstrip("\n").split("\t")
        tot[r] = tot.get(r, 0) + int(n)
        per_sample[s] += int(n)
grand = sum(tot.values())
print(f"{'region':10} {'read_pairs':>14} {'pct':>7}")
for r, n in tot.items():
    pct = 100.0 * n / grand if grand else 0.0
    print(f"{r:10} {n:14,} {pct:6.2f}%")
print(f"{'TOTAL':10} {grand:14,}")
unk = 100.0 * tot.get('unknown', 0) / grand if grand else 0.0
print()
if unk > 20:
    print(f"WARNING: {unk:.1f}% of pairs are unassigned. Try --probe to test read "
          f"orientation, or relax --error-rate / --min-overlap.")
else:
    print(f"unassigned: {unk:.1f}% (typically <20% for a healthy run)")
PY

echo
echo "Per-sample per-region counts: $COUNTS"
echo "Next: scripts/02_make_samplesheets.py --demux-dir $OUTDIR"
