#!/usr/bin/env bash
#
# Step 03 -- run nf-core/ampliseq once per primer region.
#
# Each region is an independent amplicon with its own length, its own optimal
# DADA2 truncation and (for ITS) its own reference database, so each gets its own
# run, work directory and results directory. Primers were already removed in
# step 01, hence --skip_cutadapt.
#
# Usage:
#   scripts/03_run_ampliseq.sh --all                 # every region, sequentially
#   scripts/03_run_ampliseq.sh V3V4                  # one region
#   scripts/03_run_ampliseq.sh V3V4 V4V5             # a subset
#   scripts/03_run_ampliseq.sh --all --dry-run       # print commands only
#
# Options:
#   --all               run every region in config/regions.tsv
#   --config FILE       region table               [config/regions.tsv]
#   --samplesheet-dir D from step 02               [samplesheets]
#   --outdir DIR        results root               [results_by_region]
#   --workdir DIR       nextflow work root         [work]
#   --profile NAME      nextflow -profile          [singularity]
#   --revision REV      ampliseq version, -r       [2.14.0]
#   --params-file FILE  extra ampliseq params yml  [config/ampliseq_base.yml if present]
#   --auto-trunc        ignore trunclenf/trunclenr, use quality-based truncation
#   --extra "ARGS"      extra args appended to every region's command
#   --no-resume         do not pass -resume
#   --dry-run           print the commands without running them
#   --continue-on-error keep going if one region fails
#
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"
LIB="$SCRIPT_DIR/lib/regions.py"

CONFIG="$PROJECT_DIR/config/regions.tsv"
SHEETDIR="$PROJECT_DIR/samplesheets"
OUTROOT="$PROJECT_DIR/results_by_region"
WORKROOT="$PROJECT_DIR/work"
PROFILE="singularity"
REVISION="2.14.0"
PARAMS_FILE=""
AUTO_TRUNC=0
GLOBAL_EXTRA=""
RESUME="-resume"
DRY_RUN=0
CONTINUE_ON_ERROR=0
REQUESTED=()
RUN_ALL=0

while [[ $# -gt 0 ]]; do
    case "$1" in
        --all)               RUN_ALL=1; shift ;;
        --config)            CONFIG="$2"; shift 2 ;;
        --samplesheet-dir)   SHEETDIR="$2"; shift 2 ;;
        --outdir)            OUTROOT="$2"; shift 2 ;;
        --workdir)           WORKROOT="$2"; shift 2 ;;
        --profile)           PROFILE="$2"; shift 2 ;;
        --revision)          REVISION="$2"; shift 2 ;;
        --params-file)       PARAMS_FILE="$2"; shift 2 ;;
        --auto-trunc)        AUTO_TRUNC=1; shift ;;
        --extra)             GLOBAL_EXTRA="$2"; shift 2 ;;
        --no-resume)         RESUME=""; shift ;;
        --dry-run)           DRY_RUN=1; shift ;;
        --continue-on-error) CONTINUE_ON_ERROR=1; shift ;;
        -h|--help)           sed -n '2,32p' "${BASH_SOURCE[0]}"; exit 0 ;;
        -*) echo "Unknown option: $1" >&2; exit 2 ;;
        *)  REQUESTED+=("$1"); shift ;;
    esac
done

[[ -f "$CONFIG" ]] || { echo "ERROR: config not found: $CONFIG" >&2; exit 2; }
command -v nextflow >/dev/null || { echo "ERROR: 'nextflow' not on PATH" >&2; exit 127; }

if [[ -z "$PARAMS_FILE" && -f "$PROJECT_DIR/config/ampliseq_base.yml" ]]; then
    PARAMS_FILE="$PROJECT_DIR/config/ampliseq_base.yml"
fi

# portable stand-in for bash 4's mapfile (macOS ships bash 3.2)
ALL_REGIONS=()
while IFS= read -r line; do [[ -n "$line" ]] && ALL_REGIONS+=("$line"); done \
    < <(python3 "$LIB" regions "$CONFIG")
if [[ $RUN_ALL -eq 1 ]]; then
    REGIONS=("${ALL_REGIONS[@]}")
elif [[ ${#REQUESTED[@]} -gt 0 ]]; then
    REGIONS=("${REQUESTED[@]}")
else
    echo "ERROR: name at least one region, or pass --all. Available: ${ALL_REGIONS[*]}" >&2
    exit 2
fi

# Singularity image cache -- shared across regions so images are pulled once.
if [[ "$PROFILE" == *singularity* ]]; then
    : "${NXF_SINGULARITY_CACHEDIR:=$PROJECT_DIR/singularity}"
    export NXF_SINGULARITY_CACHEDIR
    mkdir -p "$NXF_SINGULARITY_CACHEDIR"
    echo "NXF_SINGULARITY_CACHEDIR=$NXF_SINGULARITY_CACHEDIR"
fi

mkdir -p "$OUTROOT" "$WORKROOT" "$PROJECT_DIR/logs"

FAILED=()
for REGION in "${REGIONS[@]}"; do
    # validate against the config so a typo fails loudly
    found=0
    for r in "${ALL_REGIONS[@]}"; do [[ "$r" == "$REGION" ]] && found=1; done
    if [[ $found -eq 0 ]]; then
        echo "ERROR: '$REGION' is not in $CONFIG (have: ${ALL_REGIONS[*]})" >&2
        exit 2
    fi

    SHEET="$SHEETDIR/$REGION.tsv"
    if [[ ! -s "$SHEET" ]]; then
        echo "SKIP $REGION: samplesheet not found: $SHEET" >&2
        continue
    fi
    # header only == no samples survived step 02's --min-reads filter
    if [[ "$(wc -l < "$SHEET")" -le 1 ]]; then
        echo "SKIP $REGION: samplesheet has no samples: $SHEET" >&2
        continue
    fi

    TRUNCF="$(python3 "$LIB" field "$CONFIG" "$REGION" trunclenf)"
    TRUNCR="$(python3 "$LIB" field "$CONFIG" "$REGION" trunclenr)"
    REFTAX="$(python3 "$LIB" field "$CONFIG" "$REGION" ref_taxonomy)"
    EXTRA="$(python3 "$LIB" field "$CONFIG" "$REGION" extra_args)"
    [[ "$EXTRA" == "-" ]] && EXTRA=""

    TRUNC_ARGS=()
    if [[ $AUTO_TRUNC -eq 1 || "$TRUNCF" == "auto" || "$TRUNCR" == "auto" ]]; then
        # ampliseq's own quality-based truncation: keep the position where median
        # quality drops below trunc_qmin, provided trunc_rmin of reads survive
        TRUNC_ARGS=(--trunc_qmin 25 --trunc_rmin 0.75)
        TRUNC_DESC="auto (--trunc_qmin 25 --trunc_rmin 0.75)"
    else
        TRUNC_ARGS=(--trunclenf "$TRUNCF" --trunclenr "$TRUNCR")
        TRUNC_DESC="fixed ($TRUNCF / $TRUNCR)"
    fi

    N_SAMPLES=$(( $(wc -l < "$SHEET") - 1 ))
    echo
    echo "=============================================================="
    echo " region      : $REGION"
    echo " samples     : $N_SAMPLES"
    echo " samplesheet : $SHEET"
    echo " truncation  : $TRUNC_DESC"
    echo " taxonomy    : $REFTAX"
    echo " outdir      : $OUTROOT/$REGION"
    echo "=============================================================="

    CMD=(nextflow run nf-core/ampliseq
         -r "$REVISION"
         -profile "$PROFILE"
         -work-dir "$WORKROOT/$REGION"
         -ansi-log false)
    [[ -n "$RESUME" ]] && CMD+=("$RESUME")
    [[ -n "$PARAMS_FILE" ]] && CMD+=(-params-file "$PARAMS_FILE")
    CMD+=(--input "$SHEET"
          --outdir "$OUTROOT/$REGION"
          --skip_cutadapt
          --dada_ref_taxonomy "$REFTAX"
          "${TRUNC_ARGS[@]}")
    # shellcheck disable=SC2206
    [[ -n "$EXTRA" ]] && CMD+=($EXTRA)
    # shellcheck disable=SC2206
    [[ -n "$GLOBAL_EXTRA" ]] && CMD+=($GLOBAL_EXTRA)

    if [[ $DRY_RUN -eq 1 ]]; then
        printf '%q ' "${CMD[@]}"; echo
        continue
    fi

    LOG="$PROJECT_DIR/logs/ampliseq_${REGION}.log"
    if "${CMD[@]}" 2>&1 | tee "$LOG"; then
        echo "OK: $REGION"
    else
        echo "FAILED: $REGION (log: $LOG)" >&2
        FAILED+=("$REGION")
        [[ $CONTINUE_ON_ERROR -eq 1 ]] || exit 1
    fi
done

echo
if [[ ${#FAILED[@]} -gt 0 ]]; then
    echo "Regions that failed: ${FAILED[*]}"
    echo "Fix and rerun those regions; -resume will reuse completed work."
    exit 1
fi
[[ $DRY_RUN -eq 1 ]] || echo "Next: python3 scripts/04_compile_results.py --results-dir $OUTROOT"
