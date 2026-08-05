#!/usr/bin/env bash
#
# Run the whole per-region 16S/ITS workflow: demux -> samplesheets -> ampliseq -> compile.
#
# Usage:
#   scripts/run_all.sh --indir /path/to/fastq [--cores 16] [--min-reads 1000]
#
# Everything else is configured in config/regions.tsv and config/ampliseq_base.yml.
# Each step is restartable on its own, so if this stops partway you can rerun the
# individual script rather than starting over.
#
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

INDIR=""
CORES=8
MIN_READS=1000
EXTRA_AMPLISEQ=()

while [[ $# -gt 0 ]]; do
    case "$1" in
        --indir)     INDIR="$2"; shift 2 ;;
        --cores)     CORES="$2"; shift 2 ;;
        --min-reads) MIN_READS="$2"; shift 2 ;;
        --extra)     EXTRA_AMPLISEQ=(--extra "$2"); shift 2 ;;
        -h|--help)   sed -n '2,12p' "${BASH_SOURCE[0]}"; exit 0 ;;
        *) echo "Unknown option: $1" >&2; exit 2 ;;
    esac
done
[[ -n "$INDIR" ]] || { echo "ERROR: --indir is required" >&2; exit 2; }

echo "########## 1/4  demultiplex by primer region ##########"
bash "$SCRIPT_DIR/01_demux_regions.sh" --indir "$INDIR" --cores "$CORES"

echo
echo "########## 2/4  per-region samplesheets ##########"
python3 "$SCRIPT_DIR/02_make_samplesheets.py" --min-reads "$MIN_READS"

echo
echo "########## 3/4  ampliseq per region ##########"
bash "$SCRIPT_DIR/03_run_ampliseq.sh" --all --continue-on-error "${EXTRA_AMPLISEQ[@]+"${EXTRA_AMPLISEQ[@]}"}"

echo
echo "########## 4/4  compile across regions ##########"
python3 "$SCRIPT_DIR/04_compile_results.py"
