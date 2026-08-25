#!/usr/bin/env bash
#
# Step 00 -- check (and optionally install) the tools this pipeline needs.
#
# Default is CHECK ONLY: nothing is installed, nothing is changed. It prints what
# is present, what is missing, what is too old, and how to fix each one.
#
#   scripts/00_check_env.sh                 # report only
#   scripts/00_check_env.sh --install       # install what is missing, user-local
#   scripts/00_check_env.sh --install --dry-run
#
# Options:
#   --install          install missing tools (never uses sudo, never touches system dirs)
#   --method M         conda | pip | auto   how to install cutadapt/nextflow [auto]
#   --prefix DIR       install location                      [$HOME/.local/16s-tools]
#   --dry-run          with --install, print the commands instead of running them
#   --no-color         plain output (for logs)
#
# Requirements checked:
#   python3   >= 3.6      steps 01, 02, 04   (standard library only)
#   cutadapt  >= 3.4      step 01            (needs --pair-adapters)
#   nextflow  >= 24.04.2  step 03            (hard requirement of ampliseq 2.14.0)
#   java      17-22       step 03            (required by Nextflow 24.x; its launcher
#                                             rejects anything newer than 22, so an
#                                             unpinned "latest" JDK install breaks it)
#   container engine      step 03            singularity OR apptainer OR docker
#   gzip, awk, curl, tar  steps 01, 04
#
set -uo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"

DO_INSTALL=0
METHOD="auto"
PREFIX="$HOME/.local/16s-tools"
DRY_RUN=0
USE_COLOR=1

while [[ $# -gt 0 ]]; do
    case "$1" in
        --install)  DO_INSTALL=1; shift ;;
        --method)   METHOD="$2"; shift 2 ;;
        --prefix)   PREFIX="$2"; shift 2 ;;
        --dry-run)  DRY_RUN=1; shift ;;
        --no-color) USE_COLOR=0; shift ;;
        -h|--help)  sed -n '2,28p' "${BASH_SOURCE[0]}"; exit 0 ;;
        *) echo "Unknown option: $1" >&2; exit 2 ;;
    esac
done

case "$METHOD" in auto|conda|pip) ;; *) echo "ERROR: --method must be auto, conda or pip" >&2; exit 2 ;; esac

if [[ $USE_COLOR -eq 1 && -t 1 ]]; then
    G="\033[32m"; Y="\033[33m"; R="\033[31m"; B="\033[1m"; N="\033[0m"
else
    G=""; Y=""; R=""; B=""; N=""
fi

# minimum versions
MIN_PYTHON=3.6
MIN_CUTADAPT=3.4
MIN_NEXTFLOW=24.04.2
MIN_JAVA=17
MAX_JAVA=22  # Nextflow's launcher hard-rejects anything newer than this

# ---------------------------------------------------------------------------- #
# helpers
# ---------------------------------------------------------------------------- #

# first dotted number in stdin, e.g. "Docker version 24.0.6, build x" -> 24.0.6
first_version() { grep -oE '[0-9]+(\.[0-9]+)+' | head -1; }

# version_ge A B -> 0 if A >= B. Compares field by field numerically, so
# "1.8.0_501" < "17" and "24.10.4.5934" > "24.04.2" both come out right.
version_ge() {
    awk -v A="$1" -v B="$2" 'BEGIN{
        gsub(/^[^0-9]+/,"",A); gsub(/^[^0-9]+/,"",B);
        na=split(A,a,/[^0-9]+/); nb=split(B,b,/[^0-9]+/);
        n=(na>nb?na:nb);
        for(i=1;i<=n;i++){
            x=(i<=na&&a[i]!="")?a[i]+0:0; y=(i<=nb&&b[i]!="")?b[i]+0:0;
            if(x>y) exit 0; if(x<y) exit 1;
        }
        exit 0
    }'
}

REPORT=()      # "status|tool|found|required|note"
N_MISSING=0
N_OLD=0
MISSING_TOOLS=()

record() { REPORT+=("$1|$2|$3|$4|$5"); }

# check_tool <name> <min_version|-> <version command...>
check_tool() {
    local name="$1" min="$2"; shift 2
    local path found
    path="$(command -v "$name" 2>/dev/null)"
    if [[ -z "$path" ]]; then
        record MISSING "$name" "-" "$min" "not on PATH"
        N_MISSING=$((N_MISSING+1)); MISSING_TOOLS+=("$name")
        return 1
    fi
    found="$("$@" 2>&1 | first_version)"
    if [[ -z "$found" ]]; then
        record OK "$name" "present" "$min" "$path"
        return 0
    fi
    if [[ "$min" != "-" ]] && ! version_ge "$found" "$min"; then
        record OLD "$name" "$found" "$min" "$path"
        N_OLD=$((N_OLD+1)); MISSING_TOOLS+=("$name")
        return 1
    fi
    record OK "$name" "$found" "$min" "$path"
    return 0
}

# java needs both a floor and a ceiling: Nextflow's launcher script hard-rejects
# any Java newer than $MAX_JAVA, so a plain ">= MIN_JAVA" install can pick the
# latest JDK from conda/apt and silently break step 03.
check_java() {
    local path found
    path="$(command -v java 2>/dev/null)"
    if [[ -z "$path" ]]; then
        record MISSING java "-" "$MIN_JAVA-$MAX_JAVA" "not on PATH"
        N_MISSING=$((N_MISSING+1)); MISSING_TOOLS+=("java")
        return 1
    fi
    found="$(java -version 2>&1 | first_version)"
    if [[ -z "$found" ]]; then
        record OK java "present" "$MIN_JAVA-$MAX_JAVA" "$path"
        return 0
    fi
    if ! version_ge "$found" "$MIN_JAVA"; then
        record OLD java "$found" "$MIN_JAVA-$MAX_JAVA" "$path (too old)"
        N_OLD=$((N_OLD+1)); MISSING_TOOLS+=("java")
        return 1
    fi
    if ! version_ge "$MAX_JAVA" "$found"; then
        record OLD java "$found" "$MIN_JAVA-$MAX_JAVA" "$path (too new -- Nextflow's launcher caps at Java $MAX_JAVA)"
        N_OLD=$((N_OLD+1)); MISSING_TOOLS+=("java")
        return 1
    fi
    record OK java "$found" "$MIN_JAVA-$MAX_JAVA" "$path"
    return 0
}

# ---------------------------------------------------------------------------- #
# run the checks
# ---------------------------------------------------------------------------- #
echo
echo -e "${B}Environment check for the per-region 16S/ITS pipeline${N}"
echo "project: $PROJECT_DIR"
echo "host:    $(hostname 2>/dev/null || echo unknown)  ($(uname -sr 2>/dev/null))"
echo

check_tool python3  "$MIN_PYTHON"   python3 --version
check_tool cutadapt "$MIN_CUTADAPT" cutadapt --version
check_tool nextflow "$MIN_NEXTFLOW" nextflow -v
check_java

# container engine: any one of the three satisfies step 03
ENGINE=""
for e in singularity apptainer docker; do
    if command -v "$e" >/dev/null 2>&1; then ENGINE="$e"; break; fi
done
if [[ -n "$ENGINE" ]]; then
    ev="$("$ENGINE" --version 2>&1 | first_version)"
    record OK "container" "${ENGINE} ${ev:-present}" "any" "$(command -v "$ENGINE")"
else
    record MISSING "container" "-" "any" "need singularity, apptainer or docker"
    N_MISSING=$((N_MISSING+1)); MISSING_TOOLS+=("container")
fi

for t in gzip awk curl tar; do
    if command -v "$t" >/dev/null 2>&1; then
        record OK "$t" "present" "-" "$(command -v "$t")"
    else
        record MISSING "$t" "-" "-" "core utility, expected on any server"
        N_MISSING=$((N_MISSING+1)); MISSING_TOOLS+=("$t")
    fi
done

# ---------------------------------------------------------------------------- #
# report
# ---------------------------------------------------------------------------- #
printf "%-11s %-9s %-22s %-11s %s\n" "TOOL" "STATUS" "FOUND" "REQUIRED" "WHERE / NOTE"
printf '%s\n' "-------------------------------------------------------------------------------"
for row in "${REPORT[@]}"; do
    IFS='|' read -r st tool found req note <<< "$row"
    case "$st" in
        OK)      col="$G" ;;
        OLD)     col="$Y" ;;
        MISSING) col="$R" ;;
        *)       col="" ;;
    esac
    printf "%-11s ${col}%-9s${N} %-22s %-11s %s\n" "$tool" "$st" "$found" "$req" "$note"
done
echo

# ---------------------------------------------------------------------------- #
# environment sanity beyond the tools themselves
# ---------------------------------------------------------------------------- #
echo -e "${B}Environment${N}"

if [[ -n "${NXF_SINGULARITY_CACHEDIR:-}" ]]; then
    if [[ -d "$NXF_SINGULARITY_CACHEDIR" && -w "$NXF_SINGULARITY_CACHEDIR" ]]; then
        echo -e "  ${G}OK${N}      NXF_SINGULARITY_CACHEDIR=$NXF_SINGULARITY_CACHEDIR (writable)"
    else
        echo -e "  ${Y}WARN${N}    NXF_SINGULARITY_CACHEDIR=$NXF_SINGULARITY_CACHEDIR is missing or not writable"
    fi
else
    echo -e "  ${Y}WARN${N}    NXF_SINGULARITY_CACHEDIR is not set -- images will be re-downloaded per run"
    echo "          export NXF_SINGULARITY_CACHEDIR=/tank/lmapunda/ARI/work/singularity/"
fi

AVAIL="$(df -Pk "$PROJECT_DIR" 2>/dev/null | awk 'NR==2{printf "%.0f", $4/1048576}')"
if [[ -n "$AVAIL" ]]; then
    if [[ "$AVAIL" -lt 100 ]]; then
        echo -e "  ${Y}WARN${N}    only ${AVAIL} GB free here; 7 regions of work dirs plus reference DBs is easily >100 GB"
    else
        echo -e "  ${G}OK${N}      ${AVAIL} GB free on this filesystem"
    fi
fi

# HPC module systems are the usual way singularity and java are provided
if command -v module >/dev/null 2>&1 || [[ -n "${MODULEPATH:-}" ]]; then
    echo -e "  ${G}INFO${N}    a module system is available -- these may already be provided:"
    echo "          module avail 2>&1 | grep -iE 'singularity|apptainer|java|jdk|nextflow'"
fi
echo

# ---------------------------------------------------------------------------- #
# verdict
# ---------------------------------------------------------------------------- #
STEP123_OK=1
for t in "${MISSING_TOOLS[@]:-}"; do
    [[ -z "$t" ]] && continue
    STEP123_OK=0
done

if [[ $N_MISSING -eq 0 && $N_OLD -eq 0 ]]; then
    echo -e "${G}All requirements satisfied.${N} You can start with:"
    echo "    scripts/01_demux_regions.sh --indir \$FASTQ --probe"
    echo
    exit 0
fi

echo -e "${B}Problems:${N} $N_MISSING missing, $N_OLD too old -> ${MISSING_TOOLS[*]}"
echo

# Which steps can still run? Steps 01/02/04 only need python3 + cutadapt.
can_run_demux=1
for t in "${MISSING_TOOLS[@]:-}"; do
    [[ "$t" == "python3" || "$t" == "cutadapt" ]] && can_run_demux=0
done
if [[ $can_run_demux -eq 1 ]]; then
    echo "Note: steps 01, 02 and 04 only need python3 + cutadapt, and those are fine."
    echo "      You can demultiplex now and sort out the Nextflow stack before step 03."
    echo
fi

if [[ $DO_INSTALL -eq 0 ]]; then
    cat <<'EOF'
To install the missing pieces, user-local and without sudo:

    scripts/00_check_env.sh --install

Preview first with --install --dry-run. What it will and will not do:

  cutadapt, nextflow, java  installed into --prefix (conda env if conda/mamba is
                            available, otherwise a Python venv + the official
                            Nextflow installer)
  singularity / apptainer   NOT installed -- these need root. On a shared server
                            they are normally provided by the admins:
                                module load singularity   (or apptainer)
                            If neither exists, ask your admin, or use
                            -profile docker if you are in the docker group.
EOF
    echo
    exit 1
fi

# ---------------------------------------------------------------------------- #
# --install
# ---------------------------------------------------------------------------- #
echo -e "${B}Installing${N} (prefix: $PREFIX)"
[[ $DRY_RUN -eq 1 ]] && echo "(dry run -- commands are printed, not executed)"
echo

INSTALL_FAILED=0

run() {
    # printf %q so specs like cutadapt>=3.4 display quoted and stay copy-pasteable
    printf '    $ '; printf '%q ' "$@"; printf '\n'
    if [[ $DRY_RUN -eq 0 ]]; then
        "$@" || { echo -e "  ${R}FAILED${N}: $*" >&2; INSTALL_FAILED=1; return 1; }
    fi
}

needs() {
    local want="$1"
    for t in "${MISSING_TOOLS[@]:-}"; do [[ "$t" == "$want" ]] && return 0; done
    return 1
}

CONDA_BIN=""
for c in mamba micromamba conda; do
    if command -v "$c" >/dev/null 2>&1; then CONDA_BIN="$c"; break; fi
done

if [[ "$METHOD" == "auto" ]]; then
    if [[ -n "$CONDA_BIN" ]]; then METHOD="conda"; else METHOD="pip"; fi
fi
echo "  method: $METHOD${CONDA_BIN:+ (found $CONDA_BIN)}"
echo

ENV_DIR="$PREFIX/env"
BIN_DIR="$PREFIX/bin"

if [[ "$METHOD" == "conda" ]]; then
    if [[ -z "$CONDA_BIN" ]]; then
        echo -e "  ${R}ERROR${N}: --method conda but no conda/mamba/micromamba on PATH." >&2
        echo "  Use --method pip, or install Miniforge first:" >&2
        echo "    curl -L -o Miniforge3.sh https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-\$(uname)-\$(uname -m).sh" >&2
        echo "    bash Miniforge3.sh -b -p \$HOME/miniforge3 && \$HOME/miniforge3/bin/conda init bash" >&2
        exit 1
    fi
    PKGS=()
    needs cutadapt && PKGS+=("cutadapt>=$MIN_CUTADAPT")
    needs nextflow && PKGS+=("nextflow>=$MIN_NEXTFLOW")
    needs java     && PKGS+=("openjdk>=$MIN_JAVA,<=$MAX_JAVA")
    # nextflow from conda pulls its own JDK, but be explicit when java is the problem
    if [[ ${#PKGS[@]} -eq 0 ]]; then
        echo "  nothing conda can fix here (see the notes above)."
    else
        run "$CONDA_BIN" create -y -p "$ENV_DIR" -c conda-forge -c bioconda "${PKGS[@]}"
    fi
    ACTIVATE="$ENV_DIR/bin"
else
    [[ $DRY_RUN -eq 0 ]] && mkdir -p "$BIN_DIR"
    if needs cutadapt; then
        echo "  cutadapt -> Python venv at $PREFIX/venv"
        run python3 -m venv "$PREFIX/venv"
        run "$PREFIX/venv/bin/pip" install --quiet --upgrade pip
        run "$PREFIX/venv/bin/pip" install --quiet "cutadapt>=$MIN_CUTADAPT"
        if [[ $DRY_RUN -eq 0 && -x "$PREFIX/venv/bin/cutadapt" ]]; then
            ln -sf "$PREFIX/venv/bin/cutadapt" "$BIN_DIR/cutadapt"
        fi
    fi
    if needs nextflow; then
        echo "  nextflow -> official installer into $BIN_DIR"
        if [[ $DRY_RUN -eq 0 ]]; then
            ( cd "$BIN_DIR" && curl -fsSL https://get.nextflow.io | bash ) \
                && chmod +x "$BIN_DIR/nextflow" \
                || { echo -e "  ${R}FAILED${N}: nextflow install" >&2; INSTALL_FAILED=1; }
        else
            echo "    \$ (cd $BIN_DIR && curl -fsSL https://get.nextflow.io | bash)"
        fi
    fi
    if needs java; then
        echo -e "  ${Y}java${N}: not installed by --method pip. Options:"
        echo "        module load java        (or jdk / openjdk)"
        echo "        scripts/00_check_env.sh --install --method conda"
    fi
    ACTIVATE="$BIN_DIR"
fi

if needs container; then
    echo
    echo -e "  ${Y}container engine${N}: not installed -- singularity/apptainer need root."
    echo "        try:  module load singularity      (or: module load apptainer)"
    echo "        then re-run this check."
fi

# ---------------------------------------------------------------------------- #
# write an env file to source
# ---------------------------------------------------------------------------- #
ENV_FILE="$PROJECT_DIR/config/env.sh"
if [[ $DRY_RUN -eq 0 ]]; then
    {
        echo "# Generated by scripts/00_check_env.sh -- source this before running the pipeline."
        echo "#   source config/env.sh"
        echo "export PATH=\"$ACTIVATE:\$PATH\""
        if [[ -n "${NXF_SINGULARITY_CACHEDIR:-}" ]]; then
            echo "export NXF_SINGULARITY_CACHEDIR=\"$NXF_SINGULARITY_CACHEDIR\""
        else
            echo "# set this to a shared cache so images are downloaded once:"
            echo "# export NXF_SINGULARITY_CACHEDIR=/tank/lmapunda/ARI/work/singularity/"
        fi
    } > "$ENV_FILE"
    echo
    echo "Wrote $ENV_FILE"
fi

echo
if [[ $DRY_RUN -eq 1 ]]; then
    echo -e "${B}Dry run -- nothing was installed or written.${N}"
    echo "Re-run without --dry-run to apply:"
    echo "    scripts/00_check_env.sh --install --method $METHOD"
else
    echo -e "${B}Next:${N}"
    echo "    source config/env.sh"
    echo "    scripts/00_check_env.sh          # confirm everything is now OK"
fi

# non-zero when anything failed or is still outstanding, so this is safe to
# chain with && in a setup script
if [[ $INSTALL_FAILED -eq 1 ]]; then
    echo
    echo -e "${R}Some installs failed${N} -- see the messages above. Fix those, then re-run the check."
    exit 1
fi
if needs container || needs java; then
    exit 1
fi
exit 0
