# Job aid — running the per-region 16S/ITS pipeline

One page, in order. Full detail and troubleshooting live in [README.md](README.md).

**The flow:** `check env → probe → demux → samplesheets → ampliseq → compile`

---

## Before you start (once per run)

```bash
# 1. log in to the analysis server
ssh lmapunda@<analysis-server> -p <port>

# 2. go to the run directory
cd /tank/lmapunda/ARI/16s/Run_20072026

# 3. point Singularity at the shared image cache (avoids re-downloading images)
export NXF_SINGULARITY_CACHEDIR=/tank/lmapunda/ARI/work/singularity/

# 4. check every required tool in one go
scripts/00_check_env.sh
```

If anything is missing or too old, install it user-local (no sudo):

```bash
scripts/00_check_env.sh --install --dry-run   # preview
scripts/00_check_env.sh --install             # do it
source config/env.sh                          # put it on PATH
scripts/00_check_env.sh                       # confirm
```

`source config/env.sh` is needed in **every new shell** (add it to `~/.bashrc` to
make it permanent). Singularity cannot be installed without root — if it is
missing, try `module load singularity` or ask your admin.

Set `FASTQ=` to your raw read directory so you can paste the rest verbatim:

```bash
FASTQ=/tank/lmapunda/ARI/16s/Run_20072026/work_file
```

> Run long steps under `tmux` or `screen` so a dropped connection doesn't kill them.

---

## Step 0 — Probe (2 minutes, do not skip)

Confirms the reads are the right way round and the primers match your kit.

```bash
scripts/01_demux_regions.sh --indir $FASTQ --probe
```

**Look at the `ASSIGNED` line for each orientation.**

| What you see | What to do |
|---|---|
| `as_given` high (>70%), `swapped` ~0% | Normal — continue to step 1 |
| `swapped` much higher | Add `--swap-reads` to every step 1 command |
| Both very low | **Stop.** Primers don't match this library — check the kit/panel |
| One region at 0%, others fine | That region wasn't in your panel — delete its row from `config/regions.tsv` |

---

## Step 1 — Demultiplex by region

```bash
scripts/01_demux_regions.sh --indir $FASTQ --cores 16
```

Rough guide: a few minutes per sample. Safe to re-run — finished samples are
skipped (`--force` to redo).

**Check:** the summary table printed at the end.

- `unassigned` under ~20% → good
- all 7 regions getting a similar share → good
- one region near zero → likely not in your panel

Counts are saved to `demux/demux_counts.tsv`.

---

## Step 2 — Build samplesheets

```bash
python3 scripts/02_make_samplesheets.py --min-reads 1000
```

Takes seconds.

**Check:** the `Runnable regions:` line lists all your regions, then open:

```bash
cat samplesheets/excluded_samples.tsv
```

Samples below 1,000 reads are dropped **for that region only**. Expect your NTC
and the ~6 near-empty libraries here — that's correct. If a *real* sample appears
across many regions, investigate it before continuing.

---

## Step 3 — Run ampliseq (the long one)

Preview the commands first:

```bash
scripts/03_run_ampliseq.sh --all --dry-run
```

Then run one region to make sure the setup works end to end:

```bash
scripts/03_run_ampliseq.sh V3V4
```

If that finishes, run the rest:

```bash
scripts/03_run_ampliseq.sh --all --continue-on-error
```

Rough guide: hours per region, sequential. `--continue-on-error` finishes the
other regions and reports failures at the end. Re-running resumes — completed
work is not repeated.

**Check per region:**

```bash
column -t -s$'\t' results_by_region/V3V4/overall_summary.tsv | less -S
```

- `retained_percent` mostly >50% → good
- `merged` collapsing vs `denoisedF` → truncation problem, see the fix below
- logs are in `logs/ampliseq_<REGION>.log`

---

## Step 4 — Compile

```bash
python3 scripts/04_compile_results.py
```

Takes under a minute. Then read the two files that matter:

```bash
column -t -s$'\t' compiled/region_qc_summary.tsv          # which regions worked
head -20 compiled/genus_matrix_across_regions.tsv         # the biology
```

In `genus_matrix_across_regions.tsv`, `n_regions_detected` is your confidence
column: taxa found by several regions are solid, taxa found by one need checking.

Everything else in `compiled/` is described in the README table.

---

## Or run all four at once

Only once you've done the probe and are happy with the config:

```bash
scripts/run_all.sh --indir $FASTQ --cores 16
```

---

## The one fix you're most likely to need

**Symptom:** a region's `merged` count is far below `denoisedF` (this is what went
wrong in the earlier pooled run).

**Cause:** truncation left too little overlap for DADA2 to join R1 and R2.

**Fix:** edit that region's row in `config/regions.tsv` — lower `trunclenr` by
20–30, or set both `trunclenf` and `trunclenr` to `auto` — then rerun only that
region:

```bash
scripts/03_run_ampliseq.sh V3V4
python3 scripts/04_compile_results.py
```

V3V4 has the least overlap headroom, so it's the usual suspect. The overlap table
in the README shows the arithmetic.

---

## Quick reference

| Need | Command |
|---|---|
| Check required tools | `scripts/00_check_env.sh` |
| Install missing tools | `scripts/00_check_env.sh --install` |
| Check orientation / primers | `scripts/01_demux_regions.sh --indir $FASTQ --probe` |
| Redo demux from scratch | add `--force` to step 1 |
| Keep low-read samples | `--min-reads 0` in step 2 |
| See commands without running | add `--dry-run` to step 3 |
| Rerun one region | `scripts/03_run_ampliseq.sh V3V4` |
| Rebuild metadata.tsv from the AFI/ARI .dta files | `python3 scripts/build_metadata.py` |
| Check metadata covers every sample before running diversity | `python3 scripts/check_metadata.py` |
| Add diversity analysis | `scripts/03_run_ampliseq.sh --all --extra "--metadata config/metadata.tsv"` |
| Alpha/beta diversity by specimen type and site (after step 4) | `python3 scripts/05_diversity_analysis.py` |
| Taxa barplots per region (after step 4) | `python3 scripts/06_taxa_barplots.py --regions V1V2 V2V3 ...` |
| Diversity boxplots per region (after step 5) | `python3 scripts/07_diversity_boxplots.py` |
| Per-sample organism export for biostatistics (after step 3) | `python3 scripts/08_export_for_biostatistician.py --regions V3V4 V4V5` |
| Bacteria only (drop host/mito) | add `--extra "--filter_ssu bac"` to step 3 |

**Edit `config/regions.tsv` only** — `config/primers_*.fasta` are regenerated
automatically on every step 1 run and any manual edits there will be overwritten.
