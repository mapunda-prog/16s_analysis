# QIAseq 16S/ITS — per-region ampliseq workflow

Splits a QIAseq 16S/ITS multiplex library into its 7 primer regions, runs
nf-core/ampliseq independently on each, and compiles the results back together.

```
scripts/00_check_env.sh         check/install required tools
scripts/01_demux_regions.sh     reads  ->  demux/<REGION>/
scripts/02_make_samplesheets.py        ->  samplesheets/<REGION>.tsv
scripts/03_run_ampliseq.sh             ->  results_by_region/<REGION>/
scripts/04_compile_results.py          ->  compiled/
scripts/run_all.sh                     all four in order

-- optional, local-analyst-side (see Requirements for extra deps) --
scripts/build_metadata.py       clinical .dta files  ->  config/metadata.tsv
scripts/check_metadata.py       config/metadata.tsv vs samplesheets, pre-flight check
scripts/05_diversity_analysis.py       compiled/ + metadata  ->  alpha_beta_diversity/
scripts/06_taxa_barplots.py            compiled/genus_counts_by_region.tsv  ->  taxa_barplots/
scripts/07_diversity_boxplots.py       alpha_beta_diversity/  ->  alpha_beta_diversity/boxplots/
scripts/08_export_for_biostatistician.py  per-region qiime2 rel-abundance tables  ->  organisms_by_sample/
scripts/09_merge_organisms_with_metadata.py  organisms_by_sample/ + config/metadata.tsv  ->  organisms_by_sample/*_with_metadata.tsv
scripts/10_samples_missing_results.py     metadata vs organisms_by_sample/  ->  console report
scripts/11_export_by_rank.py           per-region dada2 ASV+tax tables (incl. ITS)  ->  organisms_by_sample/{family,genus,species}_level.tsv
scripts/12_rarefaction_plots.py        per-region dada2 ASV+tax tables  ->  alpha_beta_diversity/rarefaction_curves/
```

---

## Why per-region runs are necessary

The QIAseq Screening Panel amplifies 7 targets (V1V2, V2V3, V3V4, V4V5, V5V7,
V7V9, ITS1) in one library, so every sample's FASTQ is a mixture of amplicons
ranging from ~250 to ~450 bp. DADA2 assumes one amplicon of one length: a single
`--trunclenf/--trunclenr` pair cannot leave a mergeable overlap for all 7 at once,
and the error model is fitted across incompatible amplicons.

That is visible in the earlier pooled run in `../results/`:

| sample | denoisedF | merged | retained |
|--------|-----------|--------|----------|
| AIC003 | 264,366   | 3,437  | 11.4%    |
| AIC013 | 133,908   | 1,234  | 7.0%     |

~1–3% of denoised reads merged. Inspecting `dada2/ASV_seqs.fasta` from that run
confirms the cause directly: sequences begin with the V3V4 forward primer
(`CCTACGGG…`) and end in the V5V7 forward primer (`GGGATTAGATACCCGGGTAGTC`) — the
primers were never removed (`--skip_cutadapt` with no prior trimming) and reads
from different regions were denoised together.

Splitting first fixes both problems: each region gets its own truncation, its own
error model, and its own reference database.

---

## Primer provenance

`config/regions.tsv` carries the primer sequences, and they are not guesses.

QIAGEN does not publish the panel primers in the handbook, but its official
demultiplexer image (`qiaseq/qiaseq-16s`, the container behind
[github.com/qiaseq/qiaseq-16S](https://github.com/qiaseq/qiaseq-16S)) ships the
primer table it uses internally. That file is preserved here at
`reference/16S_primers_ScreeningPanel_Demultiplex_Final.txt.expanded` — 1,760 rows
covering all 7 regions.

It lists every **phase variant** of every primer with degenerate bases
pre-expanded. QIAseq uses phased primers: 0–11 extra bases prepended to the 5' end
to stagger the reads and improve base diversity. Phase *N* is exactly the phase-0
primer with *N* bases in front:

```
phase  0  CCTACGGGNGGCWGCAG            <- 357F core
phase  1  GCCTACGGGNGGCWGCAG
phase  2  AGCCTACGGGNGGCWGCAG
...
phase 11  TAGCTAGTTAGCCTACGGGNGGCWGCAG
```

So `config/regions.tsv` stores the **phase-0 primer with degenerate expansions
collapsed back to IUPAC**, and step 01 matches it with an *unanchored* 5' search.
A 5' trim removes everything up to and including the match, which strips the
phasing bases and the primer together — one sequence per region handles all 12
phases, with no need to enumerate them.

The recovered primers match their published counterparts, which is a useful
independent check:

| Region | Forward | Reverse | Known as |
|--------|---------|---------|----------|
| V1V2 | `AGRGTTTGATYMTGGCTC` | `CTGCTGCCTYCCGTA` | 27F-mod / BSR357 |
| V2V3 | `GGCGNACGGGTGAGTAA` | `WTTACCGCGGCTGCTGG` | F101-mod / 518R |
| V3V4 | `CCTACGGGNGGCWGCAG` | `GACTACHVGGGTATCTAATCC` | 341F / 805R |
| V4V5 | `GTGYCAGCMGCCGCGGTAA` | `CCGYCAATTNMTTTRAGTTT` | 515F / 926R |
| V5V7 | `GGATTAGATACCCBRGTAGTC` | `ACGTCRTCCCCDCCTTCCTC` | 785F / 1175R |
| V7V9 | `YAACGAGCGMRACCC` | `TACGGYTACCTTGTTANGACTT` | P699D / 1492R |
| ITS  | `CTTGGTCATTTAGAGGAAGTAA` | `GCTGCGTTCTTCATCGATGC` | ITS1f / ITS2 |

To regenerate the reference file yourself:

```bash
docker create --name qs qiaseq/qiaseq-16s
docker cp qs:/home/qiagen/data/16S_primers_ScreeningPanel_Demultiplex_Final.txt.expanded .
docker rm qs
```

> If your kit is the **Region Panel** with a custom subset rather than the
> Screening Panel, verify the region list against your order and drop any rows
> from `config/regions.tsv` you did not amplify.

---

## Requirements

| Tool | Used by | Notes |
|------|---------|-------|
| `cutadapt` ≥ 3.4 | step 01 | needs `--pair-adapters`; developed against 5.2 |
| `python3` ≥ 3.6 | steps 01, 02, 04 | standard library only, no pandas |
| `nextflow` ≥ 24.04.2 | step 03 | hard requirement of ampliseq 2.14.0 |
| `java` 17+ | step 03 | required by Nextflow 24.x |
| Singularity / Docker | step 03 | container engine for `-profile` |

Steps 5–12 are optional, local-analyst-side scripts (run on your own machine
against a downloaded `results_by_region/`/`compiled/`, not on the server).
Most need one extra package not covered by `00_check_env.sh`: `pandas`
(`build_metadata.py`), `numpy` (`05_diversity_analysis.py`), `matplotlib`
(`06_taxa_barplots.py`, `07_diversity_boxplots.py`, `12_rarefaction_plots.py`).
Steps 8-11 (and `check_metadata.py`) are pure standard library.

Check all of it at once, and install what is missing without sudo:

```bash
scripts/00_check_env.sh                 # report only, changes nothing
scripts/00_check_env.sh --install       # user-local install, then: source config/env.sh
```

It verifies versions rather than mere presence (Nextflow 24.04.2 is a hard
requirement of ampliseq 2.14.0, and Nextflow 24.x itself needs Java 17+), checks
that `NXF_SINGULARITY_CACHEDIR` is set and writable, warns on low disk, and exits
non-zero if anything is unresolved so it can be chained in a setup script.
Singularity/apptainer need root and are reported, not installed — on a shared
server they normally come from `module load singularity`.

Steps 01/02/04 are plain bash + Python and run anywhere; only step 03 needs the
Nextflow stack, so a missing Nextflow does not block demultiplexing.

---

## Usage

Everything below assumes you are in this directory.

### 0. Sanity-check orientation and primer match first

Cheap, and it catches the two failure modes that would otherwise waste a full run:

```bash
scripts/01_demux_regions.sh --indir /path/to/fastq --probe
```

This demultiplexes 100,000 read pairs from one sample twice — as-given and with
R1/R2 swapped — and reports the assignment rate for each. If the swapped
orientation wins, add `--swap-reads` to the real run. If *neither* orientation
assigns much, stop: the primers in the config do not match your library.

### 1. Demultiplex by region

```bash
scripts/01_demux_regions.sh --indir /path/to/fastq --cores 16
```

One cutadapt pass per sample. A read pair is assigned to a region only when R1
carries that region's forward primer **and** R2 carries the matching reverse
primer (`--pair-adapters`); both are then trimmed. Pairs that fail go to
`demux/unknown/` rather than being discarded, so you can inspect them.

Writes `demux/<REGION>/<sample>_R{1,2}.fastq.gz`, per-sample cutadapt logs and
JSON reports under `demux/logs/`, and `demux/demux_counts.tsv`. Finished samples
are marked done and skipped on rerun unless you pass `--force`.

For regions flagged `trim_readthrough=yes` (ITS by default, since ITS1 amplicons
can be shorter than a read) a second pass removes the opposite primer from the 3'
end.

Useful knobs: `--error-rate 0.15`, `--min-overlap` (defaults to the shortest
primer length, 15), `--min-length 100`.

### 2. Build per-region samplesheets

```bash
python3 scripts/02_make_samplesheets.py --min-reads 1000
```

Writes `samplesheets/<REGION>.tsv` in ampliseq's native
`sampleID/forwardReads/reverseReads` layout with absolute paths. Sample IDs are
cleaned the same way as `../make_ampliseq_samplesheet.py`
(`MLTP-AIC-126-SW` → `AIC126SW`) so they stay comparable with the earlier run.

Samples below `--min-reads` are dropped **for that region only** and recorded in
`samplesheets/excluded_samples.tsv`. This is deliberate: DADA2 fits a per-sample
error model, and one near-empty sample can fail an entire region's run. Read the
exclusion file — it is also your negative-control check.

### 3. Run ampliseq per region

```bash
scripts/03_run_ampliseq.sh --all                # all regions, sequentially
scripts/03_run_ampliseq.sh V3V4                 # just one
scripts/03_run_ampliseq.sh --all --dry-run      # print the commands first
```

Each region gets its own `-work-dir` and `--outdir`, so regions never collide and
`-resume` works per region. Truncation, reference database and extra flags come
from that region's row in `config/regions.tsv`; shared settings live in
`config/ampliseq_base.yml`. `--skip_cutadapt` is always passed, since step 01
already removed the primers.

`--continue-on-error` finishes the remaining regions and reports failures at the
end instead of stopping at the first one.

To add metadata-driven diversity analysis, first regenerate `config/metadata.tsv`
from the clinical `.dta` files (needs pandas; run locally, then commit the
result) and confirm it actually covers every sample that's about to be run:

```bash
python3 scripts/build_metadata.py
python3 scripts/check_metadata.py    # fails if a real sample has no metadata row
```

`check_metadata.py` compares `config/metadata.tsv` against the per-region
samplesheets from step 2. QIIME2's diversity plugins silently drop or error on
any sample missing from the metadata file, so this is meant to catch a gap
before a multi-hour ampliseq run rather than after. A sample already listed in
`samplesheets/excluded_samples.tsv` (dropped for low reads) is not flagged.

Then:

```bash
scripts/03_run_ampliseq.sh --all --extra "--metadata config/metadata.tsv"
```

(and remove the `skip_alpha_rarefaction` / `skip_diversity_indices` lines from
`config/ampliseq_base.yml`).

### 4. Compile

```bash
python3 scripts/04_compile_results.py
```

| Output in `compiled/` | Contents |
|---|---|
| `region_qc_summary.tsv` | one row per region: samples, reads in/out, median % retained, ASV and genus counts |
| `read_tracking_by_region.tsv` | every region's `overall_summary.tsv`, stacked with a `region` column |
| `asv_long.tsv.gz` | tidy `region, sampleID, ASV_ID, count` + full taxonomy, nonzero counts only |
| `genus_counts_by_region.tsv` | genus-level counts per sample per region |
| `genus_matrix_across_regions.tsv` | genus × region totals plus `n_regions_detected` |
| `demux_summary.tsv` | step 01 assignment counts |

`genus_matrix_across_regions.tsv` is the one to look at first. Taxa detected by
several regions are solid; taxa seen by only one are either that region's genuine
extra resolution or an artifact, and the column tells you which to check.

ASVs classified as Mitochondria or Chloroplast (host/plant DNA) are excluded
from every table here, mirroring ampliseq's own `--exclude_taxa` default --
this matters most for blood specimens, where host mtDNA can otherwise
dominate a sample's reads entirely.

### 5. Alpha/beta diversity by specimen type and collection site

```bash
python3 scripts/05_diversity_analysis.py --asv-long compiled/asv_long.tsv.gz
```

Needs numpy (run locally, like `build_metadata.py`). Reads `compiled/asv_long.tsv.gz`
and `config/metadata.tsv`, drops samples below 1,000 reads per region, rarefies
the rest to that region's lowest surviving depth, and writes to
`alpha_beta_diversity/`: per-sample richness/Shannon/Simpson/Pielou, summaries
by `type_ofsample` and `hf_name`, a Bray-Curtis distance matrix per region, and
a PERMANOVA-style permutation test per region × grouping variable. See the
`README.md` written into that output directory for the full method and headline
result.

### 6. Taxa barplots (local substitute for QIIME2's own)

```bash
python3 scripts/06_taxa_barplots.py --regions V1V2 V2V3 V3V4 V4V5 V7V9
```

Needs matplotlib. Stacked genus-composition PNG per requested region, built
from `compiled/genus_counts_by_region.tsv` -- useful on its own, and a
substitute when a region's `QIIME2_BARPLOT` process didn't produce one (as
happened for ITS/V5V7 in the first full run). Controls are excluded from the
plot; the top 8 genera get a fixed color across every region so the same
genus reads as the same color everywhere. See the `README.md` written into
`taxa_barplots/` for the full method.

### 7. Diversity boxplots

```bash
python3 scripts/07_diversity_boxplots.py
```

Needs matplotlib. Reads `alpha_diversity_per_sample.tsv` (from step 5) and
writes a 4-panel boxplot (richness/Shannon/Simpson/Pielou) per region ×
grouping variable into `alpha_beta_diversity/boxplots/`, each group
annotated with its sample count.

### 8. Per-sample export for downstream biostatistics

```bash
python3 scripts/08_export_for_biostatistician.py --regions V3V4 V4V5
```

Pure standard library. Reshapes the requested regions' own
`qiime2/rel_abundance_tables/rel-table-ASV_with-DADA2-tax.tsv` (ampliseq's
native per-ASV taxonomy + classifier confidence + relative abundance) into
one tidy row per sample × region × detected organism, written to
`organisms_by_sample/organisms_by_sample.tsv`. Unlike step 6, controls are
**included** (flagged via `is_control`), not filtered -- this is meant as a
complete handoff. See the `README.md` written into that output directory
for the full column reference.

### 9. Join the organism export with clinical metadata

```bash
python3 scripts/09_merge_organisms_with_metadata.py
```

Pure standard library. Left-joins step 8's output with `config/metadata.tsv`
on `sample_id`, adding every AFI/ARI demographic and pathogen-panel column
so downstream analysis can cross bacterial detections against clinical
results directly. Writes `organisms_by_sample/organisms_by_sample_with_metadata.tsv`
and warns if any sample has no metadata row. See that output directory's
`README.md` for column details and the ambiguous-specimen caveat.

### 10. Samples with metadata but no results

```bash
python3 scripts/10_samples_missing_results.py --regions V3V4 V4V5
```

Pure standard library. The reverse of `check_metadata.py`: finds samples
with real AFI/ARI clinical data (not a blank control row, not one of the
"no clinical metadata found" gaps) that have zero detection rows in step
8's export, and explains why per region -- excluded for low raw reads,
kept but DADA2 produced zero surviving ASVs (a merge/truncation failure
for that specific sample), or never sequenced in that region at all.

### 11. Rank-collapsed export (family / genus / species), including ITS

```bash
python3 scripts/11_export_by_rank.py --regions V3V4 V4V5 ITS --ranks family genus species
```

Pure standard library. Unlike step 8 (ASV-level, V3V4/V4V5 only, built from
ampliseq's merged QIIME2 table), this reads each region's own
`dada2/ASV_table.tsv` + `dada2/ASV_tax_species.<ref>.tsv` directly, so it
works the same way for ITS (UNITE-fungi reference, no merged QIIME2 table of
its own) as for the SILVA regions. ASVs sharing the same lineage up to the
named rank are merged (relative abundance summed, confidence
abundance-weighted, `n_asvs_collapsed` reported) -- the same semantics as
QIIME2's own `taxa collapse`. Writes `family_level.tsv`, `genus_level.tsv`,
`species_level.tsv` into `organisms_by_sample/`, controls included. See
that directory's `README.md` for the full column reference.

### 12. Rarefaction curves

```bash
python3 scripts/12_rarefaction_plots.py --regions V3V4 V4V5 ITS
```

Needs matplotlib. Expected ASV richness vs. sequencing depth, one line per
sample, using Hurlbert's exact rarefaction formula (deterministic, no
repeated subsampling needed). Lines are colored by specimen type, same fixed
mapping as step 7's boxplots. A dashed reference line marks the depth step
5 rarefied to for that region (lowest depth among samples with >=1,000
reads); no line is drawn if no sample clears that bar (as for ITS, which is
also why step 5 skips ITS entirely -- these curves make that visually
obvious). Writes `<REGION>_rarefaction.png` into
`alpha_beta_diversity/rarefaction_curves/`.

---

## Tuning truncation

This is the setting that decides whether a region yields anything, so the
defaults in `config/regions.tsv` are a starting point, not an answer.

Two constraints, both computed for this run's 2×276 bp reads:

**Budget** — a read only has `276 − 11 (max phasing) − primer_len` usable bases,
so `trunclenf` must stay under that or DADA2 discards the read.

**Overlap** — `trunclenf + trunclenr` must exceed the primer-to-primer insert by
DADA2's minimum overlap (12 bp) plus a margin.

| Region | Amplicon | Insert | trunclenf/r | Sum | Est. overlap | R1/R2 budget |
|--------|----------|--------|-------------|-----|--------------|--------------|
| V1V2 | 331 | 298 | 210 / 185 | 395 | ~97 | 247 / 250 |
| V2V3 | 418 | 384 | 225 / 200 | 425 | ~41 | 248 / 248 |
| V3V4 | 449 | 411 | 235 / 215 | 450 | ~39 | 248 / 244 |
| V4V5 | 412 | 373 | 225 / 200 | 425 | ~52 | 246 / 245 |
| V5V7 | 391 | 350 | 220 / 195 | 415 | ~65 | 244 / 245 |
| V7V9 | 393 | 356 | 225 / 195 | 420 | ~64 | 250 / 243 |
| ITS  | variable | variable | auto | — | — | — |

V3V4 has the least headroom and is the region most likely to need attention.
ITS1 is genuinely variable in length, so fixed truncation would silently delete
the long variants; it is set to `auto` (`--trunc_qmin`) instead.

**Two things to know about these reads.** They are already quality-trimmed to a
variable 35–276 bp (median 276), so any fixed `trunclen` discards every read
shorter than it — check `filtered` vs `DADA2_input` in each region's
`overall_summary.tsv`. And 6 of the 45 libraries have a median length of ~47 bp;
those are effectively empty and `--min-reads` will drop them.

If a region merges poorly, either lower `trunclenr` to buy overlap or switch that
row to `auto`, then rerun just that region — `-resume` keeps the rest.

```bash
scripts/03_run_ampliseq.sh V3V4        # after editing config/regions.tsv
```

---

## Validation

The demultiplexer was tested on a synthetic library built from QIAGEN's own
expanded primer file — 300 read pairs per region drawn from randomly chosen phase
variants (0–11), plus 200 random-sequence pairs as decoys:

- 2,300 / 2,300 pairs accounted for, **0 misassigned** across all 7 regions
- all 200 decoys correctly routed to `unknown`
- trimmed lengths matched `276 − phase − primer_len` exactly, confirming the
  phasing bases and primer are both removed; `unknown` reads stayed untrimmed

Regenerate with `--probe` on real data before trusting a new run.

Specificity note: the 7 forward primers are mutually distinguishable, as are the
7 reverse primers, and cutadapt only searches forward primers in R1 and reverse
primers in R2. So although V5V7's forward primer (785F) and V3V4's reverse primer
(805R) target the same conserved locus on opposite strands, they are never
compared against each other.

---

## Troubleshooting

**A large `unknown` fraction (>20%).** Run `--probe` to rule out swapped reads.
Failing that, relax `--error-rate` to 0.2 or lower `--min-overlap`, and check
`demux/logs/<sample>.cutadapt.log` for which primers are matching at all.
Some unassigned reads are normal — off-target amplification and primer dimers
land there by design.

**A region merges badly but others are fine.** Truncation, not demultiplexing.
See the tuning section.

**Host contamination.** The earlier pooled run's most abundant lineage was
`Rickettsiales;Mitochondria` (~150k reads) with another ~18k assigned to
`Eukaryota` — host mitochondrial and eukaryotic 16S, not bacteria. ampliseq's
`--filter_ssu bac` will exclude those if you want bacteria only.
`04_compile_results.py` already excludes Mitochondria/Chloroplast from every
compiled table by default (step 4), and prints per-region counts of what it
dropped. It's worst in Blood specimens — the diversity (step 5) and taxa
barplot (step 6) outputs quantify and visualize that split by specimen type.

**A whole region is empty.** Check that region's row in `demux_counts.tsv`. If it
got reads but produced no ASVs, it is truncation. If it never got reads, that
region likely was not in your panel — remove the row.

**Sample split across lanes.** Step 01 stops with an error rather than silently
using one lane. Concatenate per read first:
`cat X_S1_L00*_R1_001.fastq.gz > merged/X_S1_L001_R1_001.fastq.gz`.

---

## Layout

```
config/regions.tsv          the one file to edit: primers, truncation, taxonomy
config/ampliseq_base.yml    shared ampliseq params
config/primers_{fwd,rev}.fasta, readthrough.tsv
                            generated from regions.tsv on every run of step 01
config/metadata.tsv         generated by build_metadata.py, committed -- QIIME2
                            metadata for step 3's --extra "--metadata ..." and
                            for steps 5, 9, 10
metadata/                   source clinical data: AFI/ARI .dta files (Stata),
                            not derived from anything in this repo
reference/                  QIAGEN's expanded primer table (provenance)
scripts/lib/regions.py      shared config parsing, ID cleaning, FASTQ pairing
scripts/lib/samplesheets.py shared: collect sample IDs from samplesheets/,
                            read excluded_samples.tsv, control-ID detection
scripts/lib/taxonomy.py     shared: host-contaminant filtering, reading a
                            region's dada2 ASV table + taxonomy uniformly
                            (SILVA or UNITE-fungi)
demux/, samplesheets/, results_by_region/, work/, logs/, compiled/
                            generated (steps 0-4)
alpha_beta_diversity/, taxa_barplots/, organisms_by_sample/
                            generated (steps 5-12, local-analyst-side only)
```

`config/primers_*.fasta` are rewritten from `config/regions.tsv` at the start of
every step 01 run, so they cannot drift out of sync. Edit `regions.tsv`, never the
generated FASTAs.
