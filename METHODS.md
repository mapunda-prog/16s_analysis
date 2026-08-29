# Methods

## Amplicon library and primers

Samples were sequenced on the QIAseq 16S/ITS Screening Panel, which
co-amplifies seven variable-region/ITS targets per sample in one library
(V1V2, V2V3, V3V4, V4V5, V5V7, V7V9, ITS1) using phased primers (0–11
staggering bases prepended per primer, to improve base-calling diversity).
Reads were paired-end, 2×276 bp. Primer sequences (phase-0 core, degenerate
bases collapsed to IUPAC) were recovered from QIAGEN's own demultiplexer
reference file (`16S_primers_ScreeningPanel_Demultiplex_Final.txt.expanded`,
shipped inside the official `qiaseq/qiaseq-16s` container) rather than
assumed from public primer databases.

## Region-wise demultiplexing

Because all seven amplicons share one library, reads were split into
per-region FASTQ files before any downstream processing — pooling
amplicons of different lengths into one DADA2 run is invalid (a single
truncation length cannot serve incompatible amplicon lengths, and the
error model would be fit across non-comparable data).

Demultiplexing used **cutadapt** (≥3.4 required; developed/tested against
5.2) with paired-adapter matching (`--pair-adapters`): a read pair was
assigned to a region only when R1 carried that region's forward primer
**and** R2 carried the matching reverse primer, matched with an unanchored
5′ search so any phasing prefix was removed together with the primer in one
trim. Parameters:

| Parameter | Value |
|---|---|
| Error rate (`-e`) | 0.15 |
| Minimum overlap (`-O`) | 15 nt (shortest primer across the full panel; applied uniformly to every region) |
| Minimum length post-trim | 100 bp |

Read pairs matching no region's primer pair went to an `unknown` bin
(inspected, not discarded) rather than being silently dropped. For ITS,
where amplicon length can be shorter than the read, a second cutadapt pass
removed 3′ read-through of the opposite primer.

## Per-region sample sheets and inclusion threshold

Samples with fewer than **1,000 read pairs** in a given region (after
demultiplexing) were excluded from that region's ASV inference — DADA2
fits a per-sample error model, and a near-empty sample can degrade or fail
the model for an entire region's run. Exclusion is per region: a sample
absent from one region for this reason can still be present and analyzed
in another.

## ASV inference (nf-core/ampliseq / DADA2)

Each region was processed as an independent nf-core/ampliseq run (own
working directory and output, so regions never share state), primers
already removed (`--skip_cutadapt`). Truncation lengths were set per region
from empirical read-quality and per-region amplicon length, subject to two
constraints: (i) budget — `trunclenf` must leave enough usable bases after
subtracting phasing (≤11 bp) and primer length; (ii) overlap — `trunclenf +
trunclenr` must exceed the primer-to-primer insert size by DADA2's minimum
merge overlap (12 bp) plus a margin.

| Region | Amplicon (bp, incl. primers) | trunclenF | trunclenR |
|---|---|---|---|
| V1V2 | 331 | 210 | 185 |
| V2V3 | 418 | 225 | 200 |
| V3V4 | 449 | 235 | 215 |
| V4V5 | 412 | 225 | 200 |
| V5V7 | 391 | 220 | 195 |
| V7V9 | 393 | 225 | 195 |
| ITS | variable (~230–450) | auto (quality-based) | auto (quality-based) |

DADA2 steps and parameters (as recorded by the pipeline; defaults except
`truncLen`, which is per-region above):

| Step | Key parameters |
|---|---|
| `filterAndTrim` | `truncLen` per region above; `maxN = 0`; `truncQ = 2`; `maxEE = c(2, 2)`; `minLen = 50`; `rm.phix = TRUE` |
| `learnErrors` | `nbases = 1e8`, sampled per region across all its retained samples jointly |
| `dada` (ASV inference) | default parameters, `pool = FALSE` (per-sample, error model shared per region) |
| `mergePairs` | `minOverlap = 12`, `maxMismatch = 0` |
| `removeBimeraDenovo` | `method = "consensus"` (default UCHIME-based chimera removal) |

Merged, chimera-filtered ASVs shorter than **100 bp** (`min_len_asv`) were
then excluded (ampliseq's `FILTER_LEN_ASV` step).

## Taxonomic classification

Bacterial/archaeal regions (V1V2, V2V3, V3V4, V4V5, V5V7, V7V9) were
classified against **SILVA 138.2** (prokaryotic SSU). ITS was classified
against **UNITE general FASTA release for Fungi, version 10.0** (release
04.04.2024). Both use DADA2's RDP-style naive Bayesian classifier
(`assignTaxonomy`, `minBoot = 50` — i.e. a rank is only reported if it
clears 50% bootstrap support; lower-confidence ranks are left blank rather
than guessed), followed by exact-match species-level assignment
(`addSpecies`, `tryRC = FALSE`). The reported classifier **confidence**
value used throughout this analysis is this bootstrap-support proportion.

## Host/reagent contaminant exclusion

ASVs classified as **Mitochondria** or **Chloroplast** at any taxonomic
rank were excluded from all downstream compiled tables, diversity
calculations, and biostatistics exports. This mirrors nf-core/ampliseq's
own default QIIME2-stage filter (`--exclude_taxa mitochondria,chloroplast`)
and matters most for blood-derived specimens, where host mitochondrial DNA
can dominate a low-biomass sample's reads entirely.

Two genera showing a strongly blood-specimen-concentrated distribution
pattern — *Rhizobiaceae* (also present in the NTC negative control, and the
single most abundant genus overall in this dataset) and *Solitalea*
(concentrated in blood specimens but absent from all controls) — were
flagged as likely or possible background signal respectively, and treated
with caution rather than excluded outright (see Results/Discussion).

## Diversity analyses

Alpha and beta diversity were computed per region on ASV counts (post
contaminant exclusion), joined to sample metadata (specimen type: Blood /
Swab / Blood+Swab [ambiguous, see below]; collection site: 4 health
facilities).

- **Inclusion**: samples with fewer than 1,000 reads (post-exclusion) in a
  region were dropped from that region's diversity analysis. A region was
  skipped entirely if fewer than 2 samples cleared this bar.
- **Rarefaction**: retained samples were rarefied once to the lowest
  surviving depth in that region (single random draw, fixed seed 100 for
  reproducibility).
- **Alpha diversity metrics**: observed ASV richness, Shannon index
  (natural log base), Simpson index, Pielou's evenness — computed on the
  rarefied counts.
- **Beta diversity**: Bray-Curtis dissimilarity on rarefied counts, per
  region.
- **Group testing**: a PERMANOVA-style permutation test (Anderson 2001
  pseudo-F statistic; 999 permutations) for each grouping variable
  (specimen type, collection site), restricted to groups with ≥2 samples.
- **Rarefaction curves** (QC visualization, not used for the diversity
  metrics above): expected ASV richness as a function of sequencing depth,
  computed exactly via Hurlbert's (1971) formula (no repeated random
  subsampling), for the full, unrarefied per-sample ASV counts.

## Taxonomic aggregation for biostatistics export

For downstream biostatistical analysis, per-ASV results were also
collapsed to family, genus, and species level (independently, as three
separate tables) by grouping ASVs sharing an identical lineage up to the
named rank — the same semantics as QIIME2's `taxa collapse`. Relative
abundance was summed across merged ASVs; the reported confidence at a
collapsed rank is the **abundance-weighted mean** of the constituent ASVs'
own classifier confidence (a derived summary — there is no native
"confidence of a genus" in the underlying data).

## Clinical metadata integration

Two clinical panels were joined to sequenced samples by enrolment ID: AFI
(Acute Febrile Illness — blood specimens; pathogens: brucellosis,
coxiellosis, dengue, chikungunya, Zika, West Nile virus) and ARI (Acute
Respiratory Infection — swab specimens; respiratory viral/bacterial panel).
Where a sequenced sample ID carried an explicit specimen-type suffix
(`BL`/`SW` or `B`/`S`), the matching panel alone was used. Where an
enrolment had both panels available but the sequenced ID carried no such
suffix (8 samples in the current dataset), both panels were merged onto
that sample's row rather than assigning one arbitrarily, and the affected
samples are flagged (`panel_source = AFI+ARI`) for manual specimen-type
confirmation rather than treated as resolved.

## Software and tool versions

| Tool | Version | Role |
|---|---|---|
| nf-core/ampliseq | 2.14.0 | Per-region ASV pipeline orchestration |
| Nextflow | 24.10.5 (≥24.04.2 required) | Workflow execution |
| Singularity | — (container profile) | Containerized execution |
| Java (OpenJDK) | 21 (17–22 required; Nextflow's launcher rejects >22) | Nextflow runtime |
| cutadapt | ≥3.4 required (developed/tested against 5.2) | Region-wise primer demultiplexing (custom step, outside ampliseq) |
| DADA2 (R package) | 1.30.0 | Denoising, ASV inference, chimera removal, taxonomy |
| R | 4.3.2 (core DADA2 steps) | DADA2 runtime |
| QIIME2 | 2024.10.1 | Feature table / taxonomy export, relative-abundance tables |
| barrnap | 0.9 | rRNA domain summary |
| FastQC | 0.12.1 | Read quality reports |
| phyloseq (R) | 1.46.0 | Phyloseq object construction |
| TreeSummarizedExperiment (R) | 2.10.0 | SummarizedExperiment export |
| Python | ≥3.6 (pipeline scripts); 3.9.1 (ampliseq's own Python steps) | Demultiplexing/compilation scripts, downstream analysis |
| pandas | — (local analysis only; ampliseq's own Python steps used 1.1.5) | Clinical metadata parsing |
| numpy | — (local analysis only) | Diversity/rarefaction computation |
| matplotlib | — (local analysis only) | Figure generation |
| SILVA | 138.2 (prokaryotic SSU) | 16S taxonomic reference |
| UNITE | General FASTA release for Fungi, v10.0 (2024-04-04) | ITS taxonomic reference |

Exact per-process software versions (including minor auxiliary steps not
listed above) are recorded in each region's
`pipeline_info/software_versions.yml`, written automatically by
nf-core/ampliseq for full reproducibility.

## Summary of thresholds and cutoffs

| Parameter | Value | Applies to |
|---|---|---|
| Cutadapt error rate | 0.15 | Region demultiplexing |
| Cutadapt minimum overlap | 15 nt (shortest primer, applied globally) | Region demultiplexing |
| Minimum trimmed read length | 100 bp | Region demultiplexing |
| Minimum read pairs per sample per region | 1,000 | Region sample-sheet inclusion |
| filterAndTrim `maxEE` | 2, 2 (fwd, rev) | DADA2 quality filtering |
| filterAndTrim `truncQ` | 2 | DADA2 quality filtering |
| filterAndTrim `minLen` | 50 bp | DADA2 quality filtering |
| mergePairs `minOverlap` | 12 bp | DADA2 read merging |
| mergePairs `maxMismatch` | 0 | DADA2 read merging |
| Minimum ASV length | 100 bp | Post-merge ASV filtering |
| Taxonomy bootstrap confidence | 50% (`minBoot`) | Taxonomic classification |
| Contaminant taxa excluded | Mitochondria, Chloroplast (any rank) | All compiled/exported tables |
| Minimum reads for diversity inclusion | 1,000 (post-contaminant-exclusion) | Alpha/beta diversity |
| Minimum samples for a region's diversity analysis | 2 | Alpha/beta diversity |
| Rarefaction | single draw to lowest surviving depth, seed 100 | Alpha/beta diversity |
| PERMANOVA permutations | 999 | Beta diversity group testing |
| Minimum group size for PERMANOVA | 2 samples | Beta diversity group testing |
