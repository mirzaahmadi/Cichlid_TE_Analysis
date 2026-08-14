# Cichlid Transposable Element (TE) Radiation Pipeline

A small, local pipeline for comparing **transposable element activity in the
neighborhoods of single-copy orthologous genes** across cichlid species. The
current implementation covers four African cichlids; it is designed to be
extended to South American cichlids next (see
[`FUTURE_DIRECTIONS.docx`](FUTURE_DIRECTIONS.docx)).

The core question: for genes that are shared 1:1 across species (single-copy
orthologs), **how much of the DNA immediately flanking each gene is made of
TEs, what kinds of TEs are they, and how much does this vary from species to
species?** Orthologs whose flanks are highly variable across species are the
interesting ones — they may reflect lineage-specific TE churn around conserved
genes.

---

## Species in the current analysis

| Species | Short code used in files |
|---|---|
| *Metriaclima zebra* | `Metriaclima_zebra_AF` |
| *Pundamilia nyererei* | `Pundamilia_nyererei_AF` |
| *Oreochromis niloticus* | `Oreochromis_niloticus_AF` |
| *Neolamprologus brichardi* | `Neolamprologus_brichardi_AF` |

`_AF` = African. (`Oreochromis` is the best-annotated of the four and is used as
the preferred source of gene symbols in the aggregation step.)

---

## Repository layout

```
cichlid-TE-pipeline/
├── README.md                     ← you are here
├── FUTURE_DIRECTIONS.docx        ← next steps / ideas from project meetings (Word doc)
├── scripts/
│   ├── extract_orthogroup_gene_info.py     STAGE 1 (Python)
│   ├── V3_te_flanking_analysis.R           STAGE 2 (R, run once per species)
│   ├── V3_cross_species_aggregation.R      STAGE 3 (R, run once)
│   └── top10_ortholog_neighborhoods.R      STAGE 4 (R, run once)
└── docs/
    ├── FULL_PIPELINE_WORKFLOW.docx   Illustrated, plain-language walkthrough
    │                                 of every step (39 screenshots). Read this
    │                                 first if you don't read code.
    └── PIPELINE_SCRIPT_MAP.html      One-page visual "script → output" map.
                                      Open in any browser.
```

> **Note on the data.** The raw genomes, RepeatMasker TE annotations (`.out`),
> GFF3 annotations, and OrthoFinder outputs are **not** included in this repo —
> the annotation/TE files alone are several hundred MB each and exceed GitHub's
> file-size limits. See [Data & where it comes from](#data--where-it-comes-from)
> below for what each script expects and how to obtain it.

---

## The pipeline at a glance

Run the four scripts **in order**. Each stage's output feeds the next. The
`docs/PIPELINE_SCRIPT_MAP.html` file is a visual version of this same table.

| # | Script | Language | Runs | Reads | Writes (key handoff in **bold**) |
|---|--------|----------|------|-------|----------------------------------|
| 1 | `extract_orthogroup_gene_info.py` | Python 3 | once | Orthogroups.tsv, single-copy list, 4× GFF3 | **`OUTPUT_TABLE_orthogroup_gene_info.tsv`** |
| 2 | `V3_te_flanking_analysis.R` | R | **once per species (×4)** | the table from #1, that species' TE `.out`, its GFF3 | per-species `*_intersection_outputs/` incl. **`master_table.tsv`** + **`density_profile.tsv`** |
| 3 | `V3_cross_species_aggregation.R` | R | once | the 4 species folders from #2 | heatmaps, combined density, **`cross_species_sd_top10.tsv`** |
| 4 | `top10_ortholog_neighborhoods.R` | R | once | the top-10 table from #3, the table from #1, 4× TE `.out` | per-ortholog neighborhood figures |

### What each stage does

**Stage 1 — `extract_orthogroup_gene_info.py`**
Takes the list of single-copy orthogroups plus `Orthogroups.tsv` (both from
OrthoFinder) and the four species GFF3 annotation files. For each ortholog it
traces the protein → its `GeneID` (via the CDS rows) → the gene's genomic
coordinates (via the gene rows), scanning each large GFF3 exactly once. Produces
one master coordinate table (`OUTPUT_TABLE_orthogroup_gene_info.tsv`) with one
row per (orthogroup, species, protein): chromosome, start, end, strand,
gene name.

**Stage 2 — `V3_te_flanking_analysis.R`** (run 4×, once per species)
For one species: parses the RepeatMasker `.out` TE annotation into a
`GenomicRanges` object, drops sub-100 bp fragments (noise), then for every
single-copy ortholog extends the gene by ±10 kb and intersects those flanking
windows with the TEs. It removes "orthologs" whose own gene body is ≥90 % TE
(transposon-derived), computes per-ortholog TE statistics (% of flank that is
TE, TE density, per-class breakdown), and writes a `master_table.tsv` plus
several diagnostic tables and plots (density along the flank by tier, %TE
distribution, TE-class composition, flank completeness, top-20 orthologs by TE,
genome-wide vs flank comparison). It also samples two random baselines (random
20 kb genome chunks and random gene flanks) for comparison.

**Stage 3 — `V3_cross_species_aggregation.R`** (run once, after all 4 species)
Reads the four `master_table.tsv` + `density_profile.tsv` files and compares
species. Produces a combined density plot, two orthogroup × species heatmaps
(absolute and row z-score, rows hierarchically clustered), a cross-species
standard-deviation table, and — the key handoff — `cross_species_sd_top10.tsv`,
the 10 orthogroups whose flank TE coverage varies **most** across the four
species.

**Stage 4 — `top10_ortholog_neighborhoods.R`** (run once, after Stage 3)
For each of those 10 most-variable orthologs, draws the ±10 kb genomic
neighborhood in all four species stacked together: the gene as an oriented box
(re-oriented so every panel reads 5'→3') and every surrounding TE as a rectangle
coloured by TE class. This is the figure set discussed in the later meetings.

---

## Requirements

- **Python 3** (standard library only — no packages to install for Stage 1).
- **R** with:
  - Bioconductor: `GenomicRanges`, `IRanges`
    ```r
    install.packages("BiocManager")
    BiocManager::install(c("GenomicRanges", "IRanges"))
    ```
  - CRAN: `ggplot2`, `dplyr`, `data.table`
    ```r
    install.packages(c("ggplot2", "dplyr", "data.table"))
    ```

Everything runs comfortably on a normal laptop — no cluster or high-memory
machine is needed.

---

## How to run

Each R script has a **`USER INPUTS`** block at the very top (file paths and the
species name). Edit that block, then run the whole script. The scripts use
**relative paths** (`data/...`, `insights/...`), so run each one from a working
directory that contains the `data/` and `insights/` folders it expects.

### Stage 1
```bash
python3 scripts/extract_orthogroup_gene_info.py Orthogroups_SingleCopyOrthologues.txt \
    --orthogroups-tsv Orthogroups.tsv \
    --gff-dir gff3_files \
    --output OUTPUT_TABLE_orthogroup_gene_info.tsv
```

### Stage 2 (repeat for all four species)
Open `scripts/V3_te_flanking_analysis.R` and set the `USER INPUTS` block for the
species you're running:
```r
te_file         <- "data/TE_<species>_genome.fa.out"
orthogroup_file <- "data/OUTPUT_TABLE_orthogroup_gene_info.tsv"
species_name    <- "<species>"                 # e.g. "Metriaclima_zebra_AF"
output_dir      <- "insights/<species>_intersection_outputs"
gff3_file       <- "data/GENE_<species>.annotation.gff3"   # optional but recommended
```
Then run it (e.g. `Rscript scripts/V3_te_flanking_analysis.R`). Do this **four
times**, once per species.

### Stage 3
Confirm the four `output_dir` paths in the `species_dirs` block match what
Stage 2 wrote, then run `scripts/V3_cross_species_aggregation.R`.

### Stage 4
Place `cross_species_sd_top10.tsv` (from Stage 3) and
`OUTPUT_TABLE_orthogroup_gene_info.tsv` (from Stage 1) plus the four TE `.out`
files where the `USER INPUTS` block expects them, then run
`scripts/top10_ortholog_neighborhoods.R`.

### ⚠️ The pipeline is not (yet) fully automated end-to-end
This is a research pipeline, not a one-click program. There is some **manual
file shuffling between stages** — e.g. copying `OUTPUT_TABLE_orthogroup_gene_info.tsv`
into each stage's `data/` folder, and moving `cross_species_sd_top10.tsv` from
Stage 3's output into Stage 4's `data/` folder. Each script's `USER INPUTS`
block tells you exactly which files it expects and where. Making the whole thing
flow automatically is listed as a next step in
[`FUTURE_DIRECTIONS.docx`](FUTURE_DIRECTIONS.docx).

---

## Data & where it comes from

None of the large data files live in this repo. Here is what each stage needs
and its origin:

- **OrthoFinder outputs** (`Orthogroups.tsv`,
  `Orthogroups_SingleCopyOrthologues.txt`) — produced by running
  [OrthoFinder](https://github.com/davidemms/OrthoFinder) on the four species'
  `protein.fa` files. In this project they were provided ready-made by a
  collaborator (Jordana).
- **GFF3 gene annotations** (`<species>.annotation.gff3`) — the standard genome
  annotation for each species.
- **RepeatMasker TE annotations** (`TE_<species>_genome.fa.out`) — RepeatMasker
  `.out` files, one per genome. Note these have **two** header lines plus a blank
  line (the parser skips 3 lines accordingly).

To reproduce the run, obtain these files, place them in a `data/` (and
`gff3_files/`) folder next to the scripts as shown above, and follow the run
order. If you are picking this project up, ask the project owner for the exact
data bundle used.

---

## Where to read more

- **`docs/FULL_PIPELINE_WORKFLOW.docx`** — the illustrated, plain-language
  walkthrough. Every step is described in natural language with screenshots, so
  you can follow the logic without reading a line of code. **Start here.**
- **`docs/PIPELINE_SCRIPT_MAP.html`** — a one-page visual of how scripts feed
  each other and exactly what files each writes. Open in a browser.
- The scripts themselves are heavily commented line-by-line.

---

## References & tools

Method/tool references (please verify links against the current published
versions):

- **OrthoFinder** — Emms DM & Kelly S. (2019). *OrthoFinder: phylogenetic
  orthology inference for comparative genomics.* Genome Biology.
- **RepeatMasker** — Smit AFA, Hubley R & Green P. `RepeatMasker` —
  https://www.repeatmasker.org
- **GenomicRanges / IRanges** — Lawrence M *et al.* (2013). *Software for
  computing and annotating genomic ranges.* PLoS Computational Biology
  (Bioconductor).
- **African cichlid genomes / adaptive radiation** — Brawand D *et al.* (2014).
  *The genomic substrate for adaptive radiation in African cichlids.* Nature.

*(A fuller, project-specific reading list can be added here — this section is a
starting point.)*

---

## Handoff note

This repository is the transparent handoff of the cichlid TE pipeline for the
next person to continue the work (including porting it to the South American
cichlids). The code is complete and commented; the documentation in `docs/`
explains every step in plain language. See
[`FUTURE_DIRECTIONS.docx`](FUTURE_DIRECTIONS.docx) for the agreed next steps and
open ideas from the project meetings.
