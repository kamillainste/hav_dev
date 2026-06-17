# HAV Batch Analysis Pipeline

Hepatitis A Virus (HAV) genomic analysis pipeline for genotyping, lineage
assignment, and phylogenetic placement of clinical sequences. It produces an
HTML report summarising all results per batch.

The pipeline accepts **two kinds of input**:

| Mode | Input | Assembly |
|------|-------|----------|
| `sanger` | Short-fragment **Sanger** sequences (VP1/2A or VP1‑2B junction), one consensus per sample | None — consensus used directly |
| `wgs` | **Nanopore** whole-genome tiling-PCR reads (raw FASTQ per barcode) | NanoPlot QC → Chopper filtering → ViroConstrictor |

> **Key design point — everything is short-fragment.**
> Regardless of input length, **all downstream classification uses the
> short-fragment VP1‑2B resources**: the `hav-vp1-2b-lineages` NextClade dataset
> and the short-fragment local BLAST database. A whole-genome consensus is
> simply aligned/placed against the VP1‑2B region. There is no whole-genome
> classification step — the whole-genome NextClade dataset present in the repo
> is **not** used for clade/lineage assignment.

---

## Setup

The pipeline uses several conda environments. They are **all created by a
single setup script**, `scripts/setup_havdev_env.sh`, from environment files
checked into the repo.

### Conda environments used

| Environment | Defined in | Used by | Purpose |
|-------------|-----------|---------|---------|
| `HAVDEV` (project-local, `./.conda/HAVDEV`) | `environment.yml` | All modes | nextclade, blastn, mafft, iqtree, seqkit, augur, jq, Python |
| `R_shared` (`$HOME/.conda/R_shared`) | `r_environment.yml` | All modes | R + report packages; pipeline calls `$HOME/.conda/R_shared/bin/Rscript` |
| `nanoplot` (named) | `envs/nanoplot.yml` | WGS only | NanoPlot read QC |
| `chopper` (named) | `envs/chopper.yml` | WGS only | Read quality/length filtering |
| `viroconstrictor` (named) | `envs/viroconstrictor.yml` | WGS only | Amplicon assembly → consensus genomes |

`R_shared` is kept deliberately separate from the `HAVDEV` tools environment and
lives under your home directory (`$HOME/.conda/R_shared`) so it can be shared
across projects. The WGS environments are invoked via `conda run -n <env>` from
`run_pipeline.sh`, so they must exist as **named** conda environments.

### Prerequisites

- [WSL 2](https://learn.microsoft.com/en-us/windows/wsl/install) (Windows) or native Linux
- [Miniconda / Miniforge](https://docs.conda.io/en/latest/miniconda.html) on `PATH`

### Install

```bash
cd /path/to/hav_dev

# Create/update ALL environments (HAVDEV + R_shared + WGS envs)
bash scripts/setup_havdev_env.sh

# Or skip the WGS environments if you only run Sanger batches:
bash scripts/setup_havdev_env.sh --core

# Activate HAVDEV for running the pipeline
conda activate /path/to/hav_dev/.conda/HAVDEV
```

`setup_havdev_env.sh` is idempotent — re-run it to update existing environments
after changes to any of the environment files. It always sets up `HAVDEV` and
`R_shared`; the three WGS environments are added unless you pass `--core`.

---

## Quick start — `hav_analyse.sh`

`scripts/hav_analyse.sh` is the **single entry point** for a batch. It runs the
mode-specific input preparation and then the shared analysis core (BLAST +
NextClade → per-sequence trees → HTML report). Run it from the project root
with `HAVDEV` active.

```bash
# Sanger batch (no primer trimming)
bash scripts/hav_analyse.sh --mode sanger batches/Batch-1 \
  --samplesheet batches/Batch-1/samplesheet.tsv

# Sanger batch with reference-coordinate primer trimming
bash scripts/hav_analyse.sh --mode sanger batches/Batch-1 \
  --samplesheet batches/Batch-1/samplesheet.tsv \
  --primer-names "HAV_6.1_codehop,HAV_10_codehop,HAV_8.2_codehop,HAV_11_codehop"

# WGS batch — full pipeline from raw Nanopore reads
bash scripts/hav_analyse.sh --mode wgs batches/HAV-2026_01 \
  --samplesheet batches/HAV-2026_01/samplesheet.tsv

# WGS batch — assembly already done, re-run analyses/report only
bash scripts/hav_analyse.sh --mode wgs batches/HAV-2026_01 \
  --samplesheet batches/HAV-2026_01/samplesheet.tsv --skip-assembly
```

Common options: `--dataset-date <date>` (default `2026-04-10`), `--threads <n>`,
`--n-neighbors <n>` (tree neighbours per query, default 30), `--skip-trees`,
`--skip-report`. Run `bash scripts/hav_analyse.sh --help` for the full list.

The batch directory (e.g. `batches/Batch-1`) is the output root. All results are
written under `<batch_dir>/output/`, with a full run log at
`<batch_dir>/output/hav_analyse.log`.

---

## Input: the samplesheet

Both modes are driven by a **TSV samplesheet** with a header row. Sanger needs
only two columns; WGS uses the full ViroConstrictor schema (the same sheet works
for both — the extra columns are ignored in Sanger mode).

**Required columns (full schema):**

`Sample`, `Virus`, `MatchRef`, `Segmented`, `Primers`, `Reference`,
`Features`, `MinCov`, `Mismatch`, `InputDir`

- `Sample` — sequence ID used everywhere downstream (e.g. `barcode81`). Must be unique.
- `InputDir` — per-sample input directory:
  - **Sanger:** a directory containing one `.fa`, `.fasta`, or `.TXT` file
    (raw Sanger contig text). The first matching file is used; the `Sample`
    column overrides any embedded FASTA header.
  - **WGS:** a directory containing the raw `*.fastq.gz` reads for that barcode.
- `Virus`, `MatchRef`, `Segmented`, `Primers`, `Reference`, `Features`,
  `MinCov`, `Mismatch` — ViroConstrictor parameters (WGS assembly). `Primers`
  and `Reference` point at files under `data/Primers/` and `data/References/`.

See `batches/HAV-2026_01/samplesheet.tsv` for a working WGS example.

> **Sanger note:** `prepare_input_fasta.R` only needs `Sample` + `InputDir`, so a
> minimal two-column sheet is valid for `--mode sanger`. The WGS assembly path
> (`parse_samplesheet.R`) enforces the full schema.

---

## What each stage does (and which script runs it)

`hav_analyse.sh` orchestrates the stages below. You can also run any stage
directly for debugging.

### WGS-only preparation

1. **`run_pipeline.sh`** `<samplesheet> <batch_dir> [threads]`
   - `parse_samplesheet.R` — validates the samplesheet (strict schema), writes a
     per-sample loop file and a validated VC source sheet.
   - **NanoPlot** — read QC on raw reads (and again on filtered reads), per sample.
   - **Chopper** — quality/length filtering (`Q10`, min length 200 bp).
   - `write_viroconstrictor_samplesheet.R` — drops empty/0-read samples and
     renames columns to the ViroConstrictor format.
   - **ViroConstrictor** — Nanopore amplicon assembly (`--amplicon-type
     end-to-end`) → consensus genomes.
2. **Collect consensus** (in `hav_analyse.sh`) — concatenates
   `viroconstrictor/.../all_consensus.fasta` into `<batch>/<batch>.fa`, trimming
   each header to the sample name so all tools agree on the ID.

### Sanger-only preparation

1. **`prepare_input_fasta.R`** `<samplesheet> <batch>.fa` — collects one consensus
   per sample (from `.fa`/`.fasta`/`.TXT`) into the batch multi-FASTA.

### Shared analysis core (both modes)

3. **`run_all_analyses.sh`** `<batch_dir> [dataset_date] [--sanger ...]`
   - *(Sanger, optional)* **`trim_primers_by_reference.py`** — trims primer-covered
     ends using reference coordinates rather than sequence-end matching. Enabled
     via `--sanger --primer-names ...`.
   - **`blast_batch.sh`** — `blastn` of the batch FASTA against the local
     short-fragment HAV BLAST database. Writes `output/blast_results.tsv`.
   - **NextClade** — runs against `data/nextclade_datasets/hav-vp1-2b-lineages`,
     assigning genotype `clade` (IA/IB/IIA/IIIA) + `lineage_phylo`. Writes
     `output/lineages/{nextclade.tsv,nextclade.ndjson,nextclade.aligned.fasta}`.
4. **`build_per_seq_trees.sh`** `<batch_dir> [dataset_date] [n_neighbors]` — for
   each query, gathers nearest neighbours from BLAST hits + NextClade
   nearest-node sequences, aligns with MAFFT, builds an IQ-TREE ML tree.
   Writes `output/trees/<seqName>/`.
5. **`batch_report.Rmd`** — renders the HTML report
   (`output/batch_report.html`) plus `output/summary.tsv` and
   `output/comparison.tsv`.

### Report contents

- Per-sequence summary table (clade, lineage, QC status, BLAST agreement)
- Clade assignment comparison matrix
- Mutational profiles (substitutions, private mutations)
- BLAST hit details with similarity plots
- Per-sequence phylogenetic trees with SNP-distance heatmaps
- Multiple sequence alignments (VP1‑2B region)
- *(WGS only)* Amplicon coverage + primer-scheme panels from ViroConstrictor

---

## Outputs

```
<batch_dir>/
├── <batch>.fa                  # Batch multi-FASTA (Sanger consensus or WGS consensus)
├── <batch>_blast_results.tsv   # BLAST output (also copied into output/)
├── samplesheet.tsv             # Batch samplesheet
├── pipeline.log                # WGS assembly log (run_pipeline.sh)
├── Nanoplot/ filtered_reads/ viroconstrictor/   # WGS intermediates
└── output/
    ├── hav_analyse.log         # Full run log
    ├── blast_results.tsv
    ├── lineages/               # NextClade results (hav-vp1-2b-lineages)
    ├── trees/<seqName>/        # Per-sequence ML trees + alignments
    ├── summary.tsv
    ├── comparison.tsv
    └── batch_report.html       # Final report
```

Batch folders live under top-level `batches/` and are **git-ignored** (raw
reads and outputs are not committed).

---

## Dataset Management Workflow (Curator)

The local dataset has a **single purpose**: it is the reference collection
backing the **short-fragment BLAST database**. There are no longer any
clade-separated datasets or local reference trees — lineage assignment is
handled entirely by the public-derived `hav-vp1-2b-lineages` NextClade dataset.

### Refresh the local BLAST dataset

```bash
# 1. Export from the sequence database into data/local_datasets/<YYYY-MM-DD>/:
#      export.csv   — metadata (Key, Genotype, OUTBREAK_VARIANT, …)
#      2PA.fa       — all sequences (FASTA)

# 2. Run the BLAST-only update (prepare_dataset.R → metadata_corrected.tsv → make_blast_db.sh)
bash scripts/update_datasets.sh <YYYY-MM-DD>
```

To merge in external reference sequences before building the database, use
`combine_sequences.sh <date>` (merges `input.fa` with the external SQLite DB at
`data/local_datasets/externe_sek/externe.db`, managed via `R-database.R`).

### Rebuild the lineage NextClade dataset (public-only)

```bash
bash scripts/build_vp1b_lineage_dataset.sh
# Writes: data/nextclade_datasets/hav-vp1-2b-lineages/
```

This fetches the community `hav-vp1-2b` dataset fresh from the NextClade server
and designates lineages **directly on the public tree** — it has no dependency
on `data/local_datasets/` and uses no private FHI sequences. It also regenerates
the SILO `lineage_definitions/lineages.yaml` hierarchy.

---

## Lineage Designation Workflow

When sequences cluster into a candidate new lineage:

```bash
# 1. Scan the reference tree for qualifying monophyletic groups
python3 scripts/designate_lineages.py <date>     # writes data/lineage_proposals.tsv

# 2. Review proposals (mutations, dates, countries, outbreaks).
#    Proposals are pre-added to the register; no action needed to ACCEPT.
python3 scripts/review_lineage_proposals.py

# 3. To REJECT/rename: edit data/lineage_designations.tsv
#    (set deprecated=TRUE for a row, or rename as appropriate)

# 4. Rebuild the lineage dataset to activate changes
bash scripts/build_vp1b_lineage_dataset.sh
```

---

## Project Structure

```
hav_dev/
├── batches/                        # Per-batch input + output (git-ignored)
│   ├── Batch-<N>/                  # Sanger batches
│   └── HAV-2026_01/                # WGS (Nanopore) batches
├── data/
│   ├── local_datasets/
│   │   ├── <YYYY-MM-DD>/           # BLAST reference snapshot (per build date)
│   │   │   ├── 2PA.fa / export.csv         # Raw export from sequence DB
│   │   │   ├── input.fa                    # Cleaned internal sequences
│   │   │   ├── metadata.tsv / metadata_corrected.tsv
│   │   │   └── blast_db/                    # BLAST database + input_dedup.fa
│   │   └── externe_sek/externe.db          # External reference seqs (SQLite)
│   ├── nextclade_datasets/
│   │   ├── hav-vp1-2b/             # Community VP1-2B genotyping dataset
│   │   ├── hav-vp1-2b-lineages/    # Community dataset + FHI lineage labels (USED)
│   │   └── hav-whole-genome/       # Whole-genome dataset (NOT used for classification)
│   ├── lineage_definitions/lineages.yaml   # SILO/LAPIS hierarchical lineages
│   ├── lineage_designations.tsv            # Canonical lineage register
│   ├── lineage_proposals.tsv               # Candidate lineages awaiting review
│   ├── Primers/                            # PCR_primers/ + Tiling_PCR/ (FASTA, GFF3)
│   └── References/                         # Per-genotype reference genomes
├── scripts/                        # See script reference below
├── envs/                           # WGS conda env files (nanoplot, chopper, viroconstrictor)
├── environment.yml                 # HAVDEV environment (core tools + Python)
├── r_environment.yml               # R_shared environment (R + report packages)
└── NEXTCLADE_GUIDE.md              # Detailed NextClade operational guide
```

### Script reference

| Script | Role |
|--------|------|
| `hav_analyse.sh` | **Master wrapper** — dispatches by `--mode`, runs the full pipeline |
| `run_pipeline.sh` | WGS: NanoPlot QC → Chopper → ViroConstrictor assembly |
| `parse_samplesheet.R` | Validate samplesheet; emit loop file + VC source sheet |
| `write_viroconstrictor_samplesheet.R` | Build the ViroConstrictor samplesheet from valid samples |
| `prepare_input_fasta.R` | Sanger: samplesheet → batch multi-FASTA |
| `run_all_analyses.sh` | Shared core: BLAST + NextClade (+ optional Sanger primer trim) |
| `trim_primers_by_reference.py` | Coordinate-based Sanger primer trimming |
| `blast_batch.sh` | `blastn` of batch vs local short-fragment BLAST DB |
| `build_per_seq_trees.sh` | Per-sequence MAFFT + IQ-TREE phylogenies |
| `batch_report.Rmd` | HTML report (+ summary.tsv, comparison.tsv) |
| `update_datasets.sh` | Refresh local BLAST dataset (BLAST-only flow) |
| `prepare_dataset.R` | Clean internal database export → input.fa + metadata |
| `combine_sequences.sh` | Merge internal + external sequences |
| `R-database.R` | Manage the external `externe.db` SQLite database |
| `make_blast_db.sh` | Build the BLAST nucleotide database |
| `build_vp1b_lineage_dataset.sh` | Build `hav-vp1-2b-lineages` from the public tree |
| `designate_lineages.py` | Identify candidate new lineages on the reference tree |
| `review_lineage_proposals.py` | Display lineage proposals for curator review |
| `extract_auspice_tree.py` | Extract tree + tip metadata from an Auspice JSON |
| `patch_vp1b_tree.py` | Patch lineage labels onto the community tree |
| `download_community_seqs.py` | Fetch community reference sequences |
| `generate_silo_hierarchy.py` | Generate `lineages.yaml` (SILO/LAPIS hierarchy) |
| `generate_primer_coordinates.py`, `make_primer_alignment.py` | Primer coordinate/alignment utilities |
| `make_microreact_metadata.py` | Build Microreact metadata |
| `setup_havdev_env.sh` | Create/update all conda environments (HAVDEV + R_shared + WGS) |

---

## Key Datasets

| Dataset | Location | Purpose |
|---------|----------|---------|
| `hav-vp1-2b` | `data/nextclade_datasets/hav-vp1-2b/` | Community genotyping dataset (IA/IB/IIA/IIIA), fetched from NextClade. Source for the lineages dataset. |
| `hav-vp1-2b-lineages` | `data/nextclade_datasets/hav-vp1-2b-lineages/` | **Used for all classification.** Community dataset extended with FHI lineage labels (IA.1, IB.1.1, …). Built by `build_vp1b_lineage_dataset.sh`. |
| `hav-whole-genome` | `data/nextclade_datasets/hav-whole-genome/` | Present but **not** used for clade/lineage assignment. |
| Local BLAST DB | `data/local_datasets/<date>/blast_db/` | Short-fragment HAV nucleotide DB (internal + external). Backs `blast_batch.sh` and neighbour selection. |

---

## Collaboration / Handover

All team members work locally via WSL/Linux with the shared `environment.yml`.
Sync changes through Git commits and pushes.

```bash
# One-time (or after pulling env changes): set up all environments
bash scripts/setup_havdev_env.sh

# Activate the analysis environment
conda activate /path/to/hav_dev/.conda/HAVDEV
```

Handover notes:
- Run `bash scripts/setup_havdev_env.sh` once on a new machine — it builds
  `HAVDEV`, `R_shared`, and the three WGS environments from the checked-in env
  files. Use `--core` to skip the WGS environments (Sanger-only machines).
- Start a batch with `hav_analyse.sh` — it is the only command you normally run.
- R runs from `$HOME/.conda/R_shared` (not from `HAVDEV`); the scripts reference
  it via `$HOME`, so it is portable across user accounts.
- `batches/` is git-ignored; commit dataset/script changes, not batch outputs.
- See `NEXTCLADE_GUIDE.md` for NextClade dataset details and
  `docs/HAV_WGS_classification_plan.md` for the roadmap.
