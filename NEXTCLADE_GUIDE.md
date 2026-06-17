# NextClade Quick Start Guide

## What is NextClade?

NextClade is a tool for viral genome alignment, mutation calling, clade assignment, quality checks and phylogenetic placement.

## Environment Setup (WSL + Conda)

This project uses the `HAVDEV` conda environment running inside WSL on Windows.

### First-time setup

If the environment does not yet exist, create it from the provided YAML file:

```bash
# In WSL terminal
conda env create -f environment.yml
```

If the environment already exists and you want to update it:

```bash
conda env update -f environment.yml --prune
```

### Activate the environment

```bash
conda activate HAVDEV
```

Verify NextClade is available:

```bash
nextclade --version
```

### Opening a WSL terminal in VS Code

1. Open the Command Palette (`Ctrl+Shift+P`) and run **Terminal: Select Default Profile**
2. Choose **WSL Bash** (or open a new terminal with the WSL profile)
3. Activate the project-local environment:
  `conda activate /path/to/hav_dev/.conda/HAVDEV`

### Accessing network drives (e.g. N:) from WSL

WSL only auto-mounts local drives. Mapped network drives must be mounted manually.
Make sure the drive is connected in Windows first, then run:

```bash
sudo mkdir -p /mnt/n
sudo mount -t drvfs N: /mnt/n
```

You can then access paths like:

```bash
ls "/mnt/n/Virologi/Hepatitt/Hepatitt A/HAV genotyping/2026/HBSPCR1-0407-1"
```

To make this persistent across WSL sessions, add it to `/etc/fstab`:

```bash
echo 'N: /mnt/n drvfs defaults 0 0' | sudo tee -a /etc/fstab
```

---

## HAV Genotyping Workflow (VP1-2B Junction)

This project genotypes Hepatitis A virus from Sanger sequencing contigs using the
`community/masphl-bioinformatics/hav/vp1-2b-junction` dataset.

### Dataset

The VP1-2B junction fragment is the standard target for HAV genotyping.
The NextClade dataset covers clades IA, IB, IC, IIA, IIB, IIIA, IIIB.

### Step 1 â€” Download the HAV dataset (once)

```bash
conda activate HAVDEV

mkdir -p data/hav-vp1-2b
nextclade dataset get \
  --name community/masphl-bioinformatics/hav/vp1-2b-junction \
  --output-dir data/hav-vp1-2b
```

### Step 2 â€” Prepare input FASTA

Use `scripts/prepare_input_fasta.R batch_dir out_dir` to produce a single multi-FASTA from a batch directory.
The script handles `.fa`, `.fasta`, and `.TXT` (raw Sanger contigs) files and writes the
output to `data/<batch_name>/<batch_name>.fa`.

```bash
conda activate HAVDEV
Rscript scripts/prepare_input_fasta.R data/Batch-1 data/Batch-1/output
```

The resulting file (e.g. `data/Batch-1/Batch-1.fa`) is the input for Step 3.

### Step 3 â€” Run NextClade with built-in database

**Purpose:** assign the HAV genotype (clade) to each sample using the community VP1-2B dataset.

```bash
mkdir -p data/Batch-1/output

nextclade run \
  --input-dataset data/hav-vp1-2b \
  --alignment-preset high-diversity \
  --output-tsv data/Batch-1/output/nextclade.tsv \
  --output-fasta data/Batch-1/output/nextclade.aligned.fasta \
  --output-json data/Batch-1/output/nextclade.json \
  data/Batch-1/Batch-1.fa
```

The key output column is `clade`, which gives the HAV genotype (IA, IB, IIA, IIB, IIIA, IIIB).

### Run NextClade with local database

**Purpose:** assign the outbreak variant/lineage and obtain an accurate mutation profile
relative to a locally curated reference set.

#### Prepare the local dataset files (only needed if updated sequences)
##### 4.1 â€” Prepare the local dataset files

Run the preparation script once per dataset update to produce `input.fa` and `metadata.tsv`
from the exported sequence database:

```bash
cd /path/to/hav_dev
Rscript scripts/prepare_dataset.R
```

This writes two files into `data/local_datasets/<date>/`:

| File | Description |
|------|-------------|
| `input.fa` | Reference sequences with clean key-only FASTA headers |
| `metadata.tsv` | Metadata with columns `id`, `genotype`, `lineage`, `date`, `country` |

##### 4.2 â€” Build the local NextClade datasets (one per genotype)

A separate dataset is built for each HAV genotype (IA, IB, IIA, IIB, IIIA, IIIB) so that
lineage/variant assignment and mutation calling are done against genotype-matched references.

Each dataset lives in its own directory and contains six files:

```
data/local_datasets/<date>/<genotype>/
â”œâ”€â”€ reference.fasta          # one representative sequence for this genotype
â”œâ”€â”€ genome_annotation.gff3   # VP1-2B CDS annotation (shared across all genotypes)
â”œâ”€â”€ sequences.fasta          # all sequences for this genotype (used to build tree)
â”œâ”€â”€ pathogen.json            # dataset config with clade-like attributes
â”œâ”€â”€ tree.json                # Auspice JSON reference tree (built with augur)
â”œâ”€â”€ README.md
â””â”€â”€ CHANGELOG.md
```

###### 4.2.1 â€” Choose a reference sequence per genotype

From `input.fa`, pick one high-quality, complete sequence per genotype to serve as the
alignment reference. Record its key (FASTA header) and add it as `reference.fasta`.

> **Tip:** use the `metadata.tsv` to filter by `genotype` and pick a sequence with a known
> lineage and a good sample date. Prefer sequences that are centrally placed in the clade.

###### 4.2.2 â€” Copy and reuse the genome annotation

The GFF3 annotation is the same for all genotypes since the dataset covers the same
VP1-2B junction fragment. Copy from the existing community dataset:

```bash
DATASET_DIR="data/local_datasets/$(date +%Y-%m-%d)"
for GENOTYPE in IA IB IIA IIB IIIA IIIB; do
  mkdir -p "$DATASET_DIR/$GENOTYPE"
  cp data/nextclade_datasets/hav-vp1-2b/genome_annotation.gff3 \
     "$DATASET_DIR/$GENOTYPE/genome_annotation.gff3"
done
```

###### 4.2.3 â€” Build the reference trees with augur

`scripts/build_trees.sh` runs the full augur pipeline for every genotype in a
single command:

```bash
cd /path/to/hav_dev
bash scripts/build_trees.sh 2026-04-10
```

The script performs seven steps per genotype:

| Step | Tool | Output | Purpose |
|------|------|--------|---------|
| 1 | `nextclade run` | `aligned.fasta` | Align sequences to the genotype reference |
| 2 | Python filter | `aligned_filtered.fasta` | Drop sequences with < 100 non-gap bases |
| 3 | `augur tree` (IQ-TREE) | `tree_raw.nwk` | Maximum-likelihood tree |
| 4 | `augur refine` | `tree.nwk`, `branch_lengths.json` | Midpoint root + optimised branch lengths |
| 5 | `augur ancestral` | `nt_muts.json` | Ancestral mutation reconstruction |
| 6 | `annotate_clades.py` | `clades.json` | Lineage labels on tips **and** internal nodes |
| 7 | `augur export v2` | `tree.json` | Auspice v2 JSON consumed by NextClade |

> **Why midpoint rooting?** Many sequences in the local database lack
> collection dates, so time-tree inference is unreliable. Midpoint rooting
> gives a stable, reproducible root without requiring dates.

##### How clade/lineage assignment works

NextClade assigns a lineage to a query sequence by:

1. Aligning it to the reference.
2. Placing it on the reference tree using parsimony.
3. Walking up the tree to the nearest ancestor with a `clade_membership` label.
4. Reporting that label as the `clade` value in the output.

For this to work correctly, **internal nodes** of the reference tree must carry
`clade_membership` labels â€” not just the tips.

`scripts/annotate_clades.py` handles this via a strict monophyletic propagation
rule: an internal node is labelled if and only if **all** of its leaf descendants
share the same non-empty lineage label.  Sequences with missing or unknown
lineage values block propagation through their ancestors, preventing incorrect
assignments.

```
         A (V1)
        /
  NODE_1 (V1)        â† all descendants are V1 â†’ labelled
        \
         B (V1)

         C (V1)
        /
  NODE_2               â† mixed descendants â†’ NOT labelled
        \
         D (V2)
```

Queries placed in the subtree under `NODE_1` will be assigned lineage V1.
Queries placed near `NODE_2` will walk further up until they find a labelled
ancestor (or return no assignment if there is none).

##### Defining new lineages

To introduce a new lineage label:

1. Identify the sequences that belong to the new variant using the tree or
   epidemiological metadata.
2. In `data/local_datasets/<date>/metadata.tsv` (or the per-genotype copy),
   set the `lineage` column for those sequences to a consistent name, e.g.
   `NOR-2026-V1`.
3. Re-run the tree build step:
   ```bash
   bash scripts/build_trees.sh <date>
   ```
4. `annotate_clades.py` will automatically propagate the new label to all
   monophyletic internal nodes in that subtree.

> **Naming convention:** use `{Country}-{Year}-{Variant}{Sublabel}` where
> practical (e.g. `NOR-2026-V1`, `NOR-2026-V1a`). Check existing lineage names
> in `metadata.tsv` for consistency.

###### 4.2.4 â€” Create `pathogen.json`

Each genotype dataset needs its own `pathogen.json`. The `clade_membership` key maps to
`lineage` so Nextclade reports the outbreak variant in its output:

```json
{
  "$schema": "https://raw.githubusercontent.com/nextstrain/nextclade/refs/heads/release/packages/nextclade-schemas/input-pathogen-json.schema.json",
  "schemaVersion": "3.0.0",
  "alignmentParams": {
    "alignmentPreset": "high-diversity",
    "minSeedCover": 0
  },
  "files": {
    "reference": "reference.fasta",
    "pathogenJson": "pathogen.json",
    "treeJson": "tree.json",
    "genomeAnnotation": "genome_annotation.gff3",
    "examples": "sequences.fasta",
    "readme": "README.md",
    "changelog": "CHANGELOG.md"
  },
  "attributes": {
    "name": "Hepatitis A virus â€” genotype IA",
    "reference name": "HAV genotype IA local reference",
    "reference accession": "<KEY_OF_REFERENCE_SEQUENCE>"
  },
  "qc": {
    "missingData": { "enabled": true, "missingDataThreshold": 100, "scoreBias": 50 },
    "snpClusters": { "enabled": true, "windowSize": 100, "clusterCutOff": 6 },
    "mixedSites": { "enabled": true, "mixedSitesThreshold": 10 },
    "privateMutations": { "enabled": true, "typical": 5, "cutoff": 50 },
    "frameShifts": { "enabled": true },
    "stopCodons": { "enabled": true }
  }
}
```
NB! In the final `clades.json` and `tree.json` nodes that have more than one descendant variant lineage are currently only named according to the node name. This should be manually inspected and in the future a hierarchical lineage naming system should be implemented. 

##### 4.2.5 â€” Validate the dataset

```bash
# Check annotation is read correctly
nextclade read-annotation "$DIR/genome_annotation.gff3"

# Test-run on example sequences
nextclade run \
  --input-dataset "$DIR" \
  --output-tsv /tmp/test_$GENOTYPE.tsv \
  "$DIR/sequences.fasta"
```

### 5 â€” Run NextClade with the local dataset

After Step 3, use the `clade` column to route each sample to its genotype-specific dataset.

#### 5.1 â€” Split batch sequences by genotype

```bash
BATCH="data/Batch-1"
STEP3_TSV="$BATCH/output/nextclade.tsv"
DATASET_DIR="data/local_datasets/$(date +%Y-%m-%d)"

# Extract genotype assignment per sample from Step 3 output
for GENOTYPE in IA IB IIA IIB IIIA IIIB; do
  # Get list of seqNames assigned to this genotype
  awk -F'\t' -v g="$GENOTYPE" 'NR>1 && $3==g {print $1}' "$STEP3_TSV" \
    > /tmp/ids_$GENOTYPE.txt

  # Extract matching sequences from the batch FASTA
  seqkit grep -f /tmp/ids_$GENOTYPE.txt "$BATCH/Batch-1.fa" \
    > "$BATCH/output/Batch-1_$GENOTYPE.fa"
done
```

> `seqkit` is available in the HAVDEV conda environment. The column index of `clade` in the
> TSV may vary â€” adjust `$3` if needed, or use `awk -F'\t' -v c="clade"`.

#### 5.2 â€” Run Nextclade per genotype

```bash
for GENOTYPE in IA IB IIA IIB IIIA IIIB; do
  INPUT="$BATCH/output/Batch-1_$GENOTYPE.fa"
  [ -s "$INPUT" ] || continue   # skip if no samples for this genotype

  mkdir -p "$BATCH/output/local/$GENOTYPE"

  nextclade run \
    --input-dataset "$DATASET_DIR/$GENOTYPE" \
    --output-tsv    "$BATCH/output/local/$GENOTYPE/nextclade.tsv" \
    --output-fasta  "$BATCH/output/local/$GENOTYPE/nextclade.aligned.fasta" \
    --output-json   "$BATCH/output/local/$GENOTYPE/nextclade.json" \
    "$INPUT"
done
```

### Step 6 â€” Review results

The key output columns when using the local database:

| Column | Description |
|--------|-------------|
| `seqName` | Sample identifier |
| `clade` | Outbreak variant / lineage (e.g. `NOR-2024-V1a`) assigned via tree placement |
| `qc.overallStatus` | Overall quality: `good`, `mediocre`, `bad` |
| `qc.overallScore` | Numeric QC score (lower is better; 0 = perfect) |
| `totalSubstitutions` | Nucleotide mutations vs. local genotype reference |
| `totalInsertions` | Inserted bases vs. reference |

> **Note:** the `clade` column in the local NextClade output corresponds to the
> lineage/variant label (e.g. `NOR-2024-V1a`), not the HAV genotype.  The
> genotype is implicit from which local dataset the sequence was assigned to in
> Step 3.  Sequences without a labelled nearest ancestor will have an empty
> `clade` value â€” this indicates an unclassified variant or a sequence that
> failed tree placement.

Future developments:
Make scripts that compares the first NextClade run (with the built-in dataset) to the second run (with the local dataset).

---

## Step 7 â€” Updating local datasets

Run the full update pipeline whenever new sequences are added to the source
database, or when you want to assign lineage labels to previously unclassified
sequences.

### 7.1 â€” Export new data from the sequence database

Export two files from the database and place them in a dated directory:

```
data/local_datasets/YYYY-MM-DD/
â”œâ”€â”€ export.csv    â† metadata export (semicolon-delimited, Latin-1 encoded)
â””â”€â”€ 2PA.fa        â† FASTA of all sequences
```

### 5.2 â€” Run the update script

```bash
cd /path/to/hav_dev
conda activate /path/to/hav_dev/.conda/HAVDEV
bash scripts/update_datasets.sh 2026-05-01
```

This runs the full pipeline:

| Step | Script | What it does |
|------|--------|--------------|
| 1 | `prepare_dataset.R` | Cleans FASTA headers and produces UTF-8 `metadata.tsv` |
| 2 | `build_local_datasets.R` | Creates per-genotype dirs with reference, sequences, annotation, pathogen.json |
| 3 | `build_trees.sh` | Rebuilds all four genotype trees, including internal node clade annotation |
| 4 | Validation | Runs nextclade on the reference sequences and reports QC pass rate |

### 5.3 â€” Assign lineages to new sequences

After the initial build, inspect sequences that have an empty `lineage` in
`data/local_datasets/<date>/<genotype>/metadata.tsv`:

```bash
# List unclassified sequences for genotype Ia
awk -F'\t' 'NR>1 && $3==""' data/local_datasets/2026-05-01/Ia/metadata.tsv | cut -f1
```

If a cluster of new sequences forms a new variant:

1. Assign a lineage name in the **master** metadata file:
   ```
   data/local_datasets/2026-05-01/metadata.tsv
   ```
   (The per-genotype copies are regenerated from this file by `build_local_datasets.R`.)

2. Re-run just the tree step:
   ```bash
   bash scripts/build_trees.sh 2026-05-01
   ```

### 5.4 â€” Naming new lineages

Use the convention `{Country}-{Year}-V{Number}{Sublabel}`:

| Example | Meaning |
|---------|---------|
| `NOR-2026-V1` | First new Norwegian cluster in 2026 |
| `NOR-2026-V1a` | Sublineage of V1 |

Check existing names in `metadata.tsv` to avoid conflicts.  A new lineage is
warranted when a cluster is:
- Monophyletic (forms its own clade in the tree)
- Epidemiologically distinct
- Differs from related lineages by at least one informative substitution

## Output Files Explained

- **nextclade.tsv** - Tab-separated results with clades and QC metrics
- **nextclade.aligned.fasta** - Sequences aligned to the HAV VP1-2B reference
- **nextclade.json** - Detailed JSON for programmatic access

## Tips

- Use `--help` on any command for full options
- Results are reproducible â€” same input always gives same output
- Quality control: check `qc.overallStatus` and `alignmentScore` columns
- If a sequence fails QC, inspect `qc.missingData` and `qc.privateMutations`

## Resources

- NextClade documentation: https://docs.nextstrain.org/projects/nextclade/
- HAV dataset: https://github.com/nextstrain/nextclade_data
- GitHub: https://github.com/nextstrain/nextclade
