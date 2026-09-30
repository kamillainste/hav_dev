# hav-consensus

Whole-genome consensus sequences for hepatitis A virus from Oxford Nanopore
amplicon reads (HAVNET WGS primer scheme, genotypes IA, IB and IIIA).

```
raw reads ─► NanoPlot ─► chopper (Q10, ≥1000 bp) ─► NanoPlot ─► ViroConstrictor 1.6.6 ─► consensus + coverage plots
```

For each sample the pipeline writes a consensus FASTA, coverage plots, a summary
row, and all of ViroConstrictor's intermediate files (alignments, trimmed reads,
BAM, VCF, per-position depth). ViroConstrictor runs with two changes from its
defaults, chosen in an internal benchmark. See [docs/METHOD.md](docs/METHOD.md).

## Contents

1. [Install (Windows + WSL + VS Code)](#1-install-windows--wsl--vs-code)
2. [Make a samplesheet](#2-make-a-samplesheet)
3. [Run](#3-run)
4. [Output](#4-output)
5. [Before you report results](#5-before-you-report-results)
6. [Troubleshooting](#6-troubleshooting)

---

## 1. Install (Windows + WSL + VS Code)

Everything runs inside WSL (Linux on Windows). You need about 8 GB of free disk
space in WSL: roughly 5 GB for the conda environments and 2.3 GB for the
containers.

**a. WSL with Ubuntu.** In PowerShell *as administrator*:

```powershell
wsl --install -d Ubuntu
```

Restart Windows, open "Ubuntu" from the Start menu and create a Linux user.
Check that it is WSL **2**: `wsl -l -v` in PowerShell shows `VERSION 2`.

**b. VS Code.** Install the **WSL** extension (by Microsoft). From an Ubuntu
terminal, `code .` opens the current Linux folder in VS Code. Use VS Code's
terminal from then on: it is the Ubuntu shell.

**c. Conda (Miniforge)**, in the Ubuntu terminal:

```bash
wget https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh
bash Miniforge3-Linux-x86_64.sh -b
~/miniforge3/bin/conda init bash
```

Close and reopen the terminal. `conda --version` should now work.

**d. Clone this repository inside WSL**, in your Linux home folder and not
under `/mnt/c/...`:

```bash
cd ~
git clone <repository URL>      # private repo: use `gh auth login` or an SSH key first
cd hav-consensus
```

> Clone with git *from the Ubuntu terminal*. A copy made on the Windows side can
> get Windows line endings, and the scripts then fail (see Troubleshooting).

**e. Run setup once** (20–40 minutes, mostly downloads):

```bash
bash setup.sh
```

This:
- creates four conda environments: `hav_viroconstrictor`, `hav_chopper`,
  `hav_nanoplot` and `hav_plot`
- builds the patched ViroConstrictor in `.build/`
- writes `~/.ViroConstrictor_defaultprofile.ini`, unless you already have one
- downloads ViroConstrictor's containers to `~/.viroconstrictor/containers`
- runs one container as a test

It ends with `Setup complete`. It is safe to re-run.

## 2. Make a samplesheet

Tab-separated, one row per sample. Copy
[example/samplesheet_example.tsv](example/samplesheet_example.tsv) and edit it.
Excel can save it as "Text (Tab delimited)".

| column | required | meaning |
|---|---|---|
| `sample` | yes | sample ID |
| `barcode` | yes | e.g. `barcode01` |
| `run` | yes | name of the run folder; also the output subfolder |
| `path_to_run_folder` | no | folder that *contains* the run folder (see below) |
| `ct`, `gt`, `kommentar` | no | copied into the summary table |

Column order does not matter, and extra columns are ignored.

The pipeline looks for reads in the standard MinKNOW layout:

```
<path_to_run_folder>/<run>/<any>/<any>/fastq_pass/<barcode>/*.fastq.gz
```

It also accepts `<path_to_run_folder>/<run>/<barcode>/*.fastq.gz`.

**Reading data from Windows.** Windows drives appear in WSL under `/mnt/c`,
`/mnt/d` and so on. A network drive has to be mounted first, e.g. for `N:`:

```bash
sudo mkdir -p /mnt/n && sudo mount -t drvfs N: /mnt/n
```

Reading raw reads from there is fine. Write the **output** inside WSL (e.g.
`~/hav_runs/...`), which is much faster. You can open it from Windows Explorer
afterwards: run `explorer.exe .` in the output folder.

## 3. Run

```bash
cd ~/hav-consensus
bash run_hav_pipeline.sh ~/hav_runs/2026-10-01 my_samplesheet.tsv
```

If your samplesheet has no `path_to_run_folder` column, give the folder that
holds the run folders first:

```bash
bash run_hav_pipeline.sh /mnt/n/path/to/runs ~/hav_runs/2026-10-01 my_samplesheet.tsv
```

- Use a **new, empty output folder** for each run, with no spaces in its path.
- Samples are processed one after another. A sample with about 2,000 reads took
  about 1.5 minutes in testing; deep samples take longer. `THREADS=8 bash run_hav_pipeline.sh ...`
  uses more cores (default 6).
- To keep the run going after you close VS Code, run it with
  `nohup bash run_hav_pipeline.sh ... > run.log 2>&1 &` and follow it with
  `tail -f run.log`.

## 4. Output

```
<output_dir>/
  pipeline_run_info.txt           pipeline version, command, samplesheet md5, settings
  summary_table_all_runs.tsv      one row per samplesheet row
  coverage_plots/coverage_<genotype>.png   all samples of a genotype overlaid
  <run>/
    summary_table.tsv
    consensuses/<sample>_<barcode>_<run>.fasta     ◄ the consensus genomes
    coverage_plots/
    <sample>_<barcode>_<run>/
      qc_raw/ qc_filtered/        NanoPlot reports before / after chopper
      filtered/*.fastq.gz         reads after chopper (ViroConstrictor's input)
      coverage_plots/             this sample alone
      results/                    ViroConstrictor results (consensus, MultiQC, mutations)
      data/Virus~HAV/RefID~<ref>/ intermediate files:
        cleaned_fastq/aligned_rawdata/    minimap2 alignment to the reference
        cleaned_fastq/without_adapters/   reads clipped to their aligned part
        cleaned_fastq/QC_filter/          fastp output + report
        cleaned_fastq/without_primers/    primer-trimmed reads (AmpliGone)
        alignment/bam-files/, vcf-files/
        consensus/coverages/              depth per position
        primers/                          amplicon coverage
```

**summary_table.tsv columns:** `raw_reads`, `trimmed_reads` (after chopper),
`avg_depth`, `breadth_of_coverage_pct_at_mincov20` (% of the reference at ≥ 20x),
`reference_id` (the genotype reference chosen), and `status` (`OK`, `SKIPPED: …`
or `FAILED: …`).

Positions below 20x are `N` in the consensus. The scheme's outer primers do
not reach the genome ends, so even a perfect sample reaches only about 98–99%.

## 5. Before you report results

- **Negative controls (NTC).** Check them every run. An NTC can collect enough
  reads from other samples on the same flowcell to produce a partial, mixed
  consensus, often with IUPAC codes. The pipeline does **not** filter this out.
  An NTC with a consensus means the low-depth samples in that run need a closer
  look.
- **Low-depth samples (below ~1,000x mean depth)** are where results are least
  certain. Look at the coverage plot: when the odd amplicons (1, 3, 5 = pool 1)
  or the even ones (2, 4, 6 = pool 2) are all low, one primer pool usually failed
  in the lab.
- **Genotype** is the reference with the best match (`reference_id`), not the
  `gt` column you entered.

## 6. Troubleshooting

| symptom | fix |
|---|---|
| `$'\r': command not found` or `bad interpreter` | Windows line endings: delete the folder and clone again from the Ubuntu terminal |
| `ERROR: … Run setup.sh first` | run `bash setup.sh`, or re-run it if it did not end with `Setup complete` |
| ViroConstrictor asks questions (local/grid, auto-update, …) | `~/.ViroConstrictor_defaultprofile.ini` is missing or incomplete: press Ctrl-C, run `bash setup.sh` |
| setup: `apptainer cannot run …` | make sure it is WSL 2 (`wsl -l -v`), then `wsl --update` in PowerShell, restart, re-run setup |
| setup: container download failed | `ghcr.io` must be reachable (VPN/proxy?). Re-run setup; finished downloads are kept |
| `No space left on device` | free space in WSL (`df -h ~ /tmp`); or point scratch elsewhere: `HAV_TMPDIR=~/tmp bash run_hav_pipeline.sh ...` |
| a conda env is broken | `conda env remove -n hav_<name>`, then `bash setup.sh` |
| a sample says `FAILED: ViroConstrictor error` | look at `<output>/<run>/<sample>/ViroConstrictor_<time>.log` and `logs/` next to it; re-run that sample into a new output folder |

`HAV_ENV_PREFIX` (default `hav_`) changes the env names. Set it the same way
for `setup.sh` and `run_hav_pipeline.sh`.

## Repository layout

```
run_hav_pipeline.sh   the pipeline
setup.sh              one-time installation
envs/                 conda environment definitions
patches/              the change to ViroConstrictor (3 lines)
scripts/              plot_coverage.R, relabel_amplicons.py, apply_patch.py
assets/               HAV reference genomes (IA, IB, IIIA) and HAVNET primers
example/              example samplesheet
docs/METHOD.md        what each step does and why these settings
```

## Credits

ViroConstrictor and AmpliGone (RIVM, AGPL-3.0), chopper and NanoPlot
(W. De Coster), fastp, minimap2, TrueConsense, Snakemake, Apptainer.
`patches/` modifies ViroConstrictor's preset file; the modified copy is built
locally by `setup.sh` and is not redistributed here.
