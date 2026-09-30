# Method

What each step does, and why the settings are what they are. The numbers come
from an internal benchmark (September 2026): 14 samples on three flowcells,
13.8x–35,173x mean depth, all three genotypes, one positive control and one
NTC. The full report (`read_accounting_report.html`) is available from the
pipeline maintainer.

## Steps

| # | step | tool | settings |
|---|---|---|---|
| 1 | read QC | NanoPlot 1.32.1 | raw reads |
| 2 | read filter | chopper 0.12.0 | mean quality ≥ Q10, length ≥ 1,000 bp |
| 3 | read QC | NanoPlot | filtered reads |
| 4 | reference choice | ViroConstrictor `--match-ref` | best of IA (LC435031.1), IB (M14707.1), IIIA (AB279732.1) |
| 5 | alignment | minimap2 `--splice` (HEPATOVIRUS preset) | |
| 6 | clip to aligned part | ViroConstrictor clipper | drops reads with a splice gap > 200 bp |
| 7 | fastp | fastp | **`--cut_right` removed**; only drops reads < 100 bp |
| 8 | primer trimming | AmpliGone 2.0.2 | **`--amplicon-type fragmented`**, `--virtual-primers`, primer mismatch rate 0.1 |
| 9 | consensus | TrueConsense | min depth 20; below that `N` |
| 10 | plots, summary | plot_coverage.R, run_hav_pipeline.sh | |

ViroConstrictor runs each step in its own published container (Apptainer), so
the tool versions inside are fixed by ViroConstrictor 1.6.6.

## The two changes from stock ViroConstrictor

**1. `--amplicon-type fragmented` instead of `end-to-end`.** AmpliGone's two
modes use the same code and differ only in how far from a read end it will look
for a primer to cut at: 10 bp with fragmented, 10,000 bp with end-to-end. At
10,000 bp a read end is cut at the nearest primer downstream, however far
away. In the benchmark, 53.8% of the bases end-to-end removed were not the
read's own primer (19.1% with fragmented). Some reads lost ≥ 500 bp, and some
were cut to nothing.

**2. fastp `--cut_right` removed.** ViroConstrictor's default fastp settings
include a sliding-window quality trim (`--cut_right`, window 20, mean Q7). On
these reads it removed ~254 bp per read, 98% of it sample sequence rather than
primer or adapter. Reads are already quality-filtered by chopper. There is no
command-line option for this setting, so `setup.sh` builds a copy of the
ViroConstrictor package with those three lines blanked
(`patches/viroconstrictor-1.6.6_no_cut_right.diff`). The pipeline puts that
copy first on `PYTHONPATH` and checks it is in use before every run. fastp
still runs, but now only removes reads < 100 bp: 4 reads over all 14 samples.

## Benchmark result

Mean % of the genome called in the consensus (not `N`):

| configuration | ≥ 1,000x (8 samples) | 52–740x (5) | NTC (1) |
|---|---|---|---|
| stock ViroConstrictor (end-to-end, with `--cut_right`) | 97.8 | 79.2 | 32.4 |
| **this pipeline** (fragmented, no `--cut_right`) | 97.9 | 86.1 | 49.4 |
| ARTIC fieldbioinformatics 1.6.0 (for comparison) | 98.0 | 80.6 | 8.2 |

- At ≥ 1,000x every configuration gives the same genome. At 0 to 1 mismatches
  over ~22,000 called positions, ViroConstrictor and ARTIC agree.
- The gains are at low depth, where the consensus is least certain. No sample
  had a known true sequence, so the benchmark shows how much is called, not
  that the extra bases are correct.
- The scheme's outermost primers do not reach the genome ends, so the maximum
  is about 98.3–98.9% depending on genotype.

## Known limitations

- **NTC carry-over.** NTC completeness rises with the same changes that help
  low-depth samples. The NTC's HAV reads come from other samples on the same
  flowcell, and its consensus is a mixture with IUPAC codes. ARTIC's lower NTC
  number comes from how it counts depth (per primer pool), not from detecting
  contamination. No run-level contamination threshold is implemented; check
  NTCs by hand.
- **Length filter.** Chopper's ≥ 1,000 bp filter removes most raw reads in
  low-depth samples: 78% at 52–740x. Reads < 300 bp are primer dimers and
  concatemers. Reads of 300–999 bp are partial HAV amplicons, plus human
  off-target products in some samples. A lower cut-off (300–500 bp) might raise
  completeness at low depth but has not been tested in the pipeline.
- **Human off-target reads.** In some samples up to ~50–60% of filtered reads
  are human PCR products of the HAVNET primers. They do not map to HAV and drop
  out at alignment, but they lower the effective depth.
- **Primer pool failure.** In 5 of the 14 samples primer pool 1 gave only
  1.5–7.5% of amplicon reads. That leaves amplicons 1, 3 and 5 thin. It is a
  lab issue, visible in the coverage plots.
- **Primer sequence at primer sites** in the consensus has not been checked.

## Validation of this repository (2026-09-30)

- **Same result as the benchmark.** A fresh clone, installed with `setup.sh`
  (new conda envs from `envs/lock/`, a new ViroConstrictor profile, containers
  downloaded from ghcr.io), was run on one benchmark sample (genotype IA, 175x,
  2,008 reads after chopper). Its consensus was byte-identical to the
  benchmark's, and its summary row and amplicon coverage table matched.
- **Installation.** `setup.sh` took ~9 minutes on a Linux server, of which
  ~6 minutes was the 2.3 GB container download. It has not yet been run on
  WSL itself.
- **Not yet done:** the full 14-sample comparison with the benchmark.
