---
title: "A Whole-Genome Classification and Nomenclature System for Hepatitis A Virus"
subtitle: "Harmonizing HAV assignments across laboratories and countries"
author: "Norwegian Institute of Public Health (FHI) · EU Reference Laboratory for HAV"
date: "2026-06-05"
---

<style>
:root{
  --ink:#1f2933; --muted:#52606d; --accent:#1d4e89; --accent2:#137a63;
  --bg:#ffffff; --soft:#f4f6f8; --line:#d9e2ec; --warn:#b54708; --ok:#137a63;
}
html{ -webkit-text-size-adjust:100%; }
body{
  font-family: -apple-system, "Segoe UI", Roboto, Helvetica, Arial, sans-serif;
  color:var(--ink); line-height:1.55; max-width:980px; margin:2rem auto; padding:0 1.4rem;
  background:var(--bg);
}
h1{ font-size:2.0rem; line-height:1.2; color:var(--accent); border-bottom:3px solid var(--accent); padding-bottom:.4rem; }
h2{ font-size:1.4rem; color:var(--accent); margin-top:2.2rem; border-bottom:1px solid var(--line); padding-bottom:.25rem; }
h3{ font-size:1.12rem; color:var(--accent2); margin-top:1.4rem; }
.subtitle{ font-size:1.15rem; color:var(--muted); font-weight:500; }
.author, .date{ color:var(--muted); }
code{ background:var(--soft); padding:.1rem .35rem; border-radius:4px; font-size:.9em; }
pre code{ display:block; padding:.8rem 1rem; overflow:auto; }
table{ border-collapse:collapse; width:100%; margin:1rem 0; font-size:.95rem; }
th,td{ border:1px solid var(--line); padding:.5rem .65rem; text-align:left; vertical-align:top; }
th{ background:var(--soft); color:var(--ink); }
tr:nth-child(even) td{ background:#fbfcfd; }
blockquote{ border-left:4px solid var(--accent2); margin:1rem 0; padding:.4rem 1rem; background:var(--soft); color:var(--muted); }
#TOC{ background:var(--soft); border:1px solid var(--line); border-radius:8px; padding:1rem 1.4rem; margin:1.5rem 0; }
#TOC ul{ margin:.2rem 0; }
.callout{ border:1px solid var(--line); border-left:5px solid var(--accent); background:var(--soft); border-radius:6px; padding:.8rem 1rem; margin:1.2rem 0; }
.callout.warn{ border-left-color:var(--warn); }
.callout.ok{ border-left-color:var(--ok); }
.tag{ display:inline-block; font-size:.75rem; font-weight:600; padding:.1rem .5rem; border-radius:999px; background:var(--accent); color:#fff; vertical-align:middle; }
.layer-genotype{ background:#1d4e89; } .layer-lineage{ background:#137a63; } .layer-outbreak{ background:#b54708; }
hr{ border:none; border-top:1px solid var(--line); margin:2rem 0; }
.small{ font-size:.85rem; color:var(--muted); }
</style>

# A Whole-Genome Classification and Nomenclature System for Hepatitis A Virus

<p class="subtitle">Harmonizing HAV assignments across laboratories and countries, and evaluating a new <em>lineage</em> layer between genotype/clade and outbreak.</p>

<p class="small">Working planning document · v0.1 · 2026-06-05 · FHI &amp; EU Reference Laboratory for HAV</p>

---

## 1. Executive summary

Hepatitis A virus (HAV) genotyping in routine surveillance still relies on a **short subgenomic fragment** (typically the VP1/P2A or VP1/2A junction). This supports only two practical layers of classification: a **coarse genotype/clade** (IA, IB, IIA, IIB, IIIA, IIIB) and **ad-hoc outbreak labels** that are locally meaningful but **not comparable between laboratories or countries** and carry no phylogenetic meaning.

We are now moving to **whole-genome sequencing (WGS)** via a tiling-PCR protocol. WGS dramatically increases resolution but also breaks the assumptions behind current definitions: an "outbreak = identical sequence" rule is almost certainly **too strict at whole-genome resolution**, and the boundaries of genotypes/clades — and whether an **intermediate "lineage" layer is warranted** — have never been evaluated on full genomes.

This document proposes a staged programme to build a **WGS-based, mutation-defined, tree-anchored classification system** for HAV that can be **assigned automatically with Nextclade** and **harmonized across the EURL network**. It deliberately keeps the new **"lineage" layer as a hypothesis to be tested**, not a foregone conclusion, and defines a concrete **8-month contribution for an incoming master student** whose clustering work feeds directly into the key go/no-go decisions.

Two recent nomenclature efforts frame the design space and are used as templates throughout, with a third (SARS-CoV-2 Pango) as the closest *structural* analog:

- **Consortium-curated model** — Fusaro et al., *Proposal for a Global Classification and Nomenclature System for A/H9 Influenza Viruses*, Emerg Infect Dis 2024 (PMID 39043566). An international consortium of 22 laboratories, **led from the EU Reference Laboratory for Avian Influenza** and endorsed by the WOAH/FAO network, defined 3 lineages and 36 clades from 10,638 HA sequences using explicit, reusable criteria (monophyly, ≥3 sequences sampled over >3 years, bootstrap >80%, ≥1 shared amino-acid mutation, ~6% between-clade average pairwise distance, **APD**), with a public assignment tool. **This is a direct governance precedent for an EURL-led HAV consortium.**
- **Dynamic / algorithmic model** — Neher, Huddleston, Bedford et al., *Nomenclature for Tracking of Genetic Variation of Seasonal Influenza Viruses*, Influenza Other Respir Viruses 2026 (doi:10.1111/irv.70230). A semi-automatic algorithm proposes subclades from a Nextstrain reference tree by combining three scores — phylogenetic "bushiness" (size), nucleotide divergence since the parent breakpoint, and weighted amino-acid substitutions — which experts then ratify. Names are hierarchical and Pango-aliased (e.g. `G.1.3`, alias `K`), definitions are machine-readable on GitHub, and assignment is **built into Nextclade** and used routinely by GISRS.
- **Closest structural analog — SARS-CoV-2 Pango** (Rambaut et al. 2020; `autolin`, McBroome et al. 2024). Unlike the two influenza systems, which classify a *single gene/segment*, Pango is a **whole-genome, non-segmented, mutation-defined hierarchical** scheme — exactly HAV's situation. HAV should borrow Pango's whole-genome hierarchical structure, the H9 consortium's governance, and the seasonal-influenza algorithmic designation step.

---

## 2. Problem statement

### 2.1 The resolution and comparability gap

| Layer | Current basis | Strength | Weakness |
|---|---|---|---|
| **Genotype / clade** <span class="tag layer-genotype">existing</span> | Short fragment, well-established (IA…IIIB) | Stable, globally agreed, decades of data | Too coarse for outbreak attribution; says little about transmission |
| **Outbreak label** <span class="tag layer-outbreak">existing</span> | Often "identical sequence" on a short fragment; named by geography + year | Operationally useful inside one lab | Lab/country-specific names, **not comparable**, no phylogenetic meaning, definition breaks under WGS |
| **Lineage** <span class="tag layer-lineage">proposed?</span> | — | Could bridge the two (cf. SARS-CoV-2 Pango, influenza, mpox) | Does not yet exist for HAV; must be shown to be real and useful |

### 2.2 Specific problems to solve

1. **No comparability between labs.** Outbreak names encode local context, not shared, reproducible genetic definitions. The same transmission chain can have different names in different countries.
2. **No phylogenetic meaning in outbreak labels.** They cannot be placed on a tree, versioned, or reasoned about evolutionarily.
3. **Undefined WGS thresholds.** Moving from ~300–500 nt to ~7.5 kb changes the genetic-distance scale. "Identical sequence" is likely too strict; we have no agreed SNP/distance cutoff, and no estimate of within- vs between-outbreak diversity on full genomes.
4. **Undefined genotype/clade boundaries on WGS.** Clade definitions were built on a fragment; their stability and monophyly on whole genomes are untested. Recombination, if present, must be screened.
5. **A missing intermediate layer.** There is currently nothing between "clade" (too broad) and "outbreak" (too narrow). Whether HAV exhibits stable, epidemiologically meaningful intermediate structure is an open, answerable question.
6. **Sparse WGS reference data.** A Nextclade-style assignment needs a dense, representative, mutation-annotated reference tree. For HAV the short-fragment tree is good, but **whole-genome data are currently too sparse** to build one robustly.
7. **No governance / designation process.** Even a perfect technical scheme fails without an agreed naming convention, an open designation repository, versioning, and a body that ratifies changes.

<div class="callout warn">
<strong>Core tension.</strong> WGS gives us the resolution to define harmonized lineages and rational outbreak thresholds — but the reference data needed to operationalize that resolution does not yet exist at scale. The programme below is sequenced to close that gap before committing to fixed definitions.
</div>

---

## 3. Goal scenario

A HAV laboratory anywhere in the EURL network sequences an isolate (full genome where possible, fragment as fallback), runs a **single open tool (Nextclade)** against a **shared HAV reference dataset**, and receives a **comparable, versioned assignment** at up to three layers:

1. **Genotype / clade** <span class="tag layer-genotype">retained</span> — backward-compatible with the existing IA…IIIB system.
2. **Lineage** <span class="tag layer-lineage">if validated</span> — a **mutation-defined, hierarchical, tree-anchored** label (e.g. `IA.3.2`), stable over time, defined in a public repository, assignable from WGS and — where the signal exists — from the short fragment too.
3. **Outbreak** <span class="tag layer-outbreak">harmonized</span> — a **standardized, threshold-based definition** for WGS (genetic-distance cutoff + temporal/epidemiological metadata) that replaces "identical sequence" and produces the **same cluster regardless of which lab runs it**.

Supporting properties of the target system:

- **Tree-anchored & mutation-defined.** Every lineage corresponds to a clade on a mutation-annotated reference phylogeny (Nextstrain/`augur`), defined by a set of signature substitutions — so it is transparent, falsifiable, and re-assignable as data grow.
- **Open and reproducible.** Reference alignment, reference tree, lineage definitions, and assignment tool are public and versioned (GitHub, à la `pango-designation`).
- **Backward-compatible & fragment-aware.** Existing clade calls are preserved; the system degrades gracefully when only the classical fragment is available.
- **Governed by the EURL.** A defined designation/ratification process, change log, and versioning, drawing on EURL expertise and reference materials.
- **Integrated into the existing `hav_dev` pipeline** (Nextclade against community datasets + BLAST against internal datasets).

<div class="callout ok">
<strong>Decision gate built into the goal.</strong> The "lineage" layer is included <em>conditionally</em>. If WGS analysis shows stable, monophyletic, epidemiologically meaningful intermediate structure → we build and govern the lineage layer. If not → we drop it and invest only in <strong>harmonizing the outbreak layer</strong> with WGS thresholds. The plan funds the evidence to make this call (§6, §7).
</div>

---

## 4. Reference models for the lineage layer

| Dimension | Consortium-curated (H9 influenza, EID 2024) | Dynamic / algorithmic (seasonal influenza, IRV 2026) | Implication for HAV |
|---|---|---|---|
| Who defines lineages | Expert consortium consensus | Algorithm proposes; experts confirm | EURL is a natural consortium **and** can run an algorithm → **hybrid** |
| Definition basis | Monophyly + bootstrap + within/between APD + shared AA + epidemiology | Combined score: tree "bushiness" + nucleotide divergence + weighted AA substitutions | Both map cleanly onto a mutation-annotated tree |
| Clustering tool | PhyCLIP on an IQ-TREE ML tree | LBI-like algorithm on a Nextstrain tree (cf. `autolin`) | Either works; pilot both on HAV WGS |
| Naming | Lineage letter + ordinal clades (`Y1.1.1`); `X-like` for unassigned | Hierarchical, Pango-aliased labels (`G.1.3` → `K`) | Hierarchical, extensible, alias-capable naming preferred |
| Tooling | Dedicated online tool (rebuilds ML tree per query) | **Built into Nextclade** (web + CLI; 100k seqs/min) | Nextclade is already in `hav_dev` → low-friction path |
| Update cadence | Periodic expert review | ~4×/year, version-logged on GitHub | HAV evolves slowly → periodic review is realistic |

<div class="callout">
<strong>Recommended HAV approach: a hybrid.</strong> Use the <em>algorithmic</em> machinery (score-based proposals on a mutation-annotated tree, hierarchical Pango-style names, Nextclade integration) to <em>propose</em> candidate lineages, and the <em>consortium</em> machinery (EURL expert ratification, public GitHub designation repo with machine-readable definitions) to <em>confirm, name, and govern</em> them. This matches HAV's slow clock and small expert community.
</div>

### 4.1 Concrete criteria worth borrowing (and recalibrating)

<div class="callout">
<strong>APD — average pairwise (nucleotide) distance.</strong> For a set of sequences, take every pair, compute the proportion of nucleotide sites at which they differ (a <em>p-distance</em>), and average over all pairs (the H9 study computed this in MEGA X). <strong>Within-clade APD</strong> measures internal diversity (small → a tight, homogeneous group); <strong>between-clade APD</strong> measures separation between groups (large → well-separated clades). The H9 system used it as a boundary criterion (within-clade &lt; 6%, recommended between-clade &gt; 6%, with flexibility). In this plan, APD is the metric for the WGS outbreak-threshold study (within- vs between-outbreak distances) and for characterising candidate lineages — note the cut-offs must be re-measured for HAV, not copied from influenza.
</div>

The H9 system's **clade-definition criteria** are directly reusable as a starting template, *recalibrated for HAV*: a candidate group should be **monophyletic**, contain **≥3 sequences**, be **sampled over >1 year / multiple seasons** (to distinguish a sustained lineage from a point-source outbreak), have **strong nodal support**, and share **≥1 defining substitution**. The seasonal-influenza system contributes the **machine-readable definition format** (each clade = name + parent + defining substitutions, optionally representative sequences) and the **objective proposal score** (size + divergence + weighted substitutions, with a minimum group size). HAV-specific thresholds (distance cut-offs, minimum size, time span) must be *measured*, not copied — that is core master-project work (§7).

### 4.2 Why HAV is not influenza — and why that helps

<div class="callout ok">
HAV is a <strong>non-segmented, single-stranded positive-sense RNA picornavirus</strong> with a <strong>slow substitution rate</strong> and <strong>no reassortment</strong>. This is materially simpler than influenza:
<ul>
<li><strong>One genome, one tree.</strong> No per-segment nomenclature (HA vs NA); a single whole-genome phylogeny suffices — exactly the SARS-CoV-2/Pango setting.</li>
<li><strong>Recombination is rare but real.</strong> It must be screened (it is the main thing that can break a single-tree scheme), and definitions may need to be region-aware if detected.</li>
<li><strong>Slow clock = stable lineages, but weak short-term signal.</strong> Lineages will be durable (good for harmonization), but algorithmic thresholds tuned to fast-evolving influenza must be lowered, and time-scaled inference may be limited — favouring <em>substitution-based</em> over <em>time-based</em> definitions.</li>
<li><strong>Existing genotype thresholds give an anchor.</strong> Classical HAV typing places sub-genotypes (IA/IB…) at ~7.5% and genotypes at ~15% nucleotide divergence in the VP1/2A fragment. The new <em>lineage</em> layer sits <strong>below the ~7.5% sub-genotype level</strong>, and the <em>outbreak</em> layer far below that — the WGS distance scale for both must be established empirically.</li>
</ul>
</div>

---

## 5. What do we need

### 5.1 Data
- A **dense, representative whole-genome dataset**: geographically and temporally broad, spanning the relevant clades (especially IA/IB/IIIA, dominant in Europe).
    - The **master student at NIPH can create whole-genome data from already-sequenced strains** in the lab's collection, directly densifying this dataset.
- **Reference materials / characterized strains** to anchor and validate definitions.
- **Public whole genome sequences** (GenBank, ENA, etc.) curated with metadata (date, country, source, clinical/outbreak context).
- Paired **fragment + full-genome** sequences for the same isolates, to map the fragment world onto the genome world.

### 5.2 Methods & artefacts
- A **curated reference alignment** (MAFFT) and a **recombination-screened** (e.g. RDP/3SEQ/GARD), **mutation-annotated reference phylogeny** (IQ-TREE + Nextstrain `augur`/TreeTime).
- A **molecular-clock / temporal-signal** assessment (TreeTime, root-to-tip regression) to know whether time-scaled structure is reliable for HAV's slow clock.
- A **lineage-definition methodology**: a clustering step (**PhyCLIP** and/or the `autolin`-style score: size + divergence + weighted AA substitutions) plus recalibrated thresholds for between/within distance, minimum support/size, and shared substitutions; monophyly and stability criteria; the hybrid algorithm-then-expert process.
- A **naming convention** (Pango-style hierarchical + aliasing) and an open **designation repository** (GitHub, machine-readable `name + parent + defining substitutions` per lineage, versioned change log).
- A **Nextclade dataset** for HAV WGS (reference, root, clade/lineage definitions, QC rules) plus a **fragment-compatible** variant for labs still on short-fragment typing.
    - **Reference-architecture decision (must be tested, not assumed).** Nextclade aligns and calls mutations against a *single* reference per dataset; alignment degrades as query–reference divergence rises. Given HAV's cross-genotype divergence (>15% nt), the right architecture is an open question — like influenza (per-subtype datasets), not SARS-CoV-2 (one global reference). Empirically compare candidate references (per-genotype representatives, a central strain, and a **reconstructed ancestral/consensus root**) and compare a **single global dataset** vs **per-genotype datasets with an upstream router** (the existing BLAST step, or Nextclade 3's minimizer auto-selection). See §7.2(6).
- **Outbreak-threshold study**: empirical within-outbreak vs between-outbreak pairwise-distance distributions on WGS (cf. H9's within/between APD analysis) → a defensible cutoff (likely distance **+** temporal/epi metadata, not distance alone).

### 5.3 People, governance & infrastructure
- **EURL governance**: a working group that ratifies designations, owns versioning, and publishes the change log.
- **Compute & pipeline**: integration of the new reference dataset into `hav_dev` (Nextclade + BLAST), CI to rebuild the tree/dataset as data grow.
- **Validation set**: known, well-characterized outbreaks to test that the scheme recovers the right clusters.
- **Dissemination**: a methods/consortium paper (templated on the two reference papers) and EURL protocol documentation.

### 5.4 Readiness checklist
| Need | Status today | Gating action |
|---|---|---|
| Short-fragment reference tree | ✅ Good | Maintain |
| WGS tiling-PCR protocol | ✅ Developed | Scale sequencing |
| Dense WGS reference data | ❌ Too sparse | **Master project + EURL submissions** |
| WGS outbreak thresholds | ❌ Undefined | **Within/between-distance study** |
| Lineage go/no-go evidence | ❌ Missing | **Fragment-vs-WGS clustering comparison** |
| Nextclade WGS dataset | ❌ Not built | Build after tree is dense enough |
| Governance / naming | ❌ Not established | EURL working group |

---

## 6. Phased roadmap

> Phases overlap; the master project (§7) sits inside Phases 1–2 and de-risks the Phase-3 decision.

### Phase 0 — Framing & governance (months 0–3)
- Discuss in the EURL group; agree scope, success criteria, and the lineage decision gate.
- Adopt data standards (metadata schema, sequence QC, submission route).
- Lock the two reference papers as design templates.

### Phase 1 — Build the WGS evidence base (months 1–10)
- Sequence the selected strain panel (master student + EURL).
- Aggregate and curate public + partner WGS with metadata.
- Build the recombination-screened, mutation-annotated reference tree; assess temporal signal.

### Phase 2 — Characterize structure (months 4–12)
- **Fragment vs full-genome clustering comparison** (resolution gain, congruence/incongruence).
- **Within- vs between-outbreak distance distributions** → candidate WGS outbreak threshold.
- Detect candidate intermediate (lineage-level) structure; test monophyly and stability.

### Phase 3 — Decision gate: lineage layer? (≈ month 12)
- **If** stable, meaningful intermediate structure exists → define candidate lineages (hybrid method), draft naming convention, prototype Nextclade dataset.
- **If not** → finalize a harmonized **outbreak-threshold** definition only; document the negative result.

### Phase 4 — Operationalize & harmonize (months 12–24)
- Publish the reference dataset, definitions, and designation repo.
- Release the Nextclade dataset(s); integrate into `hav_dev`.
- Validate on known outbreaks; run an EURL inter-lab concordance exercise.
- Publish the methods/consortium paper; ratify governance and versioning.

---

## 7. Master-student contribution (8 months)

The student's work is the **engine of the evidence base** and feeds directly into the Phase-3 decision. It is scoped to be a coherent, publishable thesis on its own.

### 7.1 Wet-lab
- Whole-genome tiling-PCR sequencing of **selected HAV strains** (prioritized from the EURL panel and Norwegian/European surveillance), expanding the sparse WGS dataset.

### 7.2 Bioinformatics (core scientific questions)
1. **Fragment vs whole-genome clustering.** Quantify how much resolution WGS adds: extract the classical VP1/2A fragment *in silico* from each genome, cluster both, and measure where fragment-based clusters split or merge under WGS (e.g. tanglegrams, congruence metrics). → directly informs whether a lineage layer is justified.
2. **Outbreak thresholds on WGS.** Estimate within-outbreak vs between-outbreak pairwise genetic distances on full genomes (using epidemiologically confirmed outbreaks); propose a candidate cutoff to replace "identical sequence."
3. **Recombination & tree quality.** Screen for recombination (RDP/3SEQ/GARD); build a clean MAFFT alignment and IQ-TREE/`augur` tree; assess temporal signal / clock (TreeTime, root-to-tip).
4. **Pilot lineage definitions.** For one or two dominant sub-genotypes (e.g. IA/IB), apply PhyCLIP and/or an `autolin`-style score to test whether stable, monophyletic intermediate clusters exist and can be mutation-defined; characterise within/between distances per the H9 criteria.
5. **Prototype assignment.** Build a draft Nextclade dataset (reference + annotated tree + lineage definitions) from the pilot tree and benchmark assignment accuracy, including fragment-only fallback.
6. **Reference & dataset-architecture benchmarking.** Test whether a single Nextclade reference is appropriate for HAV or whether per-genotype references/datasets are needed. Compare candidate references (per-genotype representatives, a central strain, a reconstructed ancestral/consensus root) on per-genotype **alignment coverage, mutation-call accuracy** (vs a curated whole-genome MSA), **QC-flag rates**, and **assignment concordance**; and compare a **single global dataset** against **per-genotype datasets + an upstream router** (BLAST / minimizer auto-selection). Scope the full comparison to well-covered genotypes (IA/IB/IIIA); flag sparse ones as future work.

### 7.3 Suggested 8-month timeline
| Month | Focus |
|---|---|
| 1–2 | Onboarding, protocol, data assembly (public + EURL), reproduce current `hav_dev` pipeline |
| 2–5 | WGS lab work; assembly/QC; curated alignment; recombination screen |
| 4–6 | Reference tree; fragment-vs-WGS clustering comparison |
| 5–7 | Within/between-outbreak distance study → candidate thresholds; pilot lineage definitions |
| 7–8 | Prototype Nextclade dataset; benchmark; thesis write-up; draft contribution to consortium paper |

### 7.4 Deliverables
- MSc thesis.
- A curated, metadata-rich HAV **WGS dataset** (a community asset).
- A **fragment-vs-WGS clustering** analysis and a **candidate outbreak-distance threshold**.
- A **pilot lineage proposal** + prototype Nextclade dataset for ≥1 clade.
- Reproducible code/workflow integrated with `hav_dev`; likely **co-authorship** on the methods/consortium paper.

<div class="callout">
<strong>Scoping note.</strong> The student delivers the <em>evidence and prototypes</em>, not the final ratified, network-wide system — that is the EURL working group's multi-year remit. This keeps the thesis achievable while making it pivotal.
</div>

---

## 8. Risks & mitigations

| Risk | Impact | Mitigation |
|---|---|---|
| WGS reference data stay too sparse | Can't build robust tree/lineages | EURL coordinated sequencing + public-data curation; pilot on best-covered clades first |
| Recombination distorts a single-tree scheme | Lineage definitions unstable | Explicit recombination screening; consider genome-region-aware definitions |
| Weak temporal signal (slow clock) | Time-scaled structure unreliable | Test clock formally; lean on substitution-based, not time-based, definitions |
| No real intermediate structure | "Lineage" layer unjustified | Decision gate (§6.3) — gracefully fall back to harmonized outbreak thresholds |
| Lab/country adoption friction | System unused | Backward compatibility, single open tool (Nextclade), EURL governance & inter-lab exercise |
| Naming/governance disputes | Fragmentation returns | Open designation repo + versioning + EURL ratification, modeled on the two reference systems |

---

## 9. Open questions for the working group

1. Hierarchical naming scheme — extend clade names (`IA.3.2`) or a parallel namespace?
2. Algorithmic thresholds — adopt the seasonal-influenza distance/AA-substitution machinery directly, or recalibrate for HAV's slower clock?
3. Should lineages be **required** to be callable from the classical fragment, or is WGS-only acceptable with fragment fallback?
4. Where does the designation repository live, and who has merge/ratify authority?
5. Outbreak definition — pure genetic threshold, or threshold **plus** mandatory temporal/epidemiological metadata?
6. Reference architecture — a single global Nextclade reference (ideally a reconstructed ancestral root), or **per-genotype datasets with an upstream router** (BLAST / minimizer)? Decide on the §7.2(6) benchmark.

---

## References

1. Fusaro A, Pu J, Zhou Y, et al.; The International H9 Evolution Consortium. **Proposal for a Global Classification and Nomenclature System for A/H9 Influenza Viruses.** *Emerg Infect Dis.* 2024;30(8). PMID 39043566. doi:10.3201/eid3008.231176 — *consortium-curated lineage/clade model; clade-definition criteria; PhyCLIP/IQ-TREE/TreeTime; online assignment tool.*
2. Neher RA, Huddleston J, Bedford T, et al. **Nomenclature for Tracking of Genetic Variation of Seasonal Influenza Viruses.** *Influenza Other Respir Viruses.* 2026;20(2):e70230. doi:10.1111/irv.70230 — *dynamic algorithmic subclade designation; Pango-style aliasing; Nextclade integration; GISRS adoption.*
3. Rambaut A, Holmes EC, O'Toole Á, et al. **A dynamic nomenclature proposal for SARS-CoV-2 lineages to assist genomic epidemiology.** *Nat Microbiol.* 2020;5:1403–7. — *closest structural analog: whole-genome, non-segmented, hierarchical, mutation-defined.*
4. McBroome J, de Bernardi Schneider A, Roemer C, et al. **A framework for automated scalable designation of viral pathogen lineages from genomic data** (`autolin`). *Nat Microbiol.* 2024;9:550–60.
5. Aksamentov I, Roemer C, Hodcroft EB, Neher RA. **Nextclade: clade assignment, mutation calling and quality control for viral genomes.** *J Open Source Softw.* 2021;6(67):3773.

<p class="small">Document v0.1 — for discussion within FHI and the EU Reference Laboratory for HAV. Generated 2026-06-05.</p>
