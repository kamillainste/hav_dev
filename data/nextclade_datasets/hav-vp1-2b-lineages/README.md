# Nextclade dataset for HAV VP1-2B junction region — with FHI lineage designations

This is a locally enhanced copy of the community `hav-vp1-2b` Nextclade dataset
(MASPHL Bioinformatics, initial release 2025-12-02). It adds a second clade
attribute, **Lineage (phylo)** (`lineage_phylo`), carrying Pango-inspired
sub-genotype lineage designations developed at the Norwegian Institute of Public
Health (FHI). The underlying phylogeny and all broad genotype assignments
(`clade_membership`: IA, IB, IC, IIA, IIB, IIIA, IIIB) are unchanged.

**Public-only derivation.** The lineages — both *which lineages exist* and
*their defining mutations* — are designated **directly on the public community
`hav-vp1-2b` reference tree**. The build has **no dependency on any local/private
FHI dataset**: the community dataset is fetched fresh from the Nextclade dataset
server and lineages are designated on it by phylogenetic criteria alone, so the
whole system is reproducible from public sources.

## Dataset Attributes

| Attribute            | Value                                                                                                            |
| -------------------- | ---------------------------------------------------------------------------------------------------------------- |
| Name                 | Hepatitis A VP1-2B junction sequences (with FHI lineage designations)                                            |
| Upstream dataset     | community/masphl-bioinformatics/hav/vp1-2b-junction                                                              |
| Upstream authors     | Lydia A. Krasilnikova (MASPHL), Mary Godec (MASPHL), Matthew Doucette (MASPHL), Daniel J. Park (Broad Institute) |
| Lineage system       | Norwegian Institute of Public Health (FHI), Folkehelseinstituttet                                                 |
| Data source          | [GenBank](https://www.ncbi.nlm.nih.gov/genbank/) accessed 2025-02-03                                             |
| Reference            | [NC_001489.1](https://www.ncbi.nlm.nih.gov/nuccore/NC_001489.1)                                                  |

## Scope of this dataset

This dataset supports two levels of HAV classification from a single VP1-2B
junction sequence (349 bp, positions 2968–3322 of NC_001489.1):

- **Genotype** (`clade_membership`): IA, IB, IC, IIA, IIB, IIIA, IIIB — assigned
  by the upstream MASPHL tree.
- **Lineage (phylo)** (`lineage_phylo`): Pango-inspired sub-genotype lineages
  (e.g. IA.1, IA.1.2, IB.1.1.1, IIIA.1.1.1) — assigned by the FHI lineage
  system overlaid on the same tree.

Sequences that do not match any defined lineage receive `lineage_phylo =
unassigned`.

## FHI lineage system

### Naming convention

Lineage names follow a Pango-inspired hierarchical convention. The genotype
(IA, IB, IIA, IIIA) forms the root prefix. Subsequent levels are separated by
dots and represent increasingly derived monophyletic clades (e.g. IA → IA.1 →
IA.1.2 → IA.1.2.1). Aliases are introduced if and when names become unwieldy
(not currently required at present depth).

### Lineage register

Lineages are recorded in `data/lineage_designations.tsv` in the project
repository. This register is **regenerated from the public tree on every build**
(it is not hand-curated), so the set below is a snapshot — re-run the build to
refresh it. The currently designated lineages are:

| Lineage     | Genotype | Defining mutations                                            |    N |
|-------------|----------|---------------------------------------------------------------|------|
| IA.1        | IA       | A108G, A179G, A234G, C120T                                    | 1232 |
| IA.1.1      | IA       | C313A, G81A, T138C, T216C                                     | 1209 |
| IA.1.1.1    | IA       | A270G                                                         | 1202 |
| IA.1.1.2    | IA       | A87C                                                          |    4 |
| IA.1.1.3    | IA       | A201G, A45G, G285A, G57A                                      |    3 |
| IA.1.2      | IA       | A12G, G285C                                                   |   23 |
| IA.1.2.1    | IA       | G153A, G240A, T63C                                            |   22 |
| IB.1        | IB       | G231A, G285A, T69C                                            |  411 |
| IB.1.1      | IB       | G99A, T294C                                                   |  299 |
| IB.1.1.1    | IB       | A210G, T138C, T93C                                            |  298 |
| IB.1.2      | IB       | 14 mutations (highly diverged)                                |  112 |
| IB.1.2.1    | IB       | T30C                                                          |  110 |
| IIA.1       | IIA      | 17 mutations (highly diverged)                                |   13 |
| IIA.1.1     | IIA      | T24C, T93C                                                    |   10 |
| IIA.1.1.1   | IIA      | T345C                                                         |    9 |
| IIA.1.2     | IIA      | A108G, G246A, T111C, T223C                                    |    3 |
| IIIA.1      | IIIA     | 14 mutations (highly diverged)                                |   95 |
| IIIA.1.1    | IIIA     | A249G, T237A, T333C                                           |   93 |
| IIIA.1.1.1  | IIIA     | C3T, T33C                                                     |   89 |
| IIIA.1.1.2  | IIIA     | G309A, T297C                                                  |    3 |

### Note on convergent mutations and resolution limits

The VP1-2B junction is only 349 bp, so lineages defined by a single substitution
are inherently low-confidence: the same mutation can occur independently on
multiple branches across the global tree (convergence/homoplasy). Such lineages
may shift or disappear as the public `hav-vp1-2b` tree grows and is re-estimated.
For stable, high-resolution lineage boundaries the whole genome is required (see
the project's WGS classification plan).

A post-processing consistency check is applied during dataset generation: any
`lineage_phylo` label whose genotype prefix does not match the node's broad
`clade_membership` (a cross-genotype convergent misassignment) is replaced with
the parent node's label. Lineages are only designated for genotypes IA, IB, IIA
and IIIA; IC, IIB and IIIB carry too few public sequences and remain
`unassigned`.

## How lineage assignments are generated

The entire process runs on the public community tree only — no local/private
data — and is driven by `scripts/build_vp1b_lineage_dataset.sh`:

1. The community `hav-vp1-2b` dataset is fetched fresh
   (`nextclade dataset get --name community/masphl-bioinformatics/hav/vp1-2b-junction`).
2. `scripts/extract_auspice_tree.py` extracts the newick tree, branch mutations
   (`nt_muts.json`) and a tip metadata TSV (genotype from `clade_membership`,
   date and country) from the community `tree.json`.
3. `scripts/designate_lineages.py` is run **per genotype on the public tree** to
   designate monophyletic lineages by phylogenetic criteria (≥3 sequences and
   ≥1 defining branch mutation; Pango-style hierarchical names). The public tree
   has no outbreak labels, so the outbreak-seeding mode is **not** used. This
   writes `data/lineage_designations.tsv` and an `augur clades` TSV.
4. `augur clades` assigns `lineage_phylo` labels top-down on the same tree.
5. `scripts/patch_vp1b_tree.py` injects the assignments into `tree.json` and
   applies the cascade-clearance and genotype-consistency corrections.
6. `scripts/generate_silo_hierarchy.py` writes the SILO/LAPIS `lineages.yaml`.

To rebuild: `bash scripts/build_vp1b_lineage_dataset.sh` (add `--refresh` to
force re-fetch of the community dataset).



## About the VP1-2B junction region (from upstream dataset)

The HAV VP1-2B junction region covers positions 2,968 through 3,322, 1-indexed, relative to reference [NC_001489.1](https://www.ncbi.nlm.nih.gov/nuccore/NC_001489.1). It is used for analysis including clade assignment by CDC GHOST.

## Upstream dataset documentation

The following sections are reproduced from the upstream MASPHL dataset README
and describe the base phylogeny and genotype assignments that this dataset builds upon.

### Annotations

Note that current Refseq annotations are incorrect for the HAV reference genome as of October 2025. Annotations for this dataset were sourced from Table 1 in [Nainan et al 2006](https://pmc.ncbi.nlm.nih.gov/articles/PMC1360271/#t1).  For this VP1-2B junction dataset, the region includes a portion of both 1D, the complete "2A" (better known now as pX per [Shirasaki et al 2022](https://pmc.ncbi.nlm.nih.gov/articles/PMC9410543/#sec012)), and a portion of 2B. Both 1D and 2B are listed as "truncated" in the annotations. Note that per Nextclade constraints requiring CDS length to be a multiple of 3, the truncated 2B ends at position 348 (1-indexed).

### Sequences included

This dataset includes all human host VP1-2B junction region HAV GenBank sequences. First, all Hepatovirus A GenBank sequences were produced by search term [txid12092\[Organism:exp\]](https://www.ncbi.nlm.nih.gov/nuccore/?term=txid12092[Organism:exp]), accessed 2025-02-03. Sequences were filtered to those fully overlapping the VP1-2B junction region. Sequences were further filtered to human host, using the host field or, when the host field was blank, the isolation_source field and the associated publication(s) in the title field. Sequences for which the host field was blank that were passaged through animals or cell lines were excluded. The VP1-B junction region was retrieved for all included sequences and used for dataset generation.

### Clade assignment

A preliminary Nextclade dataset was generated using unpublished CDC GHOST reference sequences (VP1-2B junction region only) of clades IA, IB, IC, IIA, and IIIA along with the VP1-2B junction region of representative GenBank sequences of clades IIB ([AY644670.1](https://www.ncbi.nlm.nih.gov/nuccore/AY644670.1)) and IIIB ([AB279735.1](https://www.ncbi.nlm.nih.gov/nuccore/AB279735.1)). This preliminary dataset was used for clade assignment of sequences included in this dataset.

### Phylogeny

The phylogeny used in this dataset was generated with IQ-TREE with 1000 bootstraps using the Jukes-Cantor (JC) model. The phylogeny was rooted using the nearest non-human host outgroup with no clade assigned in the previous section, [OR452341.1](https://www.ncbi.nlm.nih.gov/nuccore/OR452341.1) (not included in this dataset). This was the tool, model, and rooting that produced branches with defining mutations for all clades IA, IB, IC, IIA, IIB, IIIA, and IIIB.

### Validation

Unpublished CDC GHOST reference sequences (VP1-2B junction region only) of clades IA, IB, IC, IIA, and IIIA along with the VP1-2B junction region of representative GenBank sequences of clades IIB ([AY644670.1](https://www.ncbi.nlm.nih.gov/nuccore/AY644670.1)) and IIIB ([AB279735.1](https://www.ncbi.nlm.nih.gov/nuccore/AB279735.1)) were analyzed using this dataset. Clade assignment matched that expected.

## What is a Nextclade dataset?

Read more about Nextclade datasets in the [Nextclade documentation](https://docs.nextstrain.org/projects/nextclade/en/stable/user/datasets.html).
