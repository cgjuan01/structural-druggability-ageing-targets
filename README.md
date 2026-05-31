# Structure-based druggability assessment of GNN-prioritised ageing targets

**Structure-based druggability screening** of computationally-prioritised drug
targets, combining **graph-neural-network target prioritisation**, **protein
structure** (experimental PDB + **AlphaFold**), and **fpocket** binding-pocket
detection. The pipeline takes targets ranked by a graph attention network and
asks the next drug discovery question: *are they structurally druggable?*

Keywords: structure-based drug discovery · druggability · AlphaFold · fpocket ·
binding-pocket detection · graph neural networks · target prioritisation ·
cheminformatics.

## Motivation

A graph attention network (GAT) prioritised exercise-responsive,
ageing-clock-supported genes as candidate targets (see *Context* below). Target
prioritisation answers *which proteins matter*; it does not answer *which can be
drugged*. This project assesses the **structural druggability** of the top-ranked
targets — i.e. whether each protein presents a ligand-binding pocket with
drug-like geometry — using structure-based pocket detection independent of prior
medicinal chemistry attention.

## Pipeline

1. `src/01_resolve_structures.py` — resolve gene symbols to UniProt accessions and
   determine structure availability (experimental PDB vs AlphaFold-only).
2. `src/02_fpocket_druggability.py` — download a representative structure per
   target and score binding pockets with **fpocket**, retaining the top
   Druggability Score (0-1; >= 0.5 conventionally druggable). Resumable.
3. `src/03_analyse.py` — summarise the druggability landscape, test druggability
   vs GNN rank, list top candidates, and render the figure.

## Headline results

Of the top-1000 GNN-ranked ageing-clock-supported genes, **698 had experimental
crystal structures**; **679** were successfully scored by fpocket.

- **~58 % of the scored (crystal-structure) targets are structurally druggable**
  (fpocket Druggability Score >= 0.5); 240 are highly druggable (>= 0.7).
- **Druggability is independent of GNN rank** — the highest-ranked
  exercise/ageing genes are no more druggable than lower-ranked ones (median
  ~0.54 vs ~0.57). Biological importance and chemical tractability are
  orthogonal here.
- A shortlist of **high-pocket, high-rank candidates** (e.g. RUNX1, IKZF1, FGF21,
  RAB27A) is reported as hypotheses for tractable targets in the exercise-ageing
  network.

![druggability landscape](results/druggability_landscape.png)

## Honest scope and limitations

This is a **hypothesis-generating structural screen**, not validated
druggability. Specifically:

- **One representative structure per target.** A given crystal form may lack the
  biologically relevant pocket (e.g. a catalytic site absent from a particular
  construct), so *individual* gene scores are hypotheses; the *aggregate*
  distribution is robust, individual calls are not. Some scores are clearly
  structure artifacts (e.g. anomalously large pocket volumes from multi-domain
  complexes). Rigorous per-target analysis would curate multiple structures.
- **Crystal-structure subset.** The 58 % applies to targets with experimental
  structures, which are disproportionately well-studied proteins; the figure is
  not claimed for the whole network. ~300 AlphaFold-only genes are handled
  separately as a lower-confidence set (predicted pockets can be model artifacts).
- **fpocket Druggability Score is a validated heuristic, not ground truth.** A
  druggable pocket means drug-like geometry exists, not that a useful drug can be
  made or will bind.

## Context

The prioritised targets derive from a graph attention network analysis of
exercise-responsive genes with biological-ageing-clock support (Horvath, Hannum,
DunedinPACE, PhenoAge, Proteomic clock). See:

> Juan CG, Ntasis L (2026). Multi-omic deep learning identifies
> exercise-responsive ageing pathways in humans. medRxiv [preprint].

This druggability assessment is an independent extension of that work.

## Run

```bash
pip install matplotlib
# fpocket required (conda: mamba install -c conda-forge -c bioconda fpocket)
python src/01_resolve_structures.py data/top1000_clock.txt resolved.json
python src/02_fpocket_druggability.py resolved.json fpocket_results.jsonl
python src/03_analyse.py fpocket_results.jsonl data/top1000_clock.txt
```
