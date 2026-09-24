# TODOs — Phylogenies_v3

Working task list. Added 24/09/2026 from the VGL-channel check (`2022_Mlei/specific_genes/VGIKs/VGL_missing_families.md`,
`2022_Mlei/Labbook.md` 2026-09-24).

## 1. Outlier removal (the biggest TODO)

The pipeline does not remove outlier sequences before or after tree building (no step in `step1.nf`/`step2.nf`,
nothing in `docs/`). Long-branch tips and fragments stay in the HG trees and possvm calls.

- Options to explore (from memory, not yet checked against the current literature):
  - **TreeShrink** (Mai & Mirarab 2018, *BMC Genomics*): per-tree removal of tips that inflate the tree diameter.
  - **Sequence-level coverage filter** before the tree: trimAl `-resoverlap` / `-seqoverlap`.
- This comes before TODO 2: a more permissive search adds more borderline sequences.

## 2. Rescore with HG-derived HMMs (recover divergent members Pfam's GA cutoff drops)

**Measured case**: `Mlei_v05_G018888` is a full-length TRPM (1,626 aa), a tandem duplicate of `G018887`
(84 % identity over the TRPM N-terminal domain, 69 % over the pore region). Its Pfam `Ion_trans` score is
**11.2 vs GA 25** (22.6 with two S1–S4 inserts removed), so the step1 search (`--cut_ga`, `genefam.csv`)
drops it and it never reaches clustering. `G009211` is lost the same way. The likely home is the
ctenophore-only `ion.Ion_trans.HG15`.

Pfam GA is a curator-set round number (`Ion_trans`: GA = TC = 25, NC = 24.9; `Ion_trans_2`: 22.5 / 22.5 / 22.4).

Idea: build a profile HMM from each HG alignment, rescore the proteomes (start with Mlei), and add hits back.

- **A gene can then belong to several HGs and alignments.** This is intended, but the pipeline may assume one
  HG per gene within a family (MCL gives a partition). Check every place that does:
  - per-species annotation tables (`workflow/gather_annotations.py`): one gene → one orthogroup call?
  - possvm / GeneRax / step4 presence–absence: a gene counted in two HGs double-counts or collides.
  - **Cross-family** multi-membership already exists and works: e.g. `Mlei_v05_G014137` is in both
    `ion.Ion_trans.HG14` and `ion.Ion_trans_2.HG4`. The new case is **within-family** multi-HG membership.
- Threshold: home-built HMMs have no GA. Calibrate per HG (e.g. above the weakest true member's score, and a
  hit's best HG must beat the others), and cross-check with E-value.
- Circularity: build the HMMs with the query species left out, or it mostly recovers itself.
- A new hit only gets a family; putting it on a tree still needs a separate step.
- Cheap first look: rescore Mlei proteins with Pfam `Ion_trans`/`Ion_trans_2` without the GA cutoff and count
  hits in the 10–25 score band.

## 3. Allow custom HMMs in step1 (simpler; also a prerequisite for TODO 2)

The code is mostly there already (read 24/09/2026):

- `phylogeny/helper/hmmsearch.py` `do_hmmsearch()` looks for `<hmm_dir>/<name>.hmm` first and fetches from
  `--pfam_db` only if that file is missing or empty.
- The threshold column of `genefam.csv` is `GA` (→ `--cut_ga`) or a number (→ `--domE <value>`).
- `phylogeny/main.py` takes `--hmm_dir` (default `hmms/`).

**FIXED in `step1.nf` (24/09/2026, uncommitted).** Before, `SEARCH` never passed `--hmm_dir`, so every model
was fetched from Pfam.
- Now `params.hmm_dir` (default `null`, in `nextflow.config`) is staged into `SEARCH` as a `path` input
  (`custom_hmms/`, so an edited `.hmm` invalidates the cache) and passed as `--hmm_dir custom_hmms`.
  `[]` is used when it is unset, which keeps the old behaviour.
- Tested locally on 3 Mlei TRPMs (`Ion_trans`, `phylo` env):
  - no `hmm_dir`: Pfam is fetched; genes.list = G004882, G018887.
  - custom `Ion_trans.hmm` with GA 40: logs `Found custom_hmms/Ion_trans.hmm`; genes.list = G004882 only
    (G018887 scores 34.7), so the custom model changes the output.
- **Still open:**
  - `step1.smk` has the same gap. Its config is the template
    `2022_Mlei/Phylogenies/workflow/templates/step1.yaml`.
  - A custom HMM usually has no `GA` line, so `--cut_ga` fails on it. Either give that family a numeric
    threshold in `genefam.csv` or write GA/TC/NC into the `.hmm` header.
  - ⚠ **Silent shadowing**: a custom file named like a Pfam model replaces Pfam's. The only trace is the
    log line `Found custom_hmms/<name>.hmm` vs `fetch`.

⚠ **Found while testing, relevant to TODO 2:** lowering the threshold alone does NOT recover G018888.
Even with a custom GA of 10, `hmmsearch` never reports it, because its acceleration filters
(MSV/Viterbi/Forward) discard it before any threshold applies. It only appears with `--max`
(score 12.7 seq / 11.2 dom). A recovery pass needs `--max` or a better model (HG-derived), not just a lower
cutoff. `do_hmmsearch()` has no `--max` option.
