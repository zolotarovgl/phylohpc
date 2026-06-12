# archive/

Dead code kept for reference (recoverable from git history). Not used by any
live pipeline.

- `old-bash-pipeline/` — the original shell-script pipeline (`s01_search.sh`,
  `run_phylo.sh`, …) superseded by the Nextflow/Snakemake steps, plus its
  `config.txt`.
- `generax.nf` — standalone GeneRax-only Nextflow entrypoint; superseded by the
  GeneRax path inside `step2.nf`. Was never wired into the documented workflow.
- `predict_resources.R` — superseded by `workflow/predict_resources.py`.
