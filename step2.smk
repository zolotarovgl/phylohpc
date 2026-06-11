# step2.smk — Snakemake port of step2.nf
#
# Quick start:
#   snakemake -s step2.smk -j4                          # local, 4 cores
#   snakemake -s step2.smk -j4 --config run_generax=1   # with GeneRax
#   snakemake -s step2.smk -n                           # dry-run
#   snakemake -s step2.smk --dag | dot -Tsvg > dag.svg  # visualise DAG
#
# Override any config key on the command line:
#   snakemake -s step2.smk -j4 --config tree_method=iqtree2 outdir=results2

import os

configfile: "config/step2.yaml"

OUTDIR        = config["outdir"]
REFSPECIES    = config["refspecies"]
REFNAMES      = config["refnames"]
SPECIES_TREE  = config["species_tree"]
TREE_METHOD   = config["tree_method"]
IQTREE2_MODEL = config["iqtree2_model"]
MAFFT_OPT     = config["mafft_opt"]
RUN_GENERAX   = bool(int(config.get("run_generax", 0)))
NCPU_GENERAX  = int(config["ncpu_generax"])
MAX_SPR       = int(config["max_spr"])
SUBS_MODEL    = config["subs_model"]

PHYLO = os.path.join(workflow.basedir, "phylogeny", "main.py")
WF    = os.path.join(workflow.basedir, "workflow")

with open(config["ids"]) as _f:
    IDS = [l.strip() for l in _f if l.strip() and not l.startswith("#")]

# POSSVM inserts the tree filename stem into its output names:
#   *.treefile       → treefile
#   *.generax.tree   → generax  (strips last .tree extension)
PVM_INFIX = "generax" if RUN_GENERAX else "treefile"


# ── target ─────────────────────────────────────────────────────────────────────

rule all:
    input:
        f"{OUTDIR}/report_step2.html",


# ── align ──────────────────────────────────────────────────────────────────────

rule aln:
    input:
        fasta = f"{OUTDIR}/clusters/{{id}}.fasta",
    output:
        aln = f"{OUTDIR}/align/{{id}}.aln.fasta",
    threads: 4
    shell:
        """
        python {PHYLO} align -f {{input.fasta}} -o {{output.aln}} -c {{threads}} -m "{MAFFT_OPT}"
        python {WF}/remove_gaponly.py {{output.aln}} {{output.aln}}.tmp
        mv {{output.aln}}.tmp {{output.aln}}
        """


# ── gene tree ──────────────────────────────────────────────────────────────────

rule phy:
    input:
        aln = f"{OUTDIR}/align/{{id}}.aln.fasta",
    output:
        tree = f"{OUTDIR}/gene_trees/{{id}}.treefile",
    threads: 4
    shell:
        """
        python {PHYLO} phylogeny \
            -f {{input.aln}} \
            --outprefix {OUTDIR}/gene_trees/{{wildcards.id}} \
            -c {{threads}} \
            --method {TREE_METHOD} \
            --iqtree2_model {IQTREE2_MODEL}
        """


# ── POSSVM on original trees ───────────────────────────────────────────────────
# Always runs; when run_generax=1 the report uses these for the toggle.

rule pvm_prev:
    input:
        tree     = f"{OUTDIR}/gene_trees/{{id}}.treefile",
        refnames = REFNAMES,
    output:
        nwk   = f"{OUTDIR}/possvm_prev/{{id}}.treefile.ortholog_groups.newick",
        csv   = f"{OUTDIR}/possvm_prev/{{id}}.treefile.ortholog_groups.csv",
        pairs = f"{OUTDIR}/possvm_prev/{{id}}.treefile.pairs_orthologs.csv",
    shell:
        """
        python {PHYLO} possvm \
            -t {{input.tree}} \
            --refsps {REFSPECIES} \
            -r {{input.refnames}} \
            -o {OUTDIR}/possvm_prev/{{wildcards.id}}.
        """


# ── GeneRax (only when run_generax=1) ─────────────────────────────────────────

if RUN_GENERAX:
    rule generax:
        input:
            aln          = f"{OUTDIR}/align/{{id}}.aln.fasta",
            tree         = f"{OUTDIR}/gene_trees/{{id}}.treefile",
            species_tree = SPECIES_TREE,
        output:
            tree = f"{OUTDIR}/generax/{{id}}.generax.tree",
            log  = f"{OUTDIR}/generax/{{id}}.generax.log",
        threads: NCPU_GENERAX
        shell:
            """
            export OMP_NUM_THREADS={{threads}}
            export OPENBLAS_NUM_THREADS={{threads}}
            python {PHYLO} generax \
                --name {{wildcards.id}} \
                --alignment {{input.aln}} \
                --gene_tree {{input.tree}} \
                --species_tree {{input.species_tree}} \
                --output_dir {OUTDIR}/generax/{{wildcards.id}}_generax \
                --subs_model {SUBS_MODEL} \
                --max_spr {MAX_SPR} \
                --cpus {{threads}} \
                --logfile {{output.log}} \
                --outfile {{output.tree}}
            """


# ── POSSVM on final trees ──────────────────────────────────────────────────────

rule pvm:
    input:
        tree     = (f"{OUTDIR}/generax/{{id}}.generax.tree"
                    if RUN_GENERAX else
                    f"{OUTDIR}/gene_trees/{{id}}.treefile"),
        refnames = REFNAMES,
    output:
        nwk   = f"{OUTDIR}/possvm/{{id}}.{PVM_INFIX}.ortholog_groups.newick",
        csv   = f"{OUTDIR}/possvm/{{id}}.{PVM_INFIX}.ortholog_groups.csv",
        pairs = f"{OUTDIR}/possvm/{{id}}.{PVM_INFIX}.pairs_orthologs.csv",
    shell:
        """
        python {PHYLO} possvm \
            -t {{input.tree}} \
            --refsps {REFSPECIES} \
            -r {{input.refnames}} \
            -o {OUTDIR}/possvm/{{wildcards.id}}.
        """


# ── HTML report ────────────────────────────────────────────────────────────────

def _report_inputs(wildcards):
    pvm = expand(f"{OUTDIR}/possvm/{{id}}.{PVM_INFIX}.ortholog_groups.newick", id=IDS)
    pvm_prev = (expand(f"{OUTDIR}/possvm_prev/{{id}}.treefile.ortholog_groups.newick", id=IDS)
                if RUN_GENERAX else [])
    return pvm + pvm_prev


rule report:
    input:
        _report_inputs,
    output:
        f"{OUTDIR}/report_step2.html",
    params:
        prev_flag = (f"--possvm_prev_dir {OUTDIR}/possvm_prev" if RUN_GENERAX else ""),
    shell:
        """
        python {WF}/report_step2.py \
            --possvm_dir  {OUTDIR}/possvm \
            {{params.prev_flag}} \
            --search_dir  {OUTDIR}/search \
            --cluster_dir {OUTDIR}/clusters \
            --family_info data/gene_families_searchinfo.csv \
            --species_tree data/species_tree.full.newick \
            --output {{output}}
        """
