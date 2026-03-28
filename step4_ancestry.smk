"""
Snakemake re-implementation of step4.ancestry.nf

Usage:
    snakemake -s step4_ancestry.smk --cores 24

Override any config value on the command line:
    snakemake -s step4_ancestry.smk --cores 24 \
        --config node_names="Metazoa,Bilateria,Euarchontoglires" ids=ancestry_ids.txt
"""

import os
from pathlib import Path

configfile: "config_ancestry.yaml"

# ── Parameters ────────────────────────────────────────────────────────────────
NODE_NAMES  = config.get("node_names",      "Metazoa,Bilateria,Euarchontoglires")
IDS_FILE    = config.get("ids",             "ancestry_ids.txt")
OUTDIR      = config.get("outdir",          "results")
SPEC_TREE   = config.get("species_tree",    "data/species_tree.full.newick")
REFSPECIES  = config.get("refspecies",      "Mmus")
REFNAMES    = config.get("refnames",        "data/Mmus_gene_names.csv")
GT_DIR      = config.get("gene_trees_dir",  os.path.join(OUTDIR, "gene_trees"))
TREE_SUFF = "generax.tree"


PDIR  = os.path.dirname(os.path.abspath(workflow.snakefile))
NODES = [n.strip() for n in NODE_NAMES.split(",") if n.strip()]

with open(IDS_FILE) as fh:
    ALL_HGS = [l.strip() for l in fh if l.strip()]

# Only keep HGs whose gene tree exists (mirrors NF .filter { treefile.exists() })
print(GT_DIR)
print(TREE_SUFF)
HG_IDS = [hg for hg in ALL_HGS if Path(f"{GT_DIR}/{hg}.{TREE_SUFF}").exists()]

if not HG_IDS:
    raise ValueError(f"No gene trees found in {GT_DIR}. Check 'gene_trees_dir' config.")


# ── Target ────────────────────────────────────────────────────────────────────
rule all:
    input:
        f"{OUTDIR}/ancestry/hog_hierarchy.html"


# ── Rule 1: extract per-clade species lists and pruned subtree ────────────────
rule extract_clade:
    input:
        tree = SPEC_TREE
    output:
        in_sp  = f"{OUTDIR}/ancestry/{{node}}/{{node}}.in_species.txt",
        ign_sp = f"{OUTDIR}/ancestry/{{node}}/{{node}}.ignore_species.txt",
        ptree  = f"{OUTDIR}/ancestry/{{node}}/{{node}}.pruned.tree"
    params:
        script    = os.path.join(PDIR, "workflow/extract_clade.py"),
        outprefix = lambda wc, output: str(Path(output.in_sp).parent / wc.node)
    shell:
        """
        python {params.script} \
            --tree       {input.tree} \
            --node       "{wildcards.node}" \
            --out_prefix {params.outprefix}
        """


# ── Rule 2: POSSVM per clade × HG ─────────────────────────────────────────────
rule pvm_clade:
    input:
        treefile = f"{GT_DIR}/{{hg}}.{TREE_SUFF}",
        ign_sp   = f"{OUTDIR}/ancestry/{{node}}/{{node}}.ignore_species.txt",
        refnames = REFNAMES
    output:
        csv = f"{OUTDIR}/ancestry/{{node}}/possvm/{{hg}}.ortholog_groups.csv"
    params:
        script     = os.path.join(PDIR, "phylogeny/main.py"),
        refspecies = REFSPECIES,
        outdir     = lambda wc, output: str(Path(output.csv).parent)
    retries: 3
    shell:
        """
        python phylogeny/submodules/possvm-orthology/possvm.py \
            -i         {input.treefile} \
            -min_support_transfer 0\
            --refsps   {params.refspecies} \
            -r         {input.refnames} \
            --outgroup {input.ign_sp} \
            -method    lpa \
            -o         {params.outdir} \
            -p         {wildcards.hg} \
            -ogprefix  {wildcards.node}.{wildcards.hg}.
        """


# ── Rule 3: link OGs across clade levels for one HG ──────────────────────────
rule link_hogs:
    input:
        csvs = lambda wc: expand(
            f"{OUTDIR}/ancestry/{{node}}/possvm/{{node}}.{wc.hg}.ortholog_groups.csv",
            node=NODES
        ),
        in_species = expand(
            f"{OUTDIR}/ancestry/{{node}}/{{node}}.in_species.txt",
            node=NODES
        )
    output:
        links = f"{OUTDIR}/ancestry/hog_links/{{hg}}.og_links.tsv",
        stats = f"{OUTDIR}/ancestry/hog_links/{{hg}}.og_stats.tsv"
    params:
        script = os.path.join(PDIR, "workflow/link_hog_levels.py"),
        levels = NODE_NAMES
    shell:
        """
        python {params.script} \
            --hg           '{wildcards.hg}' \
            --levels       '{params.levels}' \
            --csvs         {input.csvs} \
            --in_species   {input.in_species} \
            --output_links {output.links} \
            --output_stats {output.stats}
        """


# ── Rule 4: build Sankey HTML directly from POSSVM outputs ───────────────────
# Uses build_hog_report.py which inlines the link/stats computation,
# so no separate link_hogs step is needed.
rule build_report:
    input:
        csvs  = expand(
            f"{OUTDIR}/ancestry/{{node}}/possvm/{{hg}}.ortholog_groups.csv",
            node=NODES, hg=HG_IDS
        ),
        trees = expand(f"{OUTDIR}/ancestry/{{node}}/{{node}}.pruned.tree", node=NODES),
        report_script = os.path.join(PDIR, "workflow/build_hog_report.py"),
        viz_script    = os.path.join(PDIR, "workflow/visualize_hog_hierarchy.py")
    output:
        html = f"{OUTDIR}/ancestry/hog_hierarchy.html"
    params:
        ancestry_dir = f"{OUTDIR}/ancestry",
        nodes        = NODE_NAMES,
        ids          = IDS_FILE
    shell:
        """
        python {input.report_script} \
            --ancestry_dir {params.ancestry_dir} \
            --nodes        '{params.nodes}' \
            --ids          {params.ids} \
            --output       {output.html}
        """


# ── Rule 4b (kept for reference): old two-step link+visualize ─────────────────
rule visualize_hierarchy:
    input:
        links = expand(f"{OUTDIR}/ancestry/hog_links/{{hg}}.og_links.tsv", hg=HG_IDS),
        stats = expand(f"{OUTDIR}/ancestry/hog_links/{{hg}}.og_stats.tsv", hg=HG_IDS),
        trees = expand(f"{OUTDIR}/ancestry/{{node}}/{{node}}.pruned.tree",  node=NODES),
        viz_script = os.path.join(PDIR, "workflow/visualize_hog_hierarchy.py")
    output:
        html = f"{OUTDIR}/ancestry/hog_hierarchy_linked.html"
    params:
        levels = NODE_NAMES
    shell:
        """
        python {input.viz_script} \
            --links  {input.links} \
            --stats  {input.stats} \
            --trees  {input.trees} \
            --levels '{params.levels}' \
            --output {output.html}
        """
