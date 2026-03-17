#!/usr/bin/env python3
"""Generate a demo step2 HTML report without requiring ete3.

Reads demo/clusters/*.fasta and demo/search/*.genes.list,
synthesizes simple gene trees and OG assignments, then writes
demo/step2_report.html using the report_step2.py HTML template.
"""

import json
import html as _html
import random
from collections import defaultdict
from pathlib import Path
import sys

WORKFLOW   = Path(__file__).parent
REPO_ROOT  = WORKFLOW.parent
DEMO_DIR   = REPO_ROOT / "demo"

sys.path.insert(0, str(WORKFLOW))
from report_step2 import (
    HTML_TEMPLATE,
    load_family_info,
    build_family_records,
    build_hg_records,
    get_species_prefix,
)

random.seed(42)

# ── Species tree (hardcoded from demo/demo_tree.nwk, no ete3 needed) ─────────

SP_TREE_DATA = {
  "name": "Root", "dist": 0,
  "children": [
    {"name": "Bilateria", "dist": 1, "children": [
      {"name": "Vertebrata", "dist": 1, "children": [
        {"name": "Mammalia", "dist": 1, "children": [
          {"name": "Mmus", "dist": 1, "leaf": True},
          {"name": "Hsap", "dist": 1, "leaf": True},
        ]},
        {"name": "Actinopterygii", "dist": 1, "children": [
          {"name": "Drer", "dist": 1, "leaf": True},
          {"name": "Xentro", "dist": 1, "leaf": True},
        ]},
      ]},
      {"name": "Ecdysozoa", "dist": 1, "children": [
        {"name": "Diptera", "dist": 1, "children": [
          {"name": "Dmel", "dist": 1, "leaf": True},
          {"name": "Aedaeg", "dist": 1, "leaf": True},
        ]},
        {"name": "Cele", "dist": 1, "leaf": True},
      ]},
    ]},
    {"name": "Cnidaria_like", "dist": 1, "children": [
      {"name": "Nvec", "dist": 1, "leaf": True},
      {"name": "Mlei", "dist": 1, "leaf": True},
    ]},
    {"name": "Invertebrata", "dist": 1, "children": [
      {"name": "Aque", "dist": 1, "leaf": True},
      {"name": "Lophotrochozoa", "dist": 1, "children": [
        {"name": "Spur", "dist": 1, "leaf": True},
        {"name": "Cgig", "dist": 1, "leaf": True},
      ]},
    ]},
  ]
}

def sp_order(node):
    if node.get("leaf"):
        return [node["name"]]
    result = []
    for c in node.get("children", []):
        result.extend(sp_order(c))
    return result

SPECIES_ORDER = sp_order(SP_TREE_DATA)

# ── Clade groupings (manual, mirrors parse_clade_groupings logic) ─────────────

CLADE_DATA = [
    {
        "name": "Root",
        "groups": {sp: grp for grp, sps in [
            ("Bilateria",      ["Mmus","Hsap","Drer","Xentro","Dmel","Aedaeg","Cele"]),
            ("Cnidaria_like",  ["Nvec","Mlei"]),
            ("Invertebrata",   ["Aque","Spur","Cgig"]),
        ] for sp in sps},
    },
    {
        "name": "Bilateria",
        "groups": {sp: grp for grp, sps in [
            ("Vertebrata",  ["Mmus","Hsap","Drer","Xentro"]),
            ("Ecdysozoa",   ["Dmel","Aedaeg","Cele"]),
        ] for sp in sps},
    },
    {
        "name": "Vertebrata",
        "groups": {sp: grp for grp, sps in [
            ("Mammalia",         ["Mmus","Hsap"]),
            ("Actinopterygii",   ["Drer","Xentro"]),
        ] for sp in sps},
    },
]

# ── Gene-tree synthesis from FASTA files ──────────────────────────────────────

def read_fasta_genes(path):
    genes = []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                genes.append(line[1:].strip().split()[0])
    return genes


def _make_subtree(genes, dist=0.1):
    """Return a tree dict for a list of leaf gene IDs."""
    if len(genes) == 1:
        sp = get_species_prefix(genes[0])
        return {"name": genes[0], "dist": dist, "leaf": True, "species": sp}
    mid = len(genes) // 2
    return {
        "name": "",
        "dist": dist,
        "children": [
            _make_subtree(genes[:mid], dist * 0.5),
            _make_subtree(genes[mid:], dist * 0.5),
        ],
    }


def make_gene_tree(genes):
    """Build a simple gene tree: species clades → OG groups → root.

    Returns (tree_dict, ogs) where ogs = {og_name: [gene, ...]}.
    """
    if not genes:
        return {"name": "", "dist": 0, "children": []}, {}

    # Group by species
    by_sp = defaultdict(list)
    for g in genes:
        by_sp[get_species_prefix(g)].append(g)

    # Make per-species subtrees
    sp_subtrees = []
    for sp in SPECIES_ORDER:
        if sp in by_sp:
            sp_subtrees.append(_make_subtree(by_sp[sp], dist=0.05))
    # any species not in SPECIES_ORDER
    for sp, gs in by_sp.items():
        if sp not in SPECIES_ORDER:
            sp_subtrees.append(_make_subtree(gs, dist=0.05))

    if not sp_subtrees:
        return {"name": "", "dist": 0, "children": []}, {}

    # Group species subtrees into 1–3 OGs
    random.shuffle(sp_subtrees)
    n_ogs = min(max(1, len(sp_subtrees) // 3), 3)
    og_groups = [[] for _ in range(n_ogs)]
    for i, st in enumerate(sp_subtrees):
        og_groups[i % n_ogs].append(st)

    ogs = {}
    og_nodes = []
    for i, group in enumerate(og_groups):
        if not group:
            continue
        og_name = f"OG{i+1:03d}"
        og_genes = []
        def collect_leaves(n):
            if n.get("leaf"):
                og_genes.append(n["name"])
            for c in n.get("children", []):
                collect_leaves(c)
        for st in group:
            collect_leaves(st)
        ogs[og_name] = og_genes
        if len(group) == 1:
            node = group[0]
            node["name"] = og_name  # label single-sp subtree as OG
        else:
            node = {
                "name": og_name,
                "dist": 0.3,
                "children": group,
            }
        og_nodes.append(node)

    if len(og_nodes) == 1:
        root = og_nodes[0]
        root["dist"] = 0
    else:
        root = {"name": "", "dist": 0, "children": og_nodes}

    return root, ogs


# ── Build POSSVM-style records ────────────────────────────────────────────────

def build_tree_records(cluster_dir):
    records = []
    all_species = set()
    for fasta in sorted(cluster_dir.glob("*.fasta")):
        stem = fasta.stem          # e.g. adh.Collagen.HG072
        parts = stem.split(".", 2)
        if len(parts) < 3:
            continue
        prefix, family, hg = parts
        genes = read_fasta_genes(fasta)
        if not genes:
            continue
        tree_dict, ogs = make_gene_tree(genes)
        species = sorted({get_species_prefix(g) for g in genes})
        all_species.update(species)
        records.append({
            "id":        stem,
            "hg":        hg,
            "family":    family,
            "prefix":    prefix,
            "n_leaves":  len(genes),
            "species":   species,
            "og_names":  sorted(ogs.keys()),
            "n_ogs":     len(ogs),
            "tree_dict": tree_dict,
            "ogs":       ogs,
        })
    return records, sorted(all_species)


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    cluster_dir = DEMO_DIR / "clusters"
    search_dir  = DEMO_DIR / "search"
    family_info_path = REPO_ROOT / "data" / "gene_families_searchinfo.csv"
    output_path = REPO_ROOT / "docs" / "step2_report.html"

    family_info    = load_family_info(str(family_info_path))
    family_records = build_family_records(search_dir, family_info)
    hg_records     = build_hg_records(cluster_dir, family_info)

    tree_records, all_species = build_tree_records(cluster_dir)
    print(f"Built {len(tree_records)} gene trees, {len(all_species)} species.")

    # Merge species list: tree order + data
    data_sp = set(all_species)
    species_order = [s for s in SPECIES_ORDER if s in data_sp]
    for s in data_sp:
        if s not in species_order:
            species_order.append(s)

    index_records = [{k: v for k, v in r.items() if k not in ("tree_dict", "ogs")}
                     for r in tree_records]

    lazy_parts = []
    for rec in tree_records:
        detail = {"tree": rec["tree_dict"], "ogs": rec["ogs"]}
        tag_id = _html.escape(rec["id"], quote=True)
        lazy_parts.append(
            f'<script type="application/json" id="treedata-{tag_id}">'
            + json.dumps(detail, separators=(",", ":"))
            + "</script>"
        )

    html = (HTML_TEMPLATE
            .replace("%%LAZY_SCRIPTS%%",    "\n".join(lazy_parts))
            .replace("%%SPECIES_ORDER%%",   json.dumps(species_order))
            .replace("%%TREE_DATA%%",       json.dumps(SP_TREE_DATA))
            .replace("%%FAMILY_DATA%%",     json.dumps(family_records))
            .replace("%%HG_DATA%%",         json.dumps(hg_records))
            .replace("%%TREE_INDEX_JSON%%", json.dumps(index_records))
            .replace("%%SPECIES_JSON%%",    json.dumps(all_species))
            .replace("%%CLADE_DATA_JSON%%", json.dumps(CLADE_DATA)))

    output_path.write_text(html, encoding="utf-8")
    print(f"Written: {output_path}")


if __name__ == "__main__":
    main()
