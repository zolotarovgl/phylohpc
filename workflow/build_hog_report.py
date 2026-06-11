#!/usr/bin/env python3
"""
Build the interactive hOG hierarchy Sankey HTML directly from Snakemake
ancestry pipeline outputs, without requiring separate link/stats TSV files.

Logic (mirrors the R explore_hOGs.R script and link_hog_levels.py):
  For each HG and each pair of consecutive nodes, count shared genes between
  every pair of OGs (one from the broader level, one from the narrower level).
  Non-zero intersections become Sankey edges; per-OG stats are computed from
  the filtered (in-clade) gene sets.

Expected directory layout produced by step4_ancestry.smk:
  {ancestry_dir}/{node}/possvm/{hg}.ortholog_groups.csv
  {ancestry_dir}/{node}/{node}.in_species.txt
  {ancestry_dir}/{node}/{node}.pruned.tree

Usage:
  python workflow/build_hog_report.py \\
      --ancestry_dir results/ancestry \\
      --nodes        "Metazoa,Euarchontoglires" \\
      --ids          ancestry_ids.txt \\
      --output       results/ancestry/hog_hierarchy.html
"""

import argparse
import json
import sys
from collections import defaultdict
from pathlib import Path

import pandas as pd

# Re-use the HTML template, data builder, and Newick parser from the companion
# visualize script so any UI improvements there are automatically inherited.
sys.path.insert(0, str(Path(__file__).parent))
from visualize_hog_hierarchy import HTML_TEMPLATE, build_data, parse_newick_tree


# ── Helpers ───────────────────────────────────────────────────────────────────

def species_prefix(gene: str) -> str:
    """Return the species code: first token before '_' or '.'."""
    for sep in ("_", "."):
        if sep in gene:
            return gene.split(sep)[0]
    return gene


def _parse_float_field(val):
    """Parse a simple float field; returns None for missing/NA values."""
    if val is None:
        return None
    s = str(val).strip()
    if s in ("", "NA", "nan", "None"):
        return None
    try:
        return float(s)
    except ValueError:
        return None


def _parse_support_field(val):
    """Parse a support value that may be slash-separated (e.g. '87.0/87.0/87.0').
    Returns the minimum across all listed values (most conservative)."""
    if val is None:
        return None
    s = str(val).strip()
    if s in ("", "NA", "nan", "None"):
        return None
    try:
        parts = [float(x) for x in s.split("/") if x.strip()]
        return min(parts) if parts else None
    except ValueError:
        return None


def _parse_ref_ortholog(val):
    """Parse reference_ortholog; may be slash-separated when ambiguous."""
    if val is None:
        return None
    s = str(val).strip()
    return None if s in ("", "NA", "nan", "None") else s


def load_possvm_csv(path: Path) -> dict[str, dict]:
    """
    Parse a POSSVM ortholog_groups.csv → {gene: {og, og_support, ref_ortholog, ref_support}}.
    Columns: gene  orthogroup  orthogroup_support  reference_ortholog  reference_support
    Note: reference_ortholog and reference_support may be slash-separated when a gene
    maps equally well to multiple reference orthologs.
    """
    if not path.exists() or path.stat().st_size == 0:
        return {}
    try:
        # Read all as str to avoid pandas misinterpreting slash-separated values
        df = pd.read_csv(path, sep="\t", header=0, dtype=str)
        df = df.dropna(subset=["gene", "orthogroup"])
        result = {}
        for row in df.itertuples(index=False):
            result[str(row.gene)] = {
                "og":          str(row.orthogroup),
                "og_support":  _parse_float_field(getattr(row, "orthogroup_support", None)),
                "ref_ortholog": _parse_ref_ortholog(getattr(row, "reference_ortholog", None)),
                "ref_support": _parse_support_field(getattr(row, "reference_support", None)),
            }
        return result
    except Exception as exc:
        print(f"  WARNING: cannot read {path}: {exc}", file=sys.stderr)
        return {}


def load_in_species(path: Path) -> set[str]:
    if not path.exists():
        return set()
    return {ln.strip() for ln in path.read_text().splitlines() if ln.strip()}


# ── Per-HG processing ─────────────────────────────────────────────────────────

def process_hg(
    hg: str,
    nodes: list[str],
    ancestry_dir: Path,
) -> tuple[list[dict], list[dict]]:
    """
    Load POSSVM outputs for one HG, filter to in-clade genes, then compute:
      - stat_rows : per-OG gene/species counts + support values → og_stats columns
      - link_rows : cross-level gene-intersection edges → og_links columns

    Returns (stat_rows, link_rows).
    """
    # ── Load and filter per-node gene → record maps ───────────────────────────
    # level_data[node] = {gene: {og, og_support, ref_ortholog, ref_support}}
    level_data = {}  # type: dict
    level_total_species = {}  # type: dict
    for node in nodes:
        csv_path   = ancestry_dir / node / "possvm" / f"{hg}.ortholog_groups.csv"
        in_sp_path = ancestry_dir / node / f"{node}.in_species.txt"
        raw = load_possvm_csv(csv_path)
        if not raw:
            print(f"  {hg}/{node}: no POSSVM output found", file=sys.stderr)
            continue
        in_sps = load_in_species(in_sp_path)
        level_total_species[node] = len(in_sps) if in_sps else None
        if in_sps:
            filtered = {g: rec for g, rec in raw.items()
                        if species_prefix(g) in in_sps}
            dropped = len(raw) - len(filtered)
            if dropped:
                print(f"  {hg}/{node}: dropped {dropped} out-of-clade genes "
                      f"({len(filtered)} kept)", file=sys.stderr)
        else:
            filtered = raw
        if filtered:
            level_data[node] = filtered

    # ── Per-OG statistics ─────────────────────────────────────────────────────
    stat_rows: list[dict] = []
    for node in nodes:
        records = level_data.get(node)
        if not records:
            continue
        og_records: dict[str, list[dict]] = defaultdict(list)
        for gene, rec in records.items():
            og_records[rec["og"]].append({"gene": gene, **rec})
        for og, gene_recs in sorted(og_records.items()):
            genes = [r["gene"] for r in gene_recs]
            sps   = sorted({species_prefix(g) for g in genes})
            # og_support is constant per OG — take first non-None value
            og_sup = next((r["og_support"] for r in gene_recs if r["og_support"] is not None), None)
            # ref_ortholog is constant per OG — take first non-None value
            ref_ortholog = next((r["ref_ortholog"] for r in gene_recs if r["ref_ortholog"] is not None), None)
            # ref_support: minimum across genes (most conservative)
            ref_sups = [r["ref_support"] for r in gene_recs if r["ref_support"] is not None]
            ref_sup = min(ref_sups) if ref_sups else None
            # ref_genes: {species: "name1/name2"} for genes with a named reference ortholog
            ref_map: dict[str, set] = {}
            for r in gene_recs:
                if r["ref_ortholog"] is not None:
                    sp = species_prefix(r["gene"])
                    ref_map.setdefault(sp, set()).add(r["ref_ortholog"])
            ref_genes = {sp: "/".join(sorted(names)) for sp, names in ref_map.items()}
            stat_rows.append({
                "hg":          hg,
                "level":       node,
                "og":          og,
                "n_genes":     len(genes),
                "n_species":   len(sps),
                "n_total_species": level_total_species.get(node) or len(sps),
                "species":     ",".join(sps),
                "genes":       genes,
                "og_support":  og_sup,
                "ref_ortholog": ref_ortholog,
                "ref_support": ref_sup,
                "ref_genes":   ref_genes,
            })

    # ── Cross-level links (gene intersection between consecutive nodes) ────────
    link_rows: list[dict] = []
    available = [n for n in nodes if n in level_data]
    for parent_lv, child_lv in zip(available[:-1], available[1:]):
        parent_map = {g: rec["og"] for g, rec in level_data[parent_lv].items()}
        child_map  = {g: rec["og"] for g, rec in level_data[child_lv].items()}
        shared = set(parent_map) & set(child_map)
        if not shared:
            continue
        edge_count: dict[tuple[str, str], int] = defaultdict(int)
        for gene in shared:
            edge_count[(parent_map[gene], child_map[gene])] += 1
        for (p_og, c_og), n in sorted(edge_count.items()):
            link_rows.append({
                "hg":        hg,
                "parent_og": p_og,
                "child_og":  c_og,
                "n_genes":   n,
            })

    return stat_rows, link_rows


# ── Main ──────────────────────────────────────────────────────────────────────

def main() -> None:
    parser = argparse.ArgumentParser(
        description="Build hOG hierarchy Sankey HTML from Snakemake ancestry outputs."
    )
    parser.add_argument("--ancestry_dir", required=True,
                        help="Root ancestry directory (per-node subdirs inside)")
    parser.add_argument("--nodes", required=True,
                        help="Clade names, comma-separated, broad→narrow")
    parser.add_argument("--ids", required=True,
                        help="File with HG IDs, one per line")
    parser.add_argument("--output", required=True, help="Output HTML path")
    args = parser.parse_args()

    ancestry_dir = Path(args.ancestry_dir)
    nodes  = [n.strip() for n in args.nodes.split(",") if n.strip()]
    hg_ids = [ln.strip()
              for ln in Path(args.ids).read_text().splitlines() if ln.strip()]

    print(f"Nodes : {nodes}",  file=sys.stderr)
    print(f"HGs   : {hg_ids}", file=sys.stderr)

    # ── Collect stats and links across all HGs ────────────────────────────────
    all_stats: list[dict] = []
    all_links: list[dict] = []
    for hg in hg_ids:
        stat_rows, link_rows = process_hg(hg, nodes, ancestry_dir)
        all_stats.extend(stat_rows)
        all_links.extend(link_rows)
        print(f"  {hg}: {len(stat_rows)} OGs, {len(link_rows)} edges",
              file=sys.stderr)

    stats_df = pd.DataFrame(all_stats) if all_stats else pd.DataFrame()
    links_df = pd.DataFrame(all_links) if all_links else pd.DataFrame()

    # ── Load pruned species trees ─────────────────────────────────────────────
    trees: dict[str, dict] = {}
    for node in nodes:
        tree_path = ancestry_dir / node / f"{node}.pruned.tree"
        tree = parse_newick_tree(tree_path)
        if tree:
            trees[node] = tree
        else:
            print(f"  WARNING: pruned tree not found for '{node}'", file=sys.stderr)

    # ── Load POSSVM-annotated gene trees (broadest available level per HG) ────
    gene_trees: dict[str, str] = {}
    for hg in hg_ids:
        for node in nodes:   # nodes ordered broad → narrow
            nwk = ancestry_dir / node / "possvm" / f"{hg}.ortholog_groups.newick"
            if nwk.exists() and nwk.stat().st_size > 0:
                text = nwk.read_text(encoding="utf-8").strip()
                if text:
                    gene_trees[hg] = text
                    break
    print(f"Gene trees loaded: {len(gene_trees)}/{len(hg_ids)}", file=sys.stderr)

    # ── Build data JSON and render HTML ───────────────────────────────────────
    data      = build_data(stats_df, links_df, nodes, trees, gene_trees=gene_trees)
    data_json = json.dumps(data, separators=(",", ":"))
    html      = HTML_TEMPLATE.replace("%%DATA_JSON%%", data_json)
    Path(args.output).write_text(html, encoding="utf-8")
    print(f"Written: {args.output}", file=sys.stderr)


if __name__ == "__main__":
    main()
