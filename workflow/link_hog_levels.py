#!/usr/bin/env python3
"""
Link orthogroups across hierarchical taxonomic levels for a single HG.

For each pair of adjacent taxonomic levels (ordered broad→narrow via --levels),
shared genes between the POSSVM outputs at adjacent levels determine the
parent→child OG relationships.

Input CSV filename convention: {level}.{hg}.ortholog_groups.csv
POSSVM CSV columns (no header): gene_id  og_name  ref_og  ref_gene_name

Outputs:
  {hg}.og_links.tsv  — parent_level, parent_og, child_level, child_og, n_genes
  {hg}.og_stats.tsv  — level, og, n_genes, n_species, species
"""

import argparse
import sys
from collections import defaultdict
from pathlib import Path

import pandas as pd


# ── Helpers ────────────────────────────────────────────────────────────────────

def get_species_prefix(gene_id: str) -> str:
    for sep in ("_", "."):
        if sep in gene_id:
            return gene_id.split(sep)[0]
    return gene_id


def parse_possvm_csv(path: Path) -> dict[str, str]:
    """Return {gene_id: og_name} from a POSSVM ortholog_groups.csv file."""
    mapping: dict[str, str] = {}
    if not path.exists() or path.stat().st_size == 0:
        return mapping
    with path.open() as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) >= 2:
                mapping[parts[0]] = parts[1]
    return mapping


def extract_level(csv_path: str, hg: str) -> str:
    """
    Extract the clade level name from a POSSVM CSV filename.
    Filename format: {level}.{hg}.ortholog_groups.csv
    """
    name = Path(csv_path).name
    suffix = f".{hg}.ortholog_groups.csv"
    if name.endswith(suffix):
        return name[: -len(suffix)]
    # Fallback: strip known extension, then strip hg suffix from right
    stem = name.replace(".ortholog_groups.csv", "").replace(".ortholog_groups", "")
    hg_suffix = f".{hg}"
    if stem.endswith(hg_suffix):
        return stem[: -len(hg_suffix)]
    return stem


# ── Main ───────────────────────────────────────────────────────────────────────

def main() -> None:
    parser = argparse.ArgumentParser(
        description="Link OGs across hierarchical clade levels for one HG."
    )
    parser.add_argument("--hg",           required=True,
                        help="HG identifier, e.g. tfs.RFX.HG1")
    parser.add_argument("--levels",       required=True,
                        help="Ordered clade names, comma-separated (broad→narrow)")
    parser.add_argument("--csvs",         nargs="+", required=True,
                        help="POSSVM ortholog_groups CSV files for this HG")
    parser.add_argument("--output_links", required=True,
                        help="Output edges TSV: parent_og → child_og")
    parser.add_argument("--output_stats", required=True,
                        help="Output stats TSV: per-OG gene/species counts")
    args = parser.parse_args()

    ordered_levels = [lv.strip() for lv in args.levels.split(",") if lv.strip()]

    # ── Load each CSV and map it to its clade level ────────────────────────────
    level_data: dict[str, dict[str, str]] = {}   # level → {gene: og}
    for csv_path in args.csvs:
        level = extract_level(csv_path, args.hg)
        if level not in ordered_levels:
            print(
                f"  WARNING: Level '{level}' inferred from '{Path(csv_path).name}' "
                f"not in --levels ({ordered_levels}); skipping.",
                file=sys.stderr,
            )
            continue
        gene_og = parse_possvm_csv(Path(csv_path))
        if gene_og:
            level_data[level] = gene_og
        else:
            print(
                f"  WARNING: Empty/missing CSV for level '{level}': {csv_path}",
                file=sys.stderr,
            )

    # ── Per-OG statistics ──────────────────────────────────────────────────────
    stat_rows: list[dict] = []
    for level in ordered_levels:
        gene_og = level_data.get(level)
        if not gene_og:
            continue
        og_genes: dict[str, list[str]] = defaultdict(list)
        for gene, og in gene_og.items():
            og_genes[og].append(gene)
        for og, genes in sorted(og_genes.items()):
            species = sorted({get_species_prefix(g) for g in genes})
            stat_rows.append({
                "hg":       args.hg,
                "level":    level,
                "og":       og,
                "n_genes":  len(genes),
                "n_species": len(species),
                "species":  ",".join(species),
            })

    pd.DataFrame(
        stat_rows,
        columns=["hg", "level", "og", "n_genes", "n_species", "species"],
    ).to_csv(args.output_stats, sep="\t", index=False)

    # ── Cross-level OG edges ───────────────────────────────────────────────────
    # Only iterate over levels that actually have data, in declared order.
    available = [lv for lv in ordered_levels if lv in level_data]
    link_rows: list[dict] = []

    for parent_lv, child_lv in zip(available[:-1], available[1:]):
        parent_map = level_data[parent_lv]
        child_map  = level_data[child_lv]

        shared = set(parent_map) & set(child_map)
        if not shared:
            continue

        edge_count: dict[tuple[str, str], int] = defaultdict(int)
        for gene in shared:
            edge_count[(parent_map[gene], child_map[gene])] += 1

        for (p_og, c_og), n in sorted(edge_count.items()):
            link_rows.append({
                "hg":           args.hg,
                "parent_level": parent_lv,
                "parent_og":    p_og,
                "child_level":  child_lv,
                "child_og":     c_og,
                "n_genes":      n,
            })

    pd.DataFrame(
        link_rows,
        columns=["hg", "parent_level", "parent_og", "child_level", "child_og", "n_genes"],
    ).to_csv(args.output_links, sep="\t", index=False)

    print(
        f"HG {args.hg}: {len(stat_rows)} OGs across "
        f"{len(level_data)}/{len(ordered_levels)} levels, "
        f"{len(link_rows)} cross-level edges",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
