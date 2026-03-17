#!/usr/bin/env python3
"""
Build a presence/absence matrix (PAM) from POSSVM ortholog_groups.csv files.

Input:
  - One or more POSSVM *.ortholog_groups.csv files (tab-separated, no header)
    Columns: gene_id  og_name  reference_og  reference_gene_name
  - A file listing the in-clade species of interest (one prefix per line)

Output:
  - TSV where rows = species, columns = orthogroup IDs, values = 0/1
"""

import argparse
import sys
from pathlib import Path

import pandas as pd


def get_species_prefix(gene_id: str) -> str:
    """Return the species prefix (everything before the first '_' or '.')."""
    for sep in ("_", "."):
        if sep in gene_id:
            return gene_id.split(sep)[0]
    return gene_id  # single-token gene ID — return as-is


def load_gene_og_table(csv_paths: list[str], in_species: set[str]) -> pd.DataFrame:
    """
    Parse all POSSVM CSV files and return a DataFrame with columns
    [gene, og, species].  Only genes from in_species are retained.
    """
    records = []
    for path in csv_paths:
        p = Path(path)
        if not p.exists() or p.stat().st_size == 0:
            continue
        with p.open() as fh:
            for line in fh:
                line = line.rstrip("\n")
                if not line or line.startswith("#"):
                    continue
                parts = line.split("\t")
                if len(parts) < 2:
                    continue
                gene = parts[0]
                og   = parts[1]
                sps  = get_species_prefix(gene)
                if sps in in_species:
                    records.append((gene, og, sps))

    if not records:
        return pd.DataFrame(columns=["gene", "og", "species"])

    return pd.DataFrame(records, columns=["gene", "og", "species"])


def main():
    parser = argparse.ArgumentParser(
        description="Build a binary presence/absence matrix from POSSVM OG CSV files."
    )
    parser.add_argument(
        "--csvs", nargs="+", required=True,
        help="POSSVM *.ortholog_groups.csv files"
    )
    parser.add_argument(
        "--species", required=True,
        help="File listing in-clade species prefixes (one per line)"
    )
    parser.add_argument(
        "--min_presence", type=int, default=2,
        help="Minimum number of species that must carry an OG for it to be retained (default: 2)"
    )
    parser.add_argument(
        "--output", required=True,
        help="Output PAM file (TSV, species × OGs)"
    )
    args = parser.parse_args()

    # Load in-clade species list
    with open(args.species) as fh:
        in_species = [l.strip() for l in fh if l.strip()]
    in_species_set = set(in_species)

    print(f"In-clade species: {len(in_species)}", file=sys.stderr)

    # Parse all CSV files
    df = load_gene_og_table(args.csvs, in_species_set)

    if df.empty:
        print(
            "WARNING: No gene–OG mappings found for in-clade species. "
            "Writing empty PAM.",
            file=sys.stderr,
        )
        empty_pam = pd.DataFrame(index=in_species, dtype=int)
        empty_pam.index.name = "species"
        empty_pam.to_csv(args.output, sep="\t")
        return

    # Binary PAM: species × OG → 1 if species has ≥1 gene in that OG
    pam = (
        df.groupby(["species", "og"])
        .size()
        .unstack(fill_value=0)
        .gt(0)
        .astype(int)
    )

    # Reindex so every in-clade species appears (even those with no hits)
    pam = pam.reindex(in_species, fill_value=0)
    pam.index.name = "species"

    # Filter columns by minimum presence
    og_counts = pam.sum(axis=0)
    n_before = pam.shape[1]
    pam = pam.loc[:, og_counts >= args.min_presence]
    n_after = pam.shape[1]

    print(
        f"PAM: {pam.shape[0]} species × {n_after} OGs "
        f"(filtered {n_before - n_after} OGs with <{args.min_presence} species)",
        file=sys.stderr,
    )

    pam.to_csv(args.output, sep="\t")
    print(f"Written: {args.output}", file=sys.stderr)


if __name__ == "__main__":
    main()
