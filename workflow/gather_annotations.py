#!/usr/bin/env python3

import argparse
import subprocess
import os
import sys
import tempfile
from pathlib import Path


def run_checked(cmd, stdout=None):
    """Run subprocess safely and fail loudly if it errors."""
    return subprocess.run(cmd, check=True, text=True, stdout=stdout)


def parse_args():
    parser = argparse.ArgumentParser(description="Gather species annotation")

    parser.add_argument("--id", required=True, help="Species prefix")
    parser.add_argument("--outfile", help="Optional output file")
    parser.add_argument("--search-dir", help="Search directory")

    parser.add_argument(
        "--tree-dir",
        nargs="+",
        help="Directories with POSSVM outputs (searched in order of priority)"
    )

    parser.add_argument("--prefix", default=None, help="Optional family prefix")

    return parser.parse_args()


def collect_tmp_anno(tree_dirs, prefix, species_id):
    """
    Collect annotations from POSSVM outputs.
    Directories are searched in priority order.
    """

    result = {}

    for tree_dir in tree_dirs:
        tree_dir = Path(tree_dir)

        if not tree_dir.exists():
            continue

        if not prefix:
            pattern = "*groups.csv"
        else:
            pattern = f"{prefix}*groups.csv"

        for file in tree_dir.glob(pattern):

            with open(file) as f:
                for line in f:
                    if species_id in line:
                        parts = line.rstrip("\n").split("\t")

                        if len(parts) >= 4:
                            gene = parts[0]

                            # keep first occurrence (priority order)
                            if gene not in result:
                                result[gene] = parts[1:4]

    return result


def collect_ids_todo(search_dir, prefix, species_id):
    ids = set()

    if not prefix:
        pattern = "*genes.list"
    else:
        pattern = f"{prefix}*.genes.list"

    for file in Path(search_dir).glob(pattern):
        with open(file) as f:
            for line in f:
                if species_id in line:
                    ids.add(line.strip().split()[0])

    return sorted(ids)


def collect_gene2class(search_dir, prefix, species_id):
    gene2class = {}

    if not prefix:
        pattern = "*genes.list"
    else:
        pattern = f"{prefix}*.genes.list"

    for file in Path(search_dir).glob(pattern):

        pref = file.name.replace(".genes.list", "")

        with open(file) as f:
            for line in f:
                if species_id in line:
                    gene = line.strip().split()[0]
                    gene2class[gene] = pref

    return gene2class


def collect_pep2hg(cluster_dir, species_id):
    """
    Pure Python reimplementation of get_gene2cluster.sh
    """

    cluster_dir = Path(cluster_dir)

    if not cluster_dir.exists():
        raise FileNotFoundError(f"Cluster directory not found: {cluster_dir}")

    mapping = {}

    for file in cluster_dir.glob("*_cluster.tsv"):

        pref = file.name.replace("_cluster.tsv", "")

        with open(file) as f:
            for line in f:

                parts = line.rstrip("\n").split("\t")

                if len(parts) < 2:
                    continue

                cluster_id = parts[0]
                gene = parts[1]

                if gene.startswith(species_id):
                    mapping[gene] = f"{pref}.{cluster_id}"

    return mapping


def build_result(tmp_anno, ids_todo, gene2class, pep2hg):

    results = []

    for gene in ids_todo:

        if gene in tmp_anno:
            class_info = tmp_anno[gene]
        else:
            class_info = ["Unclassified"]

        hg = pep2hg.get(gene)

        if class_info[0] == "Unclassified" and hg:
            class_info = [f"{hg}:Unclassified"]

        gene_class = gene2class.get(gene)

        if class_info[0] == "Unclassified" and gene_class:
            class_info = [f"{gene_class}:Unclassified"]

        row = [gene] + class_info

        results.append("\t".join(row))

    return results


def main():

    args = parse_args()

    species_id = args.id
    prefix = args.prefix

    if not args.search_dir or not args.tree_dir:
        sys.exit("Provide both --search-dir and --tree-dir")

    search_dir = args.search_dir
    tree_dirs = args.tree_dir
    cluster_dir = "results/clusters"

    with tempfile.TemporaryDirectory(prefix="gather_anno_"):

        tmp_anno = collect_tmp_anno(tree_dirs, prefix, species_id)

        ids_todo = collect_ids_todo(search_dir, prefix, species_id)

        gene2class = collect_gene2class(search_dir, prefix, species_id)

        pep2hg = collect_pep2hg(cluster_dir, species_id)

        result = build_result(tmp_anno, ids_todo, gene2class, pep2hg)

        if args.outfile:
            with open(args.outfile, "w") as out:
                out.write("\n".join(result) + "\n")
        else:
            print("\n".join(result))


if __name__ == "__main__":
    main()