#!/usr/bin/env python3
"""
Extract species belonging to a named clade from a species tree.

Given a full species tree with named internal nodes and a target node name,
outputs:
  - {prefix}.in_species.txt     : species within the clade (one per line)
  - {prefix}.ignore_species.txt : species outside the clade (for POSSVM --ignoretips)
  - {prefix}.pruned.tree        : newick subtree rooted at the target node
"""

import argparse
import sys

from ete3 import Tree


def main():
    parser = argparse.ArgumentParser(
        description="Extract in-clade and out-of-clade species lists from a named node."
    )
    parser.add_argument("--tree", required=True, help="Species tree (newick) with named internal nodes")
    parser.add_argument("--node", required=True, help="Internal node name (clade to analyse)")
    parser.add_argument("--out_prefix", required=True, help="Prefix for output files")
    args = parser.parse_args()

    # Load tree preserving internal node names (format=1 keeps NHX + named internal nodes)
    tree = Tree(args.tree, format=1)

    # Find the target node
    hits = tree.search_nodes(name=args.node)
    if not hits:
        named = sorted(set(n.name for n in tree.traverse() if n.name and not n.is_leaf()))
        print(
            f"ERROR: Node '{args.node}' not found in tree.\n"
            f"Available named internal nodes ({len(named)}):\n  " + "\n  ".join(named),
            file=sys.stderr,
        )
        sys.exit(1)

    if len(hits) > 1:
        print(
            f"WARNING: Node name '{args.node}' is ambiguous ({len(hits)} matches). "
            "Using the first hit.",
            file=sys.stderr,
        )

    target = hits[0]

    in_species  = sorted(target.get_leaf_names())
    all_species = sorted(tree.get_leaf_names())
    out_species = [s for s in all_species if s not in set(in_species)]

    print(
        f"Node '{args.node}': {len(in_species)} in-clade, {len(out_species)} out-of-clade species",
        file=sys.stderr,
    )
    if len(in_species) < 2:
        print(
            f"ERROR: Clade '{args.node}' has fewer than 2 species — nothing to reconstruct.",
            file=sys.stderr,
        )
        sys.exit(1)

    # Write species lists
    with open(f"{args.out_prefix}.in_species.txt", "w") as fh:
        fh.write("\n".join(in_species) + "\n")

    with open(f"{args.out_prefix}.ignore_species.txt", "w") as fh:
        fh.write("\n".join(out_species) + "\n")

    # Write pruned subtree
    # Use a detached copy so the original tree is not modified
    pruned = target.copy("deepcopy")
    pruned.write(format=1, outfile=f"{args.out_prefix}.pruned.tree")

    print(
        f"Outputs: {args.out_prefix}.in_species.txt  "
        f"{args.out_prefix}.ignore_species.txt  "
        f"{args.out_prefix}.pruned.tree",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
