#!/usr/bin/env python3

import argparse
import sys
import random
from ete3 import Tree


def parse_args():
    parser = argparse.ArgumentParser(
        description="Prune a Newick tree to given IDs and ensure it is strictly binary."
    )
    parser.add_argument("infile", help="Input Newick tree")
    parser.add_argument("ids", help="File with one tip ID per line")
    parser.add_argument("outfile", help="Output Newick tree")
    parser.add_argument(
        "--random-resolve",
        action="store_true",
        help="Randomly resolve multifurcating nodes into a strictly binary tree",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=None,
        help="Random seed for reproducible random resolution",
    )
    return parser.parse_args()


def collapse_unary(node):
    """
    Collapse nodes with a single child.
    """
    for child in list(node.children):
        collapse_unary(child)

    if not node.is_root() and len(node.children) == 1:
        node.delete(prevent_nondicotomic=False)


def random_resolve_polytomies(tree):
    """
    Randomly resolve any multifurcations so all internal nodes become binary.
    """
    # Shuffle children at each multifurcating node so resolution order is random
    for node in tree.traverse("postorder"):
        if len(node.children) > 2:
            random.shuffle(node.children)
            node.resolve_polytomy(recursive=False)


def check_strict_binary(tree):
    """
    Ensure all internal nodes have exactly 2 children.
    """
    for node in tree.traverse():
        if not node.is_leaf() and len(node.children) != 2:
            raise ValueError(
                f"Tree is not strictly binary. "
                f"Node '{node.name}' has {len(node.children)} children."
            )


def main():
    args = parse_args()

    if args.seed is not None:
        random.seed(args.seed)

    # Load IDs
    with open(args.ids) as f:
        ids = {line.strip() for line in f if line.strip()}

    # Read tree
    tree = Tree(args.infile, format=1)

    tips = {leaf.name for leaf in tree.iter_leaves()}
    found = len(ids & tips)
    print(f"{found}/{len(ids)} ids found in the tree")

    missing = ids - tips
    if missing:
        print("Missing ids in tips:", ",".join(sorted(missing)))

    # Prune tree
    tree.prune(ids, preserve_branch_length=True)

    # Collapse unary nodes created by pruning
    collapse_unary(tree)

    # Randomly resolve multifurcations if requested
    if args.random_resolve:
        random_resolve_polytomies(tree)

    # Check strict binarity
    try:
        check_strict_binary(tree)
    except ValueError as e:
        print("ERROR:", e)
        sys.exit(1)

    # Write output
    tree.write(outfile=args.outfile,format = 1)
    print("Tree successfully pruned and verified as strictly binary.")


if __name__ == "__main__":
    main()