#!/usr/bin/env python3
"""Clean POSSVM's orthogroup_support column, which reports a placeholder as if it were data.

POSSVM sets a group's support to the support of the MRCA of its member genes
(possvm.py::find_support_cluster). Two cases never touch an internal node, so ete3 hands back
a DEFAULT of 1.0 and it is indistinguishable from a real, perfect bootstrap:

  * singleton group (n == 1)  -- there is no ancestor; the "MRCA" is the leaf itself
  * MRCA == tree root         -- the group spans the root split, so the value is the root's
                                 own default

Measured over 400 families / 3021 groups of results/possvm_prev (2026-08-28):
  singleton                    466  (15.4 %)
  MRCA == root                 144  ( 4.8 %)
  measured on an internal node 2411  (79.8 %), of which 19 are GENUINELY 1.0

Those 19 are why this works off the recomputed MRCA and never off the string "1.0": a blanket
rewrite of 1.0 would destroy 19 real (very poor) bootstrap values.

What it does to <prefix>.ortholog_groups.csv:
  * DROPS rows belonging to a singleton group -- a one-gene orthogroup is not an orthology
    statement, and its support is undefined
  * sets orthogroup_support = -1 where the group's MRCA is the tree root
  * leaves every measured value untouched, including a real 1.0

⚠ Only the CSV is rewritten. <prefix>.ortholog_groups.newick still carries the dropped
singleton leaves, so gene counts taken from the tree and from the table can differ.

Usage:  possvm_postprocess.py <prefix>.ortholog_groups.csv [--newick FILE] [--dry-run]
"""
import argparse
import collections
import os
import sys

ROOT_SENTINEL = "-1"


def load_tree(path):
    from ete3 import Tree
    t = Tree(path, format=1)
    # POSSVM decorates leaf labels as "<gene> | <og> | <ref> | "
    for lf in t.get_leaves():
        lf.name = lf.name.split(" | ")[0].strip()
    return t


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("csv")
    ap.add_argument("--newick", default=None,
                    help="default: the .newick beside the csv")
    ap.add_argument("--dry-run", action="store_true")
    a = ap.parse_args()

    nwk = a.newick or (a.csv[:-4] + ".newick")
    if not os.path.exists(a.csv):
        sys.exit(f"no such file: {a.csv}")
    if not os.path.exists(nwk):
        sys.exit(f"no newick beside the csv: {nwk}")

    with open(a.csv) as fh:
        lines = [l.rstrip("\n").split("\t") for l in fh]
    if not lines:
        sys.exit(f"empty: {a.csv}")
    header, rows = lines[0], [r for r in lines[1:] if len(r) >= 3]
    try:
        i_gene, i_og, i_sup = (header.index(c) for c in
                               ("gene", "orthogroup", "orthogroup_support"))
    except ValueError:
        sys.exit(f"unexpected header in {a.csv}: {header}")

    groups = collections.defaultdict(list)
    for r in rows:
        groups[r[i_og]].append(r[i_gene])

    t = load_tree(nwk)
    leaves = set(t.get_leaf_names())
    genes = {r[i_gene] for r in rows}
    missing = genes - leaves
    # House rule: assert the intersection, never intersect silently.
    if len(missing) > 0.01 * max(len(genes), 1):
        sys.exit(f"gene ids do not match the tree: {len(missing)}/{len(genes)} absent "
                 f"({sorted(missing)[:3]}) -- check the leaf-label format in {nwk}")

    singleton = {og for og, m in groups.items() if len(m) == 1}
    rooted = set()
    for og, mem in groups.items():
        if og in singleton or len(set(mem) & leaves) < 2:
            continue
        try:
            anc = t.get_common_ancestor([g for g in mem if g in leaves])
        except Exception:
            continue
        if anc is t:
            rooted.add(og)

    out = [header]
    dropped = flagged = 0
    for r in rows:
        if r[i_og] in singleton:
            dropped += 1
            continue
        if r[i_og] in rooted:
            r = list(r)
            r[i_sup] = ROOT_SENTINEL
            flagged += 1
        out.append(r)

    print(f"{os.path.basename(a.csv)}: groups={len(groups)} "
          f"singleton_groups={len(singleton)} (rows dropped={dropped}) "
          f"root_MRCA_groups={len(rooted)} (rows set to {ROOT_SENTINEL}={flagged}) "
          f"rows {len(rows)} -> {len(out)-1}")

    if a.dry_run:
        return
    tmp = a.csv + ".tmp"
    with open(tmp, "w") as fh:
        for r in out:
            fh.write("\t".join(r) + "\n")
    os.replace(tmp, a.csv)


if __name__ == "__main__":
    main()
