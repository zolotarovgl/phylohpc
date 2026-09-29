#!/usr/bin/env python3
"""Drop singleton orthogroups from a POSSVM ortholog_groups.csv.

A one-gene orthogroup is not an orthology statement, so its rows are removed.

Support placeholders are no longer handled here. POSSVM itself (zolotarovgl fork,
fix/root-mrca-indexerror, a41f982) now writes orthogroup_support = NA whenever the
group's MRCA has no measured support: a leaf (singleton), the root, or a node added
by polytomy resolution / rerooting. Earlier versions of this script rewrote the
root case to -1; that sentinel is gone -- read NA as "not measured".

For scale, in results/possvm_prev (28/08/2026, 400 families / 3021 groups):
singletons 466 (15.4 %), MRCA == root 144 (4.8 %).

⚠ Only the CSV is rewritten. <prefix>.ortholog_groups.newick still carries the dropped
singleton leaves, so gene counts taken from the tree and from the table can differ.

Usage:  possvm_postprocess.py <prefix>.ortholog_groups.csv [--dry-run]
"""
import argparse
import collections
import os
import sys


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("csv")
    ap.add_argument("--dry-run", action="store_true")
    a = ap.parse_args()

    if not os.path.exists(a.csv):
        sys.exit(f"no such file: {a.csv}")

    with open(a.csv) as fh:
        lines = [l.rstrip("\n").split("\t") for l in fh]
    if not lines:
        sys.exit(f"empty: {a.csv}")
    header, rows = lines[0], [r for r in lines[1:] if len(r) >= 3]
    try:
        i_gene, i_og = (header.index(c) for c in ("gene", "orthogroup"))
    except ValueError:
        sys.exit(f"unexpected header in {a.csv}: {header}")

    groups = collections.defaultdict(list)
    for r in rows:
        groups[r[i_og]].append(r[i_gene])
    singleton = {og for og, m in groups.items() if len(m) == 1}

    out = [header] + [r for r in rows if r[i_og] not in singleton]

    print(f"{os.path.basename(a.csv)}: groups={len(groups)} "
          f"singleton_groups={len(singleton)} rows {len(rows)} -> {len(out)-1}")

    if a.dry_run:
        return
    tmp = a.csv + ".tmp"
    with open(tmp, "w") as fh:
        for r in out:
            fh.write("\t".join(r) + "\n")
    os.replace(tmp, a.csv)


if __name__ == "__main__":
    main()
