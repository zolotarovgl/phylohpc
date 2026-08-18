#!/usr/bin/env python3
"""Print the number of multifurcating (>2 children) nodes in a newick file.

GeneRax requires a strictly binary species tree and says only "Cannot parse species tree
file" when it is not -- which cost four days and 475 dead tasks on 2026-08-14. This is the
single implementation of that check, used by workflow/preflight.sh, by step2.nf and by the
test suite, so the three cannot disagree.
"""
import sys


def count_multifurcations(newick: str) -> int:
    kids, bad = [], 0
    for ch in newick:
        if ch == '(':
            kids.append(1)
        elif ch == ',' and kids:
            kids[-1] += 1
        elif ch == ')':
            n = kids.pop() if kids else 0
            if n > 2:
                bad += 1
    return bad


if __name__ == '__main__':
    print(count_multifurcations(open(sys.argv[1]).read().strip()))
