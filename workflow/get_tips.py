import sys
from Bio import Phylo

if len(sys.argv) < 2:
    print("Usage: tips.py treefile")
    sys.exit(1)

tree = Phylo.read(sys.argv[1], "newick")
for tip in tree.get_terminals():
    print(tip.name)

