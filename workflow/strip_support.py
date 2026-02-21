import sys
from Bio import Phylo

if len(sys.argv) < 3:
    print("Usage: strip_bootstrap.py in.tree out.tree")
    sys.exit(1)

infile = sys.argv[1]
outfile = sys.argv[2]

tree = Phylo.read(infile, "newick")
for clade in tree.find_clades():
    clade.confidence = None
    clade.name = None
Phylo.write(tree, outfile, "newick")
