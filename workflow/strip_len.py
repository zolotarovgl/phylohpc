import sys
from Bio import Phylo
from io import StringIO

if len(sys.argv) < 3:
    print('Usage: strip_len.py in.tree out.tree')
    sys.exit(1)

tree = Phylo.read(sys.argv[1], 'newick')
for clade in tree.find_clades():
    clade.branch_length = None

buf = StringIO()
Phylo.write(tree, buf, 'newick')
text = buf.getvalue()
text = text.replace(':0.00000', '')
text = text.replace(':0.0', '')
text = text.replace(':0', '')

with open(sys.argv[2], 'w') as f:
    f.write(text)

