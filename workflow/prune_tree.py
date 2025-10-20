import sys
from Bio import Phylo

if len(sys.argv) < 4:
    print('script.py in.tree ids out.tree')
    sys.exit(1)

infile = sys.argv[1]
ids_file = sys.argv[2]
outfile = sys.argv[3]

with open(ids_file) as f:
    ids = {x.strip() for x in f if x.strip()}

tree = Phylo.read(infile, 'newick')
tips = [t.name for t in tree.get_terminals()]
found = sum(1 for i in ids if i in tips)
print(f'{found}/{len(ids)} ids found in the tree')
missing = [i for i in ids if i not in tips]
if missing:
    print('Missing ids in tips:', ','.join(missing))

to_remove = [t for t in tree.get_terminals() if t.name not in ids]
for t in to_remove:
    tree.prune(t)

# collapse any nodes with only one child
def collapse_unary(clade):
    while len(clade.clades) == 1:
        child = clade.clades[0]
        clade.name = child.name
        clade.branch_length = child.branch_length
        clade.clades = child.clades
    for sub in clade.clades:
        collapse_unary(sub)

collapse_unary(tree.root)
Phylo.write(tree, outfile, 'newick')

