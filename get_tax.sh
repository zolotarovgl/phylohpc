#!/usr/bin/env bash
if [ $# -ne 1 ]; then
    echo "Usage: $0 species_list.txt"
    exit 1
fi
LIST=$1
taxonkit name2taxid < $LIST \
| cut -f2 \
| taxonkit lca \
| taxonkit reformat --format newick --show-name
