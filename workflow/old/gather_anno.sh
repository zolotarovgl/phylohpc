#!/usr/bin/env bash

SEARCH_DIR="results/search"
TREE_DIR="results/gene_trees"
TREE_DIR="results/generax"
ID=""
OUTFILE=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        --id) ID="$2"; shift 2;;
        --outfile) OUTFILE="$2"; shift 2;;
        --search-dir) SEARCH_DIR="$2"; shift 2;;
        --tree-dir) TREE_DIR="$2"; shift 2;;
        *) echo "Unknown argument: $1"; exit 1;;
    esac
done

if [[ -z "$ID" ]]; then
    echo "Usage: $0 --id ID [--outfile OUTFILE] [--search-dir DIR] [--tree-dir DIR]"
    exit 1
fi

cat "${TREE_DIR}"/*groups.csv | grep "$ID" > tmp_anno
cat "${SEARCH_DIR}"/*.genes.list | grep "$ID" | sort | uniq > ids_todo

# Gene to Class Mapping 
> gene2class
for FILE in ${SEARCH_DIR}/*.genes.list; do
    PREF=$(basename "$FILE" | sed 's/.genes.list//')
    grep "$ID" "$FILE" | awk -v PREF=$PREF '{print $1"\t"PREF}'
done >> gene2class
#RESULT=$(awk 'NR==FNR{dict[$1]=$2"\t"$3"\t"$4;next}{if($1 in dict) print $1"\t"dict[$1]; else print $1"\tUnclassified"}' tmp_anno ids_todo | awk 'BEGIN{OFS="\t"}FNR==NR{d[$1]=$2;next}{if($2=="Unclassified"){$2=d[$1]":Unclassified"};print $0}' gene2class -)


bash workflow/get_gene2cluster.sh configs/config.txt pep2hg
grep -E "^${ID}" pep2hg > tmp_pep2hg; mv tmp_pep2hg pep2hg




#RESULT=$(awk 'NR==FNR{dict[$1]=$2"\t"$3"\t"$4;next}{if($1 in dict) print $1"\t"dict[$1]; else print $1"\tUnclassified"}' tmp_anno ids_todo | awk 'BEGIN{OFS="\t"}FNR==NR{d[$1]=$2;next}{if($2=="Unclassified"){$2=d[$1]":Unclassified"};print $0}' pep2hg -)

RESULT=$(awk 'NR==FNR{dict[$1]=$2"\t"$3"\t"$4;next}{if($1 in dict) print $1"\t"dict[$1]; else print $1"\tUnclassified"}' tmp_anno ids_todo | awk 'BEGIN{OFS="\t"}FNR==NR{d[$1]=$2;next}{if($2=="Unclassified"){$2=d[$1]":Unclassified"};print $0}' pep2hg - | awk 'FNR==NR{d[$1]=$2;next}{if($2==":Unclassified"){$2 = d[$1]":Unclassified"};print $1"\t"$2"\t"$3"\t"$4}' gene2class -)


if [[ -n "$OUTFILE" ]]; then
    echo "$RESULT" > "$OUTFILE"
else
    echo "$RESULT"
fi

rm gene2class pep2hg tmp_anno ids_todo
exit 0

