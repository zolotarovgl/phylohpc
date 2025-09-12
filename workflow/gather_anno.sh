#!/usr/bin/env bash

# defaults
SEARCH_DIR="results/search"
TREE_DIR="results/gene_trees"
ID=""
OUTFILE=""

# parse named args
while [[ $# -gt 0 ]]; do
    case "$1" in
        --id) ID="$2"; shift 2;;
        --outfile) OUTFILE="$2"; shift 2;;
        --search-dir) SEARCH_DIR="$2"; shift 2;;
        --tree-dir) TREE_DIR="$2"; shift 2;;
        *) echo "Unknown argument: $1"; exit 1;;
    esac
done

if [[ -z "$ID" || -z "$OUTFILE" ]]; then
    echo "Usage: $0 --id ID --outfile OUTFILE [--search-dir DIR] [--tree-dir DIR]"
    exit 1
fi

cat "${TREE_DIR}"/*groups.csv | grep "$ID" > tmp_anno
cat "${SEARCH_DIR}"/*.genes.list | grep "$ID" | sort | uniq > ids_todo

> gene2class; for FILE in $(ls ${SEARCH_DIR}/*.genes.list); do PREF=$(basename $FILE | sed 's/.genes.list//'); grep $ID $FILE | awk -v PREF=$PREF '{print $1"\t"PREF}'; done >> gene2class

awk 'NR==FNR{dict[$1]=$2"\t"$3"\t"$4;next}{if($1 in dict) print $1"\t"dict[$1]; else print $1"\tUnclassified"}' tmp_anno ids_todo | awk 'BEGIN{OFS="\t"}FNR==NR{d[$1]=$2;next}{if($2=="Unclassified"){$2=d[$1]":Unclassified"};print $0}' gene2class - > "$OUTFILE"

rm gene2class tmp_anno ids_todo
exit 0
