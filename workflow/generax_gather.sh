#!/usr/bin/env bash
set -euo pipefail

# Usage message
usage() {
  echo "Usage: $0 <LIST_DONE> <LIST_RUN>"
  echo
  echo "  LIST_DONE  Path to file that will list prefixes of completed HG runs"
  echo "  LIST_RUN   Path to file that will list prefixes to run generax on"
  exit 1
}

# Check for exactly two arguments
if [ $# -ne 2 ]; then
  usage
fi

source workflow/functions.sh
read_config configs/config.txt


# Assign command-line arguments
LIST_DONE=$1
LIST_RUN=$2

# ensure the directory for LIST_DONE exists
dir_done=$(dirname "$LIST_DONE")
mkdir -p "$dir_done"

# turn on nullglob so patterns with no matches expand to empty arrays
shopt -s nullglob

# clear out any old lists
> "$LIST_DONE"



# loop over every .newick (if any) and extract the 3rd "/"-field
for file in ${GENERAX_DIR}/*/results/*/*.newick; do
  # file is like "results/generax/HG123/results/.../file.newick"
  # extract the "HG123" prefix
  prefix="${file%%/results/*}"    # strip from "/results/..." onward
  prefix="${prefix##*/}"          # strip up to last "/"
  echo "$prefix" >> "$LIST_DONE"
done

# build LIST_RUN by taking all iqtree prefixes minus those done
grep_opts=(-F -x -v -f "$LIST_DONE")
ls ${TREE_DIR}/*.iqtree \
  | xargs -n1 basename -s .iqtree \
  | grep "${grep_opts[@]}" \
  > "$LIST_RUN"

# count and report
N_DONE=$(wc -l < "$LIST_DONE")
N_RUN=$(wc -l < "$LIST_RUN")
echo "$0: $LIST_DONE: $N_DONE HGs done"
echo "$0: $LIST_RUN: $N_RUN HGs to run"
