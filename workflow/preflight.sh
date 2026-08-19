#!/usr/bin/env bash
# Pre-flight for a phylohpc run. Reads the params YAML, checks every input, and fails
# with a NAMED cause in about a second.
#
# It exists because the 2026-08-14 run staged a multifurcating species tree, GeneRax said
# only "Cannot parse species tree file", the errorStrategy relabelled that as a family
# parsing error, and the run reported success with an empty output directory.
set -uo pipefail
YAML="${1:?usage: preflight.sh <params.yaml> <step>}"
STEP="${2:-step2}"
cd "$(dirname "$0")/.."

get() { python3 -c "import sys,yaml;print(yaml.safe_load(open(sys.argv[1])).get(sys.argv[2],'') or '')" "$YAML" "$1"; }

fail=0
say()  { printf '   %s\n' "$*"; }
err()  { printf 'PREFLIGHT ERROR: %s\n' "$*" >&2; fail=1; }

[[ -f "$YAML" ]] || { err "params file not found: $YAML"; exit 1; }
say "params file : $YAML"

for key in infasta family_info refnames; do
  v=$(get "$key")
  [[ -n "$v" ]] || { err "'$key' is not set in $YAML"; continue; }
  [[ -e "$v" ]] || err "'$key' points at a missing file: $v"
done

tree=$(get species_tree)
generax=$(get run_generax)
if [[ "$generax" == "True" || "$generax" == "true" || "$generax" == "1" ]]; then
  [[ -e "$tree" ]] || err "run_generax is on but species_tree is missing: $tree"
  if [[ -e "$tree" ]]; then
    n=$(python3 workflow/count_multifurcations.py "$tree")
    if [[ "$n" != "0" ]]; then
      err "species tree '$tree' has $n multifurcating node(s); GeneRax needs a strictly binary tree.
  Resolve it:  python workflow/check_tree.py '$tree' <ids> <out.newick> --random-resolve --seed 1
  then set species_tree to that file. (This is what killed the 2026-08-14 run.)"
    else
      say "species tree: $tree (strictly binary)"
    fi
  fi
  command -v generax >/dev/null || say "WARNING: generax not on PATH -- the bioconda build is non-MPI; load an MPI build"
  command -v mpirun  >/dev/null || say "WARNING: mpirun not on PATH"
fi

[[ $fail -eq 0 ]] || exit 1
say "preflight OK ($STEP)"
