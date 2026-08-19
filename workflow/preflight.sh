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

for key in infasta refnames; do
  v=$(get "$key")
  [[ -n "$v" ]] || { err "'$key' is not set in $YAML"; continue; }
  [[ -e "$v" ]] || err "'$key' points at a missing file: $v"
done

# IMPORTANT 7 (2026-08-19 review): step2.nf accepts EITHER 'family_info' or its alias
# 'genefam_info' (params.family_info's fallback chain), but this pre-flight required
# 'family_info' only and named just that key on failure -- so an existing params file
# using genefam_info was rejected by an unbypassable gate whose own error message named
# the wrong key. Accept the same aliases step2.nf does, and name both in the message.
fi_val=$(get family_info)
gi_val=$(get genefam_info)
if [[ -n "$fi_val" ]]; then
  [[ -e "$fi_val" ]] || err "'family_info' points at a missing file: $fi_val"
elif [[ -n "$gi_val" ]]; then
  [[ -e "$gi_val" ]] || err "'genefam_info' points at a missing file: $gi_val"
else
  err "neither 'family_info' nor its alias 'genefam_info' is set in $YAML"
fi

tree=$(get species_tree)
generax=$(get run_generax)
if [[ "$generax" == "True" || "$generax" == "true" || "$generax" == "1" ]]; then
  [[ -e "$tree" ]] || err "run_generax is on but species_tree is missing: $tree"
  if [[ -e "$tree" ]]; then
    # MINOR 8 (2026-08-19 review): if count_multifurcations.py errors, $n was previously
    # left empty and the user was told the tree has "" multifurcating nodes -- check the
    # exit status and report that failure by name instead.
    cm_errf=$(mktemp)
    n=$(python3 workflow/count_multifurcations.py "$tree" 2> "$cm_errf")
    cm_rc=$?
    if [[ $cm_rc -ne 0 ]]; then
      err "count_multifurcations.py failed (exit $cm_rc) on $tree: $(cat "$cm_errf")"
      rm -f "$cm_errf"
    elif [[ "$n" != "0" ]]; then
      # ONE message, defined in workflow/messages/species_tree_multifurcating.txt and
      # substituted here -- step2.nf substitutes the SAME file the SAME way, so the
      # failure reads identically whichever way the pipeline was launched. Do not
      # hand-maintain this sentence in two places again (see the 2026-08-19 review that
      # caught the two copies drifting).
      msg=$(sed -e "s|__TREE__|$tree|g" -e "s|__N__|$n|g" \
                "workflow/messages/species_tree_multifurcating.txt")
      err "$msg"
    else
      say "species tree: $tree (strictly binary)"
    fi
  fi
  command -v generax >/dev/null || say "WARNING: generax not on PATH -- the bioconda build is non-MPI; load an MPI build"
  command -v mpirun  >/dev/null || say "WARNING: mpirun not on PATH"
fi

[[ $fail -eq 0 ]] || exit 1
say "preflight OK ($STEP)"
