#!/usr/bin/env bash
set -uo pipefail
source "$(dirname "$0")/lib.sh"
cd "$(dirname "$0")/.."
tmp=$(mktemp -d)

cat > "$tmp/ok.yaml" <<Y
infasta: data/input.fasta
family_info: data/gene_families_searchinfo.csv
refnames: data/Mmus_gene_names.csv
species_tree: tests/data/binary.newick
run_generax: false
outdir: $tmp/results
Y
assert_ok bash workflow/preflight.sh "$tmp/ok.yaml" step2

cat > "$tmp/bad_tree.yaml" <<Y
infasta: data/input.fasta
family_info: data/gene_families_searchinfo.csv
refnames: data/Mmus_gene_names.csv
species_tree: tests/data/multifurcating.newick
run_generax: true
outdir: $tmp/results
Y
assert_fails bash workflow/preflight.sh "$tmp/bad_tree.yaml" step2
out=$(bash workflow/preflight.sh "$tmp/bad_tree.yaml" step2 2>&1 || true)
assert_contains "$out" "multifurcating" "names the actual problem"
assert_contains "$out" "check_tree.py"  "tells you how to fix it"

# The wrapper and the in-workflow guard must produce the SAME sentence, not just each
# contain the same two substrings -- that weaker check is exactly what let the two
# messages drift apart undetected (2026-08-19 review). Both read from the same template
# (workflow/messages/species_tree_multifurcating.txt); the only thing allowed to differ
# is how the tree path is spelled (preflight.sh keeps it relative, step2.nf's `file()`
# resolves it absolute), so normalise only that before comparing.
nf_out=$(nextflow run step2.nf -preview -params-file "$tmp/bad_tree.yaml" 2>&1 || true)
pf_msg=$(echo "$out"    | grep -o 'species tree.*multifurcating node.*' | head -1)
nf_msg=$(echo "$nf_out" | grep -o 'species tree.*multifurcating node.*' | head -1)
pf_msg_norm=$(echo "$pf_msg" | sed -E 's#[^ ]*multifurcating\.newick#TREEPATH#g')
nf_msg_norm=$(echo "$nf_msg" | sed -E 's#[^ ]*multifurcating\.newick#TREEPATH#g')
assert_ok test -n "$pf_msg_norm"
assert_eq "$pf_msg_norm" "$nf_msg_norm" "preflight.sh and step2.nf report the IDENTICAL multifurcating-tree message (tree path normalised)"

cat > "$tmp/missing.yaml" <<Y
infasta: does/not/exist.fasta
family_info: data/gene_families_searchinfo.csv
refnames: data/Mmus_gene_names.csv
species_tree: tests/data/binary.newick
run_generax: false
outdir: $tmp/results
Y
assert_fails bash workflow/preflight.sh "$tmp/missing.yaml" step2
rm -rf "$tmp"
