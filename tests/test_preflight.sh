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
