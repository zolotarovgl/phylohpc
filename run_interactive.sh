#!/usr/bin/env bash
# Launch phylohpc from an INTERACTIVE HPC session (no sbatch, no wrapper).
#
#   salloc -p genoa64 -c 8 --mem 16G -t 8:00:00     # or srun --pty bash
#   ./run_interactive.sh step2 --run_generax
#
# Written 2026-08-18 after the GeneRax post-mortem: the failure was invisible because the
# run was launched through SLURM + an errorStrategy that relabelled the error, and the
# real message sat in a work dir on /no_backup. This script runs the pipeline in the
# foreground, validates its inputs BEFORE launching, and tells you where the logs are.
set -euo pipefail

STEP="${1:-step2}"; shift || true

# ---- defaults, all overridable from the environment -------------------------------
: "${PROFILE:=local,fast}"                 # slurm,precise when submitting instead
: "${WORKDIR:=work/}"                      # /no_backup/... for real runs
: "${SPECIES_TREE:=data/species_tree.newick}"
: "${FAMILY_INFO:=genefam.csv}"
: "${NXF_OPTS:=-Xms500M -Xmx2G}"
export NXF_OPTS

cd "$(dirname "$0")"
mkdir -p reports logs

command -v module >/dev/null 2>&1 && module load Java || echo "note: no 'module' command, assuming java is on PATH"
command -v nextflow >/dev/null || { echo "ERROR: nextflow not on PATH"; exit 1; }
command -v java     >/dev/null || { echo "ERROR: java not on PATH (module load Java)"; exit 1; }

# ---- PRE-FLIGHT: the checks whose absence cost four days ---------------------------
echo "== pre-flight"
[[ -f "$STEP.nf" ]]        || { echo "ERROR: $STEP.nf not found"; exit 1; }
[[ -f "$FAMILY_INFO" ]]    || { echo "ERROR: family info '$FAMILY_INFO' not found"; exit 1; }
[[ -f "$SPECIES_TREE" ]]   || { echo "ERROR: species tree '$SPECIES_TREE' not found"; exit 1; }

# GeneRax needs a STRICTLY BINARY species tree. step2.nf now refuses to start otherwise,
# but check here too so the answer arrives in one second rather than after staging.
if printf '%s' "$*" | grep -q -- '--run_generax'; then
  python3 - "$SPECIES_TREE" <<'PY'
import sys
nwk = open(sys.argv[1]).read().strip()
kids, bad = [], 0
for ch in nwk:
    if   ch == '(': kids.append(1)
    elif ch == ',' and kids: kids[-1] += 1
    elif ch == ')':
        n = kids.pop() if kids else 0
        if n > 2: bad += 1
print(f"   species tree {sys.argv[1]}: {nwk.count('(')} internal nodes, {bad} multifurcating")
if bad:
    sys.exit("ERROR: GeneRax needs a strictly binary species tree. Resolve it:\n"
             "  python workflow/check_tree.py <in.newick> <ids> <out.newick> --random-resolve --seed 1\n"
             "  then re-run with SPECIES_TREE=<out.newick>")
PY
  command -v generax >/dev/null && echo "   generax: $(command -v generax)" \
      || echo "   WARNING: generax not on PATH -- --run_generax will fail (the bioconda build is non-MPI; load an MPI build)"
  command -v mpirun  >/dev/null && echo "   mpirun:  $(command -v mpirun)"  || echo "   WARNING: mpirun not on PATH"
fi
echo "   profile=$PROFILE  workdir=$WORKDIR  species_tree=$SPECIES_TREE"

# ---- run --------------------------------------------------------------------------
TS=$(date +%Y%m%d-%H%M%S)
LOG="logs/${STEP}.${TS}.log"
echo "== running $STEP (log: $LOG)"
set -x
nextflow run "$STEP.nf" \
    -profile "$PROFILE" \
    -resume -w "$WORKDIR" \
    -ansi-log false \
    --family_info  "$FAMILY_INFO" \
    --species_tree "$SPECIES_TREE" \
    -with-report   "reports/report.${STEP}.html" \
    -with-trace    "reports/trace.${STEP}.txt" \
    "$@" 2>&1 | tee "$LOG"
set +x

# ---- post-mortem help -------------------------------------------------------------
echo "== done. failed tasks (if any), with the REAL error:"
awk '/terminated with an error/ {print}' "$LOG" | head -5
echo "   read a failing task's stderr with:  cat $WORKDIR/<hash>/.command.err"
echo "   (the errorStrategy label in $STEP.nf is a GUESS -- trust .command.err and the tool log)"
