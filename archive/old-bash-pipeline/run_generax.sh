#!/bin/bash

# Function to display usage instructions
usage() {
  echo "Usage: $0 --alignment <alignment.fasta> --genetree <gene_tree.newick> --iqtreefile <iqtree.file> --speciestree <species_tree.newick> [--outdir <output_directory>] [--tmpdir <temp_directory>] [--prefix <prefix>] [--maxspr <max_spr_radius>] [--ncpu <number of cores>]"
  exit 1
}

# Default optional arguments
OUTDIR="results_annotation/generax"
GENERAX_TMP="tmp/generax_configs"
PREF="Family1"
MAXSPR=5
NCPU=1
# Parse named arguments
while [[ "$#" -gt 0 ]]; do
  case $1 in
    --alignment) ALIGNMENT="$2"; shift ;;
    --genetree) GENETREE="$2"; shift ;;
    --iqtreefile) IQTREEFILE="$2"; shift ;;
    --speciestree) SPECIESTREE="$2"; shift ;;
    --outdir) OUTDIR="$2"; shift ;;
    --tmpdir) GENERAX_TMP="$2"; shift ;;
    --prefix) PREF="$2"; shift ;;
    --maxspr) MAXSPR="$2"; shift ;;
    --ncpu) NCPU="$2"; shift ;;
    *) echo "Unknown parameter passed: $1"; usage ;;
  esac
  shift
done

# Check if required arguments are provided
if [[ -z $"PREFIX" || -z "$ALIGNMENT" || -z "$GENETREE" || -z "$IQTREEFILE" || -z "$SPECIESTREE" ]]; then
  echo "Error: --prefix, --alignment, --genetree, --iqtreefile, and --speciestree are required."
  usage
fi

# Check if required files exist
if [[ ! -f "$ALIGNMENT" ]]; then
  echo "Error: Alignment file '$ALIGNMENT' does not exist."
  exit 1
fi
if [[ ! -f "$GENETREE" ]]; then
  echo "Error: Gene tree file '$GENETREE' does not exist."
  exit 1
fi
if [[ ! -f "$IQTREEFILE" ]]; then
  echo "Error: IQ-TREE file '$IQTREEFILE' does not exist."
  exit 1
fi
if [[ ! -f "$SPECIESTREE" ]]; then
  echo "Error: Species tree file '$SPECIESTREE' does not exist."
  exit 1
fi

# Create output and temporary directories
mkdir -p "${OUTDIR}"
mkdir -p "${GENERAX_TMP}"

# Generate the family file
MODEL=$(grep 'Model of substitution:' "$IQTREEFILE" | awk '{print $NF}')
FAMFILE="${GENERAX_TMP}/${PREF}.famfile"
echo "[FAMILIES]" > "$FAMFILE"
echo "- ${PREF}" >> "$FAMFILE"
echo "alignment = ${ALIGNMENT}" >> "$FAMFILE"
echo "starting_gene_tree = ${GENETREE}" >> "$FAMFILE"
echo "subst_model = ${MODEL}" >> "$FAMFILE"
echo "Config file: ${FAMFILE} has been generated."

# Run GeneRax
RUN_OUTDIR="${OUTDIR}/${PREF}"
mkdir -p "$RUN_OUTDIR"
mpirun --use-hwthread-cpus --oversubscribe -n $NCPU  generax -s "$SPECIESTREE" -f "$FAMFILE" --per-family-rates -r UndatedDL -p "$RUN_OUTDIR" --max-spr-radius "$MAXSPR" --strategy SPR

