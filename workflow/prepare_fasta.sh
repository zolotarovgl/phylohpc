#!/bin/bash

# Ensure the required arguments are provided
if [ "$#" -lt 2 ] || [ "$#" -gt 3 ]; then
    echo "Usage: $0 <SPECIES_PREFS> <OUTPUT_FILE> [DB_DIR]"
    exit 1
fi

PREFS=$1
OUTFILE=$2
DB=${3:-/users/asebe/xgraubove/genomes/data/}

TEMPDIR=tmp_prep

if [ -d $TEMPDIR ]; then
    echo "ERROR: $TEMPDIR already exists."
    exit 1
fi

mkdir -p $TEMPDIR
for SP in $(cat $PREFS); do 
    echo $SP
    F=$DB/${SP}_long.pep.fasta
    if [ -f $F ]; then
        cp $F $TEMPDIR
    else
        echo "ERROR: ${F} doesn't exist. Skipping ..."
        exit 1
    fi
done

# Concatenate the results
cat $(ls $TEMPDIR/*fasta | grep -f $PREFS) > $OUTFILE
rm -rf $TEMPDIR
samtools faidx $OUTFILE
echo "Created: ${OUTFILE}"
exit 0
