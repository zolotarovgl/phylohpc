#!/bin/bash
set -e 
# Ensure all arguments are provided
if [ "$#" -ne 4 ]; then
    echo "Error: You must provide exactly 4 arguments: FAMILY, GENEFAM_INFO, INFASTA, OUTDIR"
    echo "Usage: $0 <FAMILY> <GENEFAM_INFO> <INFASTA> <OUTDIR>"
    exit 1
fi

FAMILY=$1 # Family from genefam file, e.g., Insulin
GENEFAM_INFO=$2
INFASTA=$3
OUTDIR=$4

#
echo "$(date '+%Y-%m-%d %H:%M:%S') $0 start"
source workflow/functions.sh
CONFIG=configs/config.txt
read_config $CONFIG

NCPU=10
# Run the Python script
mkdir -p $OUTDIR
echo "python phylogeny/main.py hmmsearch -f $INFASTA -g $GENEFAM_INFO $FAMILY -o $OUTDIR --pfam_db $PFAMDB --domain_expand 30"
python phylogeny/main.py hmmsearch -f $INFASTA -g $GENEFAM_INFO $FAMILY -o $OUTDIR --pfam_db $PFAMDB --domain_expand $DOMAIN_EXPAND --ncpu $NCPU
echo "$(date '+%Y-%m-%d %H:%M:%S') $0 done"
exit 0
