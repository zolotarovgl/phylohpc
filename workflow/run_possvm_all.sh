#!/bin/sh
source workflow/functions.sh
read_config configs/config.txt
# logs with finished trees 
grep -l 'Analysis results written to' ${TREE_DIR}/*.log | xargs -n1 basename |  sed 's/.log//g'  | sort | uniq  > tmp/tree_done
# logs with all trees
find ${TREE_DIR}/ -name '*log' | xargs -n1 basename |  sed 's/.log//g'  | sort | uniq > tmp/tree_all

./workflow/get_hg_status.sh  > hg_status.tab
comm <(cat hg_status.tab | grep 1110 | cut -f 1 | sort | uniq) tmp/tree_done -12 > tmp/torun

N=$(wc -l tmp/torun | awk '{print $1}')
echo "Running POSSVM for ${N} HGs..."

for ID in $(cat tmp/torun); do 
PREF=${ID%%.*}
FAMILY=${ID#*.};FAMILY=${FAMILY%%.*}
HG=${ID##*.}
echo $PREF $FAMILY $HG
TREE_FILE=${TREE_DIR}/${PREF}.${FAMILY}.${HG}.treefile
python phylogeny/main.py possvm -t $TREE_FILE --refsps $REFSPECIES -r $REFNAMES -o ${PREF}.${FAMILY}.${HG}"."
done 
