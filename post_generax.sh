source workflow/functions.sh
read_config configs/config.txt
# run possvm
comm <(ls results/generax/*groups.csv | xargs -n1 basename | cut -f 1-3 -d . | sort) <(ls results/generax/*treefile | xargs -n1 basename | cut -f 1-3 -d . | sort) -31 > tmp/torun.2
N=$(wc -l tmp/torun.2 | awk '{print $1}')
echo "$N jobs to run"

for PREF in $(cat tmp/torun.2); do
    python phylogeny/main.py possvm -t results/generax/${PREF}.treefile -r $REFNAMES -s $REFSPECIES -o ${PREF}"."

done
