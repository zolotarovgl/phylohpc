source workflow/functions.sh
read_config configs/config.txt
# run possvm
for F in $(find results/generax/ -name '*.treefile' | sort); do             
    PREF=$(basename $F  | cut -f 1-3 -d .)
if [ ! -f "$F.ortholog_groups.csv" ]; then
    echo $PREF
    python phylogeny/main.py possvm -t $F -r $REFNAMES -s $REFSPECIES -o ${PREF}"."
fi
done

