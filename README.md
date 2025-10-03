
# Installation  

```bash
git clone --recurse-submodules https://github.com/zolotarovgl/phylohpc.git
```


# Configuration  

`configs/config.txt` - a config file with pipeline settings.   
`species_list` - the list of species prefixes to use for the pipeline  

# Phylogeny pipeline   

The following command will run the phylogeny pipeine for the families defined in `genefam.csv`:  

```bash 
# prepare input files:
bash workflow/prepare_fasta.sh species_list data/input.fasta
```


The following command will launch the pipeline (first pass). 
```bash
bash pipeline.sh   
```
Some jobs will inevitably fail do to the amount of allocated resources. Thus, the status of each homology group should be checked:
```bash
bash workflow/get_hg_status.sh
``` 
`[HG]\t[clustering][alignment][phylogeny][possvm]`

# TODOs:   

- make sure the phylogenies are rerun (updated) not overwritten  
- better annotation gathering with the ortholog support etc
- inflation parameter should be specified in the genefam file
- don't run for the HGs without the species of interest!  
- add GeneRax 
- add GeneRax explainer   
- add updated species tree and an example one
- ~~make sure to NOT run possvm on the phylogenies that did not finish - relaunch instead~~ 
- ~~gather anno should account for the families without clusterings~~



# Changes   

The following HMM models no longer seem to be a part of the PFAM-A: removed
```
dsrm_Ferlin
GKAPp
KAT11
NCam-PTHR12231
RIMSbp-PTHR14234
Synaptotagmin-PTHR10024
```


# Useful commands 

## Family status 

```bash
source ./workflow/functions.sh
check_all_families genefam.csv

```

```bash
./workflow/get_hg_status.sh  > hg_status.tab
comm <(cat hg_status.tab | grep 1110 | cut -f 1 | sort | uniq) tree_done -12 > tmp/torun

read_config configs/config.txt
for ID in $(cat tmp/torun); do 
PREF=${ID%%.*}
FAMILY=${ID#*.};FAMILY=${FAMILY%%.*}
HG=${ID##*.}
echo $PREF $FAMILY $HG
TREE_FILE=${TREE_DIR}/${PREF}.${FAMILY}.${HG}.treefile
python phylogeny/main.py possvm -t $TREE_FILE --refsps $REFSPECIES -r $REFNAMES -o ${PREF}.${FAMILY}.${HG}"."
done 
```
 
