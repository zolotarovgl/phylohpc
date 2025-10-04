
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
- inflation parameter should be specified in the genefam file
- don't run for the HGs without the species of interest!  
- add GeneRax 
- add GeneRax explainer   
- add updated species tree and an example one
- ~~make sure to NOT run possvm on the phylogenies that did not finish - relaunch instead~~ 
- ~~gather anno should account for the families without clusterings~~
- ~~better annotation gathering with the ortholog support etc~~



# Useful commands 


## Clever job submission   

`submit_hg.py` now has an argument `--json` - if exists, this file will provide the info on the last job launched for a particular `PREF.FAMILY.HG`. 
This allows to increase the resources, if the job has timed out or OOM:  `--mem_increase`, `--time_increase`  


For example, to submit alignment and phylogeny jobs:

```bash
# run only for HGs with Clacla
grep -l '>Clacla' results/clusters/*fasta | xargs -n1 basename | sed 's/.fasta//g' | sort | uniq > ids
./workflow/get_hg_status.sh > hg_status.tab
cat hg_status.tab | grep -w -f ids  | grep 1000 | cut -f 1  > tmp/torun

# submit alignment jobs:
TIME=00:30:00
MEM=500M
NCPU=4
JSON=aln.info.json
while IFS="." read -r PREF FAMILY HG; do
    python submit_hg.py configs/config.txt --time $TIME --cpus $NCPU --mem $MEM --mode align --pref $PREF --family $FAMILY --hg $HG --mafft ""  --json $JSON 
done < tmp/torun

./workflow/get_hg_status.sh > hg_status.tab
cat hg_status.tab | grep 1100 | cut -f 1  | grep -w -f ids > tmp/torun

# remove gap-only sequences
for ID in $(cat tmp/torun); do
echo $ID
clean_fasta results/align/${ID}.aln.fasta tmp_fa 20; mv tmp_fa results/align/${ID}.aln.fasta
done

TIME=01:00:00
MEM=500M
NCPU=4
JSON=phy.info.json
while IFS="." read -r PREF FAMILY HG; do
    python submit_hg.py configs/config.txt --time $TIME --cpus $NCPU --mem $MEM --mode phylogeny --pref $PREF --family $FAMILY --hg $HG  --json $JSON
done < tmp/torun
```

## Run possvm for the hgs wth finished phylogenies 

```bash
sbatch --output=output.log --error=error.log --mem=1G --time=01:00:00 workflow/run_possvm_all.sh 
```


## Family status 

Check the pipeline status of the families defined in `genefam.csv`:  


```bash
source ./workflow/functions.sh
check_all_families genefam.csv
```

Submit clustering jobs for the families with missing clusterings:

```bash
bash workflow/submit_family.sh configs/config.txt sig.GPCRrhod --mem_s1 500M --mem_s2 30G  --dry
```  





 
# Changes   

The following HMM models no longer seem to be a part of the PFAM-A: removed
```
# missing HMM files:
for ID in $(awk '{print $2}' genefam.csv  | sed 's/,/\n/g' | sort | uniq ); do
hmmfetch $PFAMDB $ID > out.hmm 2> err.log
if [ $? -ne 0 ]; then
    echo "Error fetching $ID"
fi
done
```
