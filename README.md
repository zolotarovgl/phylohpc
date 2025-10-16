
# Installation  

```bash
git clone --recurse-submodules https://github.com/zolotarovgl/phylohpc.git
cd phylohpc
mamba env create -f workflow/environment.yaml
```




# Configuration  

`configs/config.txt` - a config file with pipeline settings.   
`species_list` - the list of species prefixes to use for the pipeline  

# Phylogeny pipeline   


Pull species' proteomes from Xavi's database. If you don't know who Xavi is, then you'd have to prepare a joint protomes fasta yourself: 
```bash 
# prepare input files:
bash workflow/prepare_fasta.sh species_list data/input.fasta
```


The following command will launch the pipeline (first pass). 
```bash
bash pipeline.sh   
```
Some jobs will inevitably fail do to the amount of allocated resources.
The resources are specified in the `configs/config.txt`. For instance:  

```
MEM_S1=1G
MEM_S2=1G
TIME_S1=00:10:00
TIME_S2=00:30:00
TIME_S3=00:01:00
TIME_S3_SHORT=1:00:00
TIME_S3_LONG=12:00:00
MEM_SHORT=1G
MEM_LONG=1G
MEM_S5=10G
TIME_S5=1:00:00
```

    - `S1` - search step   
    - `S2` - clustering step   
    - `S3` - alignment and phylogeny. The job itself is a job that submits 2 job arrays - long and short. The long job array is intended for big gene families. The difference between the big and short family is determined by the following line in the config:  

```bash
# separate between small groups
MAX_SEQ=300
``` 
    - `S5` - GeneRax 



Thus, the status of each homology group should be checked:
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
#./workflow/get_hg_status.sh > hg_status.tab
for ID in $(cat ids); do ./workflow/check_status.sh $ID; done | grep -v '#' | awk '{print $1"\t"$2$3$4$5}' > hg_status.tab
cat hg_status.tab  | grep 1000 | cut -f 1  > tmp/torun


# submit alignment jobs:
TIME=00:30:00
MEM=500M
NCPU=4
while IFS="." read -r PREF FAMILY HG; do
    python submit_hg.py configs/config.txt --time $TIME --cpus $NCPU --mem $MEM --mode align --pref $PREF --family $FAMILY --hg $HG --mafft ""  --json aln.info.json
done < tmp/torun

./workflow/get_hg_status.sh > hg_status.tab
cat hg_status.tab | grep 1100 | cut -f 1  | grep -w -f ids > tmp/torun

# remove gap-only sequences
for ID in $(cat tmp/torun); do
echo $ID
clean_fasta results/align/${ID}.aln.fasta tmp_fa 20; mv tmp_fa results/align/${ID}.aln.fasta
done

TIME=1-00:00:00
MEM=500M
NCPU=6
while IFS="." read -r PREF FAMILY HG; do
    python submit_hg.py configs/config.txt --time $TIME --cpus $NCPU --mem $MEM --mode phylogeny --pref $PREF --family $FAMILY --hg $HG  --json phy.info.json
done < tmp/torun
```

## Check phylogeny jobs status:

```bash
for ID in $(cat phy.info.json  | grep jobid | cut -f 2 -d :  | sed -E 's/ |"|,//g'); do sacct -j $ID | awk 'NR==3' ; done
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


# Gather annotations per species

```bash
ID=Clacla
mkdir -p annotations/${ID}/
bash workflow/gather_anno.sh --id $ID > annotations/$ID/$ID.all.tab
cat annotations/Clacla/Clacla.all.tab  | cut -f 2 | cut -f 1 -d . | sort | uniq -c  | sort -rn


F=annotations/$ID/$ID.all.tab
cat $F | awk -F'\t' 'BEGIN{print "#Prefix\tTotal\tClassified\tPerc_Classified"}{split($2,a,".");PREF=a[1];counter[PREF]+=1;if($2!~/Unclass/){class[PREF]+=1}}END{for(k in counter){print k"\t"counter[k]"\t"class[k]"\t"class[k]/counter[k]}}' 

for PREF in $(cat $F | cut -f 2 | cut -f 1 -d . | sort | uniq -c  | sort -rn  | awk '$1>=5 {print $2}'); do 
echo $PREF
cat $F | awk -v PREF=$PREF '{split($2,a,".")}{if(a[1]==PREF){print $0}}' > annotations/${ID}/${ID}.${PREF}.tab
done

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
