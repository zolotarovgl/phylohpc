
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

- don't run for the HGs without the species of interest!  
- add GeneRax 
- add GeneRax explainer   
- add updated species tree and an example one



# Changes   


 
