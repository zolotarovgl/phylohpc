
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


```bash
bash pipeline.sh   
```

# TODOs:   

- add GeneRax 
- add GeneRax explainer   
- add updated species tree and an example one



# Changes   


 
