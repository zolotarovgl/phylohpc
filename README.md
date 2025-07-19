# Phylogeny pipeline   

The following command will run the phylogeny pipeine for the families defined in `genefam.csv`:  

```bash 
# prepare input files:
bash prepare_fasta.sh species_list data/input.fasta
```


```bash
bash pipeline.sh   
```


# TODOs:   

- add GeneRax 
- add GeneRax explainer   
- add updated species tree and an example one 

Since GeneRax takes so long, it should be launched after the phylogenies have been done
