# Gene family information 

The full family information is stored in the `data/gene_families_searchinfo.csv`.   
This file contains the following fields:   

	 1. Family name   
	 2. Comma-separated list of PFAM domain names fetcheable with `hmmfetch` from Pfam A database.   
	 3. Inflation parameter for MCL clustering   
	 4. Legacy: Minimum number of sequences
	 5. Legacy: `GA`  
	 6. Annotation name - e.g. 'neural'  
	 7. Short prefix (e.g. 'neu') to be used in the file names.  

## Family descriptions   

* `tfs` - transcription factors
* `neu` - gene families relevant for neuronal proteins (legacy from placozoa paper)  
* `chr` - chromatin    
* `rbp` - RNA-binding proteins  
* `sig` - signalling  
* `adh` - adhesion  
* `ion` - ion channels  
* `mus` - muscle proteins  
* `cil` - cilia proteins  
* `npe` - neuropeptides  
* `cal` - Ca2+ 
* `xxx` - miscellaneous  
