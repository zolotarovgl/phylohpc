library(stringr)
open = function(x){system(sprintf('open %s',x))}
# nextflow trace file parsing 
.trace_time_convert=function(x){
  sapply(x,function(z){
    z=gsub(" ","",z)
    m=regmatches(
      z,
      regexec("^(?:(\\d+(?:\\.\\d+)?)m)?(?:(\\d+(?:\\.\\d+)?)s)?(?:(\\d+(?:\\.\\d+)?)ms)?$",z)
    )[[1]]
    
    if(length(m)==0) return(NA_real_)
    
    mins=ifelse(m[2]=="",0,as.numeric(m[2]))
    secs=ifelse(m[3]=="",0,as.numeric(m[3]))
    # milliseconds captured but ignored (treated as 0)
    
    mins*60+secs
  })
}
.trace_mem_convert <- function(x) {
  sapply(x, function(z) {
    z <- gsub(" ", "", z)
    
    m <- regmatches(
      z,
      regexec("^(\\d+(?:\\.\\d+)?)(MB|GB)$", z)
    )[[1]]
    
    if (length(m) == 0) return(NA_real_)
    
    value <- as.numeric(m[2])
    unit  <- m[3]
    
    if (unit == "GB") {
      value * 1024
    } else {
      value
    }
  })
}
# nextflow values parsing 

.nf_convert_time = function(x){
	tmap = c('s' = 1,'min' = 60,'h' = 60^2)
	as.numeric(str_split(x,'\\.',simplify = T)[,1])*tmap[str_split(x,'\\.',simplify = T)[,2]]
}
.nf_convert_mem = function(x){
	as.integer(str_split(d$mem_res,'\\.',simplify = T)[,1])
}

# GeneRax parsing 
parse_generax_log = function(fn){
	o = readLines(fn)
	ncpu = as.integer(str_extract(o[10],'[0-9]+'))
	n_genes = as.integer(unlist(str_split(o[33],' '))[length(unlist(str_split(o[33],' ')))])
	n_species = as.integer(unlist(str_split(o[32],' '))[length(unlist(str_split(o[32],' ')))])	
	b = o[grepl('radius',o)&grepl('Optimizing',o)]
	id = gsub('_generax','',unlist(str_split(o[7],' '))[4])
	out = data.frame(id = id,ncpu = ncpu,n_genes = n_genes,n_species = n_species,time = gsub('\\[|\\]','',str_split(b,' ',simplify = T)[,1]),radius = as.integer(gsub('radius=|\\.\\.\\.','',str_split(b,' ',simplify = T)[,9])))
	out$time = as.numeric(as.difftime(out$time, format = "%H:%M:%S",units = 'secs'))
	# Likelihoods
	ll = str_split(gsub('\t','',o[grepl('JointLL',o)]),' ',simplify = T)
	out$JointLL = as.numeric(gsub('JointLL=','',ll[,1]))
	out$RecLL = as.numeric(gsub('RecLL=','',ll[,2]))
	out$LibpllLL = as.numeric(gsub('LibpllLL=','',ll[,3]))
	
	# Total execution time
	out$time_total = as.numeric(as.difftime(gsub('\\[|\\]','',str_split(o[length(o)],' ',simplify = T)[,1]),format = "%H:%M:%S",units = 'secs'))	
	out
}

# count the sequnces 
count_seq = function(x){
	o = as.integer(str_split(system(sprintf('grep -c ">" %s',paste0(x,collapse = ' ')),intern = TRUE),':',simplify = T)[,2])
	if(!is.null(names(x))){
		names(o) = names(x)
	}else{
		names(o) = basename(x)
	}
	o
}
settheme = function(){theme_set(theme_cowplot())}
