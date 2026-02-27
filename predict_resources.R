#############################################
# Predict resources 
#############################################
library(stringr)
library(dplyr)

ids_fn = 'ids.txt'
cluster_dir = 'results/clusters'
models_rds = 'workflow/models/models.rds'
outfile = 'resources.tsv'

# default values per job
defaults = setNames(c(500,200,300,5,30,1),paste0(rep(c('aln','phy','pvm'),2),rep(c('_mem','_time'),each = 3)))
min_mem = 100
min_time = 1
max_mem = 10000
max_time = 48*60
increase = 0.1 # to increase predicted resources by 10 percent
# the minimum values for the jobs 
###############################################
ids = readLines(ids_fn)
f = list.files(cluster_dir,full = TRUE,pattern = 'fasta')
names(f) = gsub('.fasta','',basename(f))
f = f[intersect(ids,names(f))]
message(sprintf('%s ids. %s .fastas found',length(ids),length(f)))
#############################################
# Predict
# the input data should contain mlen and nseq
#############################################
predict_res = function(model_m,model_t,input){
	pred_m = exp(predict(model_m, newdata=input))
	pred_t = exp(predict(model_t, newdata=input))
	data.frame(mem = pred_m,time = pred_t)
}
.get_nseq = function(x) setNames(length(grep('>',readLines(x))),'nseq')
.get_mlen = function(x) setNames(median(nchar(readLines(x)[!grepl('>',readLines(x))])),'mlen')
input = t(sapply(f,FUN = function(x) c(.get_nseq(x),.get_mlen(x))))
input = as.data.frame(input)
models = readRDS(models_rds)
o = lapply(names(models),FUN = function(job_name){
	model_m = models[[job_name]][['mem']]
	model_t = models[[job_name]][['time']]
	u = predict_res(model_m,model_t,input)
	u$mem[u$mem<min_mem] = min_mem
	u$time[u$time<min_time] = min_time
	u
})
names(o) = names(models)
names(o) = tolower(names(o))
for(i in seq_along(o)){
	colnames(o[[i]]) = paste0(names(o)[i],'_',colnames(o[[i]])) 
}
o = do.call(cbind,setNames(o,NULL))
#############################################
# Round and convert to values usable by the nextflow: 
#############################################
pred = o
for(i in seq_along(pred)){
	pred[,i] = ceiling(pred[,i]*(1+increase))
}
for(x in names(defaults)){
	if(!x %in% colnames(pred)){
		pred[,x] = defaults[x]
	}
}
timcols = grep('_time',colnames(pred),value = t)
memcols = grep('_mem',colnames(pred),value = TRUE)
pred = as.data.frame(apply(pred,2,ceiling))
o = t(apply(pred,2,range))
b = paste0(sapply(1:nrow(o),FUN = function(i) sprintf('%s: %s - %s,%s',rownames(o)[i],o[i,1],o[i,2],ifelse(grepl('time',rownames(o)[i]),'min','MB'))),collapse = '\n')
message(b)
b1 = t(apply(pred[,memcols]>=max_mem,1,FUN = function(x) c(sum(x),length(x))))
b2 = t(apply(pred[,timcols]>=max_time,1,FUN = function(x) c(sum(x),length(x))))
if(colSums(b1)[1]>0){
	message(sprintf('%s / %s jobs with MEM >= max_mem (%s)',colSums(b1)[1],sum(b1),max_mem))
}
if(colSums(b2)[1]>0){
	message(sprintf('%s / %s jobs with TIME>= max_time (%s)',colSums(b2)[1],sum(b2),max_time))
}

for(i in memcols){
	pred[,i][pred[,i]>=max_mem] = max_mem
}
for(i in timcols){
	pred[,i][pred[,i]>=max_time] = max_time
}
round_base = function(x,base = 60){ifelse(x>=base,ceiling(x/base)*base,x)}
# round to hours for long jobs
for(i in timcols){
	pred[,i] = round_base(pred[,i],base = 60)
}
colSums(pred[,timcols]>=1024)
# round up for memory - GB
for(i in memcols){
	pred[,i] = round_base(pred[,i],base = 1024)
}
colSums(pred[,timcols]/60>=12)
colSums(pred[,memcols]>=1024)

# now, round and convert 
convert_mem = function(x) paste0(ceiling(x/100)*100,'.MB')
convert_time = function(x,min_time = 5){paste0(ceiling(pmax(x,min_time)),'.min')}

for(i in grep('_mem',colnames(pred))){
	pred[,i] = convert_mem(pred[,i])
}
for(i in grep('_time',colnames(pred))){
	pred[,i] = convert_time(pred[,i])
}
pred = cbind(data.frame(id = rownames(pred)),pred)
rownames(pred) = NULL
write.table(pred,outfile,sep = '\t',quote = FALSE,row.names = FALSE)
message(sprintf("Created: %s",outfile))
