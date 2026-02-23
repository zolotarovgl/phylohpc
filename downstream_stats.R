# Summarize the job running stats 
# Number of sequences per job 
library(ggplot2)
library(data.table)
library(stringr)
options(max.print = 30)
#WORKDIR=/no_backup/asebe/gzolotarov/work/
#HASH=dc/dda4f9
#FASTA=$(ls ${WORKDIR}/${HASH}*/*fasta)
#NSEQ=$(grep -c '>' $FASTA)

# gather the number of sequences 
get_job_workdir = function(workdir,hash){
	hash_l = unlist(str_split(hash,'/'))
	o = list.files(paste0(workdir,'/',hash_l[1]),pattern = hash_l[2])
	return(paste0(workdir,'/',hash_l[1],'/',o))
}
.count_fasta = function(workdir,hash){
	f = list.files(get_job_workdir(workdir,hash),pattern = 'fasta',full = TRUE)
	o = sapply(f,FUN = function(x){
		system(sprintf('grep -c ">" %s',x),intern = TRUE)
	})
	max(as.integer(o))
}
read_command_err = function(workdir,hash){
	readLines(list.files(get_job_workdir(workdir,hash),all.files = TRUE,pattern = 'command.err',full = TRUE))
}
.get_slurm_stats = function(ids){

	cmd = sprintf(
		"sacct -j %s --format=JobIDRaw,ReqMem,Timelimit,Submit --noheader -P",
		paste(ids, collapse=",")
	)

	res = read.table(
		text = system(cmd, intern=TRUE),
		sep="|",
		stringsAsFactors=FALSE,
		col.names=c("native_id","ReqMem","Timelimit","Submit")
	)

	res = res[res$ReqMem != "",]

	res$ReqMem_MB = ifelse(
		grepl("G$", res$ReqMem),
		as.numeric(sub("G","",res$ReqMem)) * 1024,
		as.numeric(sub("M","",res$ReqMem))
	)

	res$Timelimit_min =
		as.numeric(substr(res$Timelimit,1,2))*60 +
		as.numeric(substr(res$Timelimit,4,5)) +
		as.numeric(substr(res$Timelimit,7,8))/60

	res$Submit_dt = as.POSIXct(res$Submit, format="%Y-%m-%dT%H:%M:%S", tz="UTC")

	return(res)
}

#############################
# Sequence counts
#############################
f = list.files('results/clusters/',pattern = 'fasta',full = TRUE)
names(f) = gsub('.fasta','',basename(f))
f = f[names(f) %in% d$hg_id]
nseq = sapply(f,FUN = function(x){
	system(sprintf('grep -c ">" %s',x),intern = TRUE)
})
names(nseq) = gsub(".fasta","",basename(names(nseq)))


# Seqeunce lengths
f = list.files('results/clusters/',pattern = 'fasta',full = TRUE)
names(f) = gsub('.fasta','',basename(f))
f = f[names(f) %in% d$hg_id]
lens = lapply(f,FUN = function(x) nchar(readLines(x)[!grepl('>',readLines(x))]))
len = sapply(lens,median)
#############################
# Work outputs
#############################
workdir = '/no_backup/asebe/gzolotarov/work/'
d = read.table('trace.txt',sep = '\t',fill = TRUE,header = TRUE)
rownames(d) = d$hash
d$submit_dt = as.POSIXct(d$submit, format="%Y-%m-%d %H:%M:%OS", tz="UTC")
span = difftime(max(d$submit_dt), min(d$submit_dt), units="mins")
span_sec = as.numeric(difftime(max(d$submit_dt), min(d$submit_dt), units="secs"))
days    = span_sec %/% 86400
hours   = (span_sec %% 86400) %/% 3600
minutes = (span_sec %% 3600) %/% 60
message(sprintf("%s - %s\n%d days %d hours %d minutes",min(d$submit_dt),max(d$submit_dt), days, hours, minutes))
d$job_name = str_split(d$name,' ',simplify = T)[,1]
d$peak_rss_mb = as.numeric(sub(" MB","", d$peak_rss))
d$hg_id = gsub('\\(|\\)','',str_split(d$name,' ',simplify = T)[,2])
d$minutes = sapply(d$realtime, function(x){
	x = trimws(x)
	m = if(grepl("m", x)) as.numeric(sub("m.*", "", x)) else 0
	s = if(grepl("s", x)) as.numeric(sub(".*m\\s*", "", sub("s", "", x))) else as.numeric(sub("s", "", x))
	if(is.na(s)) s = 0
	m + s/60
})
d$n = as.integer(nseq[d$hg_id])
d$median_len = as.numeric(len[d$hg_id])
d = d[d$exit == 0,]

# Slurm stats
res = .get_slurm_stats(ids = unique(d$native_id))
d = merge(d, res, all.x=TRUE)

# fraction of requensted memory used 
d$mem_perc = d$peak_rss_mb/d$ReqMem_MB
d$time_perc = d$minutes/d$Timelimit_min
saveRDS(d,'results/downstream/job_info.rds')
d = readRDS('results/downstream/job_info.rds')



par(mfrow = c(1,2))
plotname = 'test.pdf'
pdf(plotname,height = 5,width = 10)
par(mfrow = c(1,2))
boxplot(split(d$mem_perc,d$job_name),ylab = 'mem perc',outline = FALSE)
boxplot(split(d$time_perc,d$job_name),ylab = 'time perc',outline = FALSE)
dev.off()
##########################################
options(tibble.width = Inf)
status = "COMPLETED"
status = c("COMPLETED","CASHED")
table(d$status)
library(dplyr)
table(is.na(d$peak_rss_mb))
d[,]%>%filter(!is.na(peak_rss_mb))%>%group_by(job_name)%>%summarize(
	n = n(),
	mem_requ = sum(ReqMem_MB)/1024,
	mem_used = sum(peak_rss_mb)/1024,
	time_requ = sum(Timelimit_min),
	time_used = sum(minutes),
	mem_perc = mem_used/mem_requ,
	time_perc = time_used/time_requ,
	mem_waste = mem_requ * (1-mem_perc),
	time_waste = time_requ * (1-time_perc))%>%arrange(-mem_waste)
summary(d[dd$job_name == 'ALN',]$peak_rss_mb)
summary(d[d$job_name == 'ALN',]$ReqMem_MB)
summary(d[d$job_name == 'PHY',]$minutes)
# average of 3 minutes? is this real?
d[d$job_name == 'PHY',]
hash = '20/b7b09c'
get_job_workdir(workdir,hash)
table(d$status)
# huge time waste for the phylogenies 
# Median is 1Gb but it only uses around 300 
# ALN job requests too much memory 	
# I overask for the time clearly 
####################################
# Resource scaling 
####################################
plotname = 'test.pdf'
pdf(plotname,height = 5,width = 10)
par(mfrow = c(1,2))
for(job_name in c("ALN","PHY")){
	p = d[d$job_name == job_name,]
	plot(p$n,p$minutes,pch = 16,log = 'xy',xlab = '# sequences', ylab = 'Duration, min',font.main = 1,main = job_name)
	plot(p$n,p$peak_rss_mb,pch = 16,log = 'xy',xlab = '# sequences', ylab = 'Memory, MB',font.main = 1,main = job_name)
}
dev.off()
options(max.print = 100)
dir.create('results/downstream')
saveRDS(p,'results/downstream/job_info.rds')
####################################
# Predict runtimes
####################################

d = readRDS('results/downstream/job_info.rds')
table(d$job_name)
dd = d[,c('hg_id','job_name','hash','minutes','Timelimit_min','peak_rss_mb','ReqMem_MB')]
colnames(dd) = c('hg_id','job_name','hash','time_used','time_requ','mem_used','mem_requ')
dd$mem_perc = dd$mem_used/dd$mem_requ
dd$time_perc = dd$time_used/dd$time_requ
dd$nseq = as.integer(nseq[dd$hg_id])
dd$mlen = len[dd$hg_id]


models = list()
for(id in c("ALN","PHY")){
	m1 = lm(log(mem_used) ~ log(nseq) * log(mlen), data=dd[dd$job_name == id,])
	m2 = lm(log(time_used) ~ log(nseq) * log(mlen), data=dd[dd$job_name == id,])
	models[[id]] = list(mem = m1,time = m2)
}


predict_res = function(model_m,model_t,new){
	pred_m = exp(predict(model_m, newdata=new))
	pred_t = exp(predict(model_t, newdata=new))
	data.frame(mem = pred_m,time = pred_t)
}

job_name = 'ALN'


# Predict job resources per file 
ids = intersect(names(len),names(nseq))
b = data.frame(mlen = len[ids],nseq = as.integer(nseq[ids]))
ids = readLines('ids.txt')
b = b[rownames(b) %in% ids,]
options(max.print = 30)
o = b
for(job_name in names(models)){
	o = cbind(o,predict_res(models[[job_name]]$mem,models[[job_name]]$time,o))
	colnames(o)[(ncol(o)-1):ncol(o)] = paste0(job_name,'_',colnames(o)[(ncol(o)-1):ncol(o)])
}
rownames(o) = rownames(b)
o%>%arrange(-ALN_time)
o%>%arrange(-ALN_mem)

plotname = 'test.pdf'
w = 6
pdf(plotname,height = w, width = w)

plot(p$mem_used,predict_res(m3,t2,p)$mem,pch = 16,cex = 0.5,ylab = 'predicted',xlab = 'used');abline(0,1,col = 'blue')
dev.off()


