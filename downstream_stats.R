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
d$hg_id = gsub('test_','',d$hg_id)
d$minutes = sapply(d$realtime, function(x){
	x = trimws(x)
	m = if(grepl("m", x)) as.numeric(sub("m.*", "", x)) else 0
	s = if(grepl("s", x)) as.numeric(sub(".*m\\s*", "", sub("s", "", x))) else as.numeric(sub("s", "", x))
	if(is.na(s)) s = 0
	m + s/60
})


# Sequence counts
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

d$n = as.integer(nseq[d$hg_id])
d$median_len = as.numeric(len[d$hg_id])
d = d[d$exit == 0,]

# Slurm stats
d$native_id
res = .get_slurm_stats(ids = unique(d$native_id))
d = merge(d, res, all.x=TRUE)

# fraction of requensted memory used 
d$mem_perc = d$peak_rss_mb/d$ReqMem_MB
d$time_perc = d$minutes/d$Timelimit_min
saveRDS(d,'results/downstream/job_info.rds')
###########################################
d = readRDS('results/downstream/job_info.rds')
summary(d$ReqMem_MB)
plotname = 'test.pdf'
pdf(plotname,height = 5,width = 10)
par(mfrow = c(1,2))
boxplot(split(d$mem_perc,d$job_name),ylab = 'Fraction of MEM used',outline = FALSE)
boxplot(split(d$time_perc,d$job_name),ylab = 'Fraction of TIME used',outline = FALSE)
dev.off()

library(cowplot)
theme_set(theme_cowplot())
pl = ggplot(d,aes(x = ReqMem_MB+1, y = mem_perc))+
geom_point(aes(col = job_name),size = 0.5)+
scale_x_log10()
ggsave(plotname)

####################################
# Resource scaling 
####################################
plotname = 'results/downstream/scaling.pdf'
pdf(plotname,height = 10,width = 10)
par(mfrow = c(3,2))
for(job_name in c("ALN","PHY","PVM")){
	p = d[d$job_name == job_name,]
	plot(p$n,p$minutes,pch = 16,log = 'xy',xlab = '# sequences', ylab = 'Duration, min',font.main = 1,main = job_name)
	plot(p$n,p$peak_rss_mb,pch = 16,log = 'xy',xlab = '# sequences', ylab = 'Memory, MB',font.main = 1,main = job_name)
}
dev.off()
options(max.print = 100)
dir.create('results/downstream')
####################################
# Predict runtimes
####################################
d = readRDS('results/downstream/job_info.rds')
cls = c('hg_id','job_name','hash','minutes','Timelimit_min','peak_rss_mb','ReqMem_MB')
cls[!cls %in% colnames(d)]
dd = d[,cls]
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
saveRDS(models,'workflow/models/models.rds')

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
id = 'tfs.Forkhead.HG5'
id = 'tfs.NFYB_NFYC.HG2'
id = 'tfs.T-box.HG1'
id = 'tfs.Forkhead.HG1'

o[id,]
# smth is off with the time scaling here 
convert_mem = function(x) paste0(ceiling(x/10)*10,'.MB')
convert_time = function(x,min_time = 5){paste0(ceiling(pmax(x,min_time)),'.min')}

#`id	aln_mem	aln_time	phy_mem	phy_time	pvm_mem	pvm_time`
out = cbind(data.frame(id = rownames(o)),o[,c('ALN_mem','ALN_time','PHY_mem','PHY_time')])
# add 10 %
out$ALN_time = out$ALN_time*1.1
out$ALN_mem = out$ALN_mem*1.1
out$PHY_time = out$PHY_time*1.1
out$PHY_mem = out$PHY_mem*1.1

out$ALN_time = ceiling(out$ALN_time/10)*10
out$PHY_time = ceiling(out$PHY_time/10)*10
summary(out$PHY_time/60)
out%>%arrange(-PHY_time)
max_time = 48
table(out$ALN_time>=max_time)
table(out$PHY_time>=max_time)
out$ALN_time[out$ALN_time>=max_time] = max_time
out$PHY_time[out$PHY_time>=max_time] = max_time


max_mem = 1024*10
table(out$ALN_mem>=max_mem)
table(out$PHY_mem>=max_mem)
out$ALN_mem[out$ALN_mem>=max_mem] = max_mem
out$PHY_mem[out$PHY_mem>=max_mem] = max_mem


out$ALN_mem = convert_mem(out$ALN_mem)
out$PHY_mem = convert_mem(out$PHY_mem)
min_time = 5
out$ALN_time = convert_time(out$ALN_time)
out$PHY_time = convert_time(out$PHY_time)
rownames(out) = NULL
colnames(out) = tolower(colnames(out))
write.table(out,'resources.tsv',sep = '\t',quote = FALSE,row.names = FALSE)
# PHY - needs about 38 mins of time 

o%>%arrange(-ALN_time)
o%>%arrange(-ALN_mem)

plotname = 'test.pdf'
w = 6
pdf(plotname,height = w, width = w)

plot(p$mem_used,predict_res(m3,t2,p)$mem,pch = 16,cex = 0.5,ylab = 'predicted',xlab = 'used');abline(0,1,col = 'blue')
dev.off()


