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
o = sapply(f,FUN = function(x){
	system(sprintf('grep -c ">" %s',x),intern = TRUE)
})
names(o) = basename(names(o))
table(is.na(d$realtime))
names(o) = gsub(".fasta","",names(o))
nseq = o
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
print(min(d$submit_dt))
print(max(d$submit_dt))
sprintf("%d days %d hours %d minutes", days, hours, minutes)

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
df = d
d = df[df$exit == 0,]
# Slurm stats
res = .get_slurm_stats(ids = unique(d$native_id))
d = merge(d, res, all.x=TRUE)
# fraction of requensted memory used 
d$mem_perc = d$peak_rss_mb/d$ReqMem_MB
d$time_perc = d$minutes/d$Timelimit_min
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
d[d$status %in% status,]%>%filter(!is.na(peak_rss_mb))%>%group_by(job_name)%>%summarize(
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
# Plots
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

####################################
# Predict runtimes
####################################
# can we predicted the runtimes? 

# Seqeunce lengths
f = list.files('results/clusters/',pattern = 'fasta',full = TRUE)
names(f) = gsub('.fasta','',basename(f))
f = f[names(f) %in% d$hg_id]
lens = lapply(f,FUN = function(x) nchar(readLines(x)[!grepl('>',readLines(x))]))
len = sapply(lens,median)
# average sequence lengths? 
p = d[d$job_name == "ALN",]
summary(p$minutes)
p%>%arrange(-minutes)
p$median_len = len[p$hg_id]
m1 = lm(log(minutes) ~ log(n) + log(median_len), data=p)
m2 = lm(log(minutes) ~ I(2*log(n)) + log(median_len), data=p)
m3 = lm(log(minutes) ~ log(n) * log(median_len), data=p)
AIC(m1, m2,m3)
# interaction coefficients 
summary(m3)
pdf('test.pdf')
plot(predict(m3), log(p$minutes),pch = 16)
abline(0,1,col = 'blue')
dev.off()
# That's very good!
summary(predict(m3))
summary(exp(predict(m3)))
# less than a minute? 


p = p[!is.na(p$peak_rss_mb),]
m1 = lm(log(peak_rss_mb) ~ log(n) + log(median_len), data=p)
m2 = lm(log(peak_rss_mb) ~ I(2*log(n)) + log(median_len), data=p)
m3 = lm(log(peak_rss_mb) ~ log(n) * log(median_len), data=p)
AIC(m1, m2,m3)
# interaction coefficients 
summary(m3)
pdf('test.pdf')
plot(predict(m3), log(p$peak_rss_mb),pch = 16)
abline(0,1,col = 'blue')
dev.off()
# That's very good!
summary(predict(m3))
summary(exp(predict(m3)))
# less than a minute? 

#################################

fit_aln_models = function(p, time_col="minutes", rss_col="peak_rss_mb", n_col="n", L_col="median_len"){
	p = p[is.finite(p[[time_col]]) & p[[time_col]]>0 & is.finite(p[[rss_col]]) & p[[rss_col]]>0 & is.finite(p[[n_col]]) & p[[n_col]]>0 & is.finite(p[[L_col]]) & p[[L_col]]>0,]
	p$ln_n = log(p[[n_col]])
	p$ln_L = log(p[[L_col]])
	p$ln_t = log(p[[time_col]])
	p$ln_m = log(p[[rss_col]])
	m_time = lm(ln_t ~ ln_n*ln_L, data=p)
	m_mem  = lm(ln_m ~ ln_n*ln_L, data=p)
	list(time=m_time, mem=m_mem)
}
predict_aln_resources = function(models, n, median_len, time_mult=1.5, mem_mult=1.5, min_time_min=5, min_mem_mb=512){
	ln_n = log(n)
	ln_L = log(median_len)
	new = data.frame(ln_n=ln_n, ln_L=ln_L)
	new$t_min = pmax(min_time_min, exp(predict(models$time, newdata=new)) * time_mult)
	new$mem_mb = pmax(min_mem_mb, exp(predict(models$mem, newdata=new)) * mem_mult)
	new$mem_gb = new$mem_mb/1024
	new
}
