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
dim(df)
table(df$job_name)

# ALN
d = df[df$exit == 0,]
#d = d[d$job_name == "ALN",]
#d$mafft = sapply(rownames(d),FUN = function(hash){
#	grep('mafft',read_command_err(workdir,hash),value = TRUE)[2]
#})
# count the number of sequences in the input fastas 
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
d%>%group_by(job_name)%>%summarize(
	mem = sum(ReqMem_MB)/1024,time = sum(Timelimit_min)/60,
	mem_perc = sum(mem_perc/n(),na.rm = TRUE),
	time_perc = sum(time_perc/n(),na.rm = TRUE),
	mem_waste = mem * (1-mem_perc),
	time_waste = time * (1-time_perc))

summary(d[d$job_name == 'ALN',]$peak_rss_mb)
summary(d[d$job_name == 'ALN',]$ReqMem_MB)
# Median is 1Gb but it only uses around 300 
# ALN job requests too much memory 	
# I overask for the time clearly 
####################################
# Plots
####################################
plotname = 'test.pdf'
pdf(plotname,height = 5,width = 10)
par(mfrow = c(1,2))
plot(d$n,d$minutes,pch = 16,log = 'xy',xlab = '# sequences', ylab = 'Duration, min',font.main = 1,main = 'ALN')
plot(d$n,d$peak_rss_mb,pch = 16,log = 'xy',xlab = '# sequences', ylab = 'Memory, MB',font.main = 1,main = 'ALN')
dev.off()
options(max.print = 100)
d%>%filter(n>=1000)
i = 1
rownames(d) = d$hash

table(is.na(d$mafft))

list.files(get_job_workdir(workdir,hash),full = TRUE,pattern = 'command')
# ok, in the fast mode, of course there are almost no resources requested 

####################################
# Fraction of requested resources 
####################################



# mem
b1 = d$peak_rss_mb/d$ReqMem_MB
# time
b2 = d$minutes/d$Timelimit_min


plotname = 'test.pdf'
pdf(plotname,height = 5,width = 5)
par(mfrow = c(1,1))
plot(density(b1[!is.na(b1)]),font.main = 1,main = 'ALN',ylab = 'Density',xlab = '% of requested memory used')
plot(density(b2[!is.na(b2)]),font.main = 1,main = 'ALN',ylab = 'Density',xlab = '% of requested time')
plot(b1*100,b2*100,pch = 16,xlab = '% mem',ylab = '% time')
dev.off()

####################################
# PHY
####################################
table(df$job_name)
d = df[df$job_name == 'PHY',]

ids = unique(d$native_id)
cmd = sprintf(
	"sacct -j %s --format=JobIDRaw,ReqMem,Timelimit --noheader -P",
	paste(ids, collapse=",")
)
res = read.table(
	text = system(cmd, intern=TRUE),
	sep="|",
	stringsAsFactors=FALSE,
	col.names=c("JobIDRaw","ReqMem","Timelimit")
)
res = res[res$ReqMem!="",]
res$ReqMem_MB = as.numeric(sub("M","",sub("G","000",res$ReqMem)))
d = merge(d, res, by.x="native_id", by.y="JobIDRaw", all.x=TRUE)
d$Timelimit_min = as.numeric(substr(d$Timelimit,1,2))*60 +
                  as.numeric(substr(d$Timelimit,4,5)) +
                  as.numeric(substr(d$Timelimit,7,8))/60
d$mem_perc = d$peak_rss_mb/d$ReqMem_MB
d$time_perc = d$minutes/d$Timelimit_min

# plots 
plotname = 'test.pdf'
pdf(plotname,height = 5,width = 5)
plot(d$n,d$minutes,pch = 16,main = 'PHY',font.main = 1,log = 'xy',xlab = '# sequences',ylab = 'Time, min')
plot(density(d$time_perc[!is.na(d$time_perc)]),xlab =  '% time')
plot(density(d$mem_perc[!is.na(d$mem_perc)]),xlab =  '% time')
plot(d$mem_perc,d$time_perc,xlab = 'mem perc',ylab = 'time perc',pch = 16)
dev.off()
d$peak_rss

summary(d$Timelimit_min)
summary(d$minutes/d$Timelimit_min)
