##########################
# HG sequence info - number of seqs, median length
##########################
source('helper.R')
trace_fn = 'reports/trace.step2.txt'
ids = read.table(trace_fn,sep = '\t')[[4]]
ids = gsub("\\(|\\)","",str_split(ids,' ',simplify = T)[,2])
f = list.files('results/clusters',full = TRUE,pattern = 'fasta')
names(f) = gsub('.fasta','',basename(f))
#count_seq = function(x){as.integer(system(sprintf('grep -c ">" %s',x),intern = TRUE))}
f = f[intersect(ids,names(f))]


n_seq = count_seq(f)
library(Biostrings)
# Median sequence length
mlen = sapply(f,FUN = function(x) median(width(readAAStringSet(x))))
length(mlen)

##########################
# Memory
##########################
source('helper.R')
trace_fn = 'reports/trace.step2.txt'
trace = read.table(trace_fn,sep = '\t',header = TRUE)
trace$job_name = str_split(trace$name,' ',simplify = T)[,1]
trace$hg_id = gsub('\\(|\\)','',str_split(trace$name,' ',simplify = T)[,2])
#d = trace[trace$job_name == job_name,]
d = trace
d = d[d$status == 'COMPLETED',]
dim(d)
summary(.trace_time_convert(d$duration))
d$time = .trace_time_convert(d$realtime)
d$mem = .trace_mem_convert(d$peak_rss)
d$nseq = n_seq[d$hg_id]
d$mlen = mlen[d$hg_id]
d$Class = str_split(d$hg_id,'\\.',simplify = T)[,1]
levs = c('ALN','PHY','GR_watcher','PVM')
pal = setNames(brewer.pal(length(levs),'Spectral'),levs)
d$job_name = factor(d$job_name,levels = levs)

res_fn = 'resources.tsv'
res = read.table(res_fn,header = TRUE)
library(RColorBrewer)
library(tidyr)
library(dplyr)

r = res %>%
  pivot_longer(
    cols = -id,
    names_to = c("job_name", ".value"),
    names_pattern = "(.+)_(mem|time)"
  )
r$job_name = toupper(r$job_name)
r$job_name[r$job_name == 'GR_WATCHER'] = 'GR_watcher'
colnames(r)[colnames(r) == 'mem'] = 'mem_res'
colnames(r)[colnames(r) == 'time'] = 'time_res'
d$id = gsub("\\(|\\)","",str_split(d$name,' ',simplify = T)[,2])
d = left_join(d,r,by = c("job_name","id"))


d$nseq = n_seq[d$id]
# resources contain the predictions 
d$time_pred = .nf_convert_time(d$time_res)
d$mem_pred = .nf_convert_mem(d$mem_res)

library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
pl1=ggplot(d,aes(x = nseq,y = time/60))+
geom_point(aes(col = job_name),size = 0.5)+
scale_x_log10()+
scale_y_log10()+
xlab('# nseq')+
ylab('Time, min')+
scale_color_manual(values = pal,breaks = names(pal))

pl2=ggplot(d[d$mem>=100,],aes(x = nseq,y = mem))+
geom_point(aes(col = job_name),size = 0.5)+
scale_x_log10()+
scale_y_log10()+
xlab('# nseq')+
ylab('mem, MB')+
scale_color_manual(values = pal,breaks = names(pal))
library(patchwork)
pl = pl1+pl2
plotname = 'downstream/scaling.pdf'
ggsave(plotname,width = 10,height = 5)
open(plotname)
plotname = gsub(".pdf",".png",plotname)
ggsave(plotname,width = 10,height = 4,bg = 'white')
open(plotname)

options(max.print = 30)

# Prediction efficiency - compare the resources data frame  
pl1 = ggplot(d[,],aes(x = job_name,y = time/time_pred))+
geom_jitter(aes(col = job_name))+
ylab('Time efficiency')+
scale_color_manual(values = pal,breaks = names(pal))+
geom_hline(yintercept = 1)


pl2 = ggplot(d[,],aes(x = job_name,y = mem/mem_pred))+
geom_jitter(aes(col = job_name))+
ylab('Memory efficiency')+
scale_color_manual(values = pal,breaks = names(pal))+
geom_hline(yintercept = 1)
pl1+pl2

# Generax vs Phylogeny time 
library(reshape2)
u = dcast(d,hg_id + Class + nseq + mlen~ job_name,value.var = 'time')
colnames(u)[grep("GR_",colnames(u))] = 'GR'
u$GR = u$GR/60
u$PHY = u$PHY/60
u = u[!is.na(u$GR)&!is.na(u$PHY),]
lims =range(c(u$PHY, u$GR), na.rm = TRUE)
ggplot(u, aes(x = PHY, y = GR)) +
  geom_point() +
  scale_x_log10(limits = lims) +
  scale_y_log10(limits = lims) +
  coord_fixed()+
  geom_abline(col = 'blue')+
  ylab('GeneRax runtime, min')+
  xlab('IQTREE2 runtime, min')
u%>%arrange(-log2(GR/PHY))



#########################################
# SLURM job stats - memory efficiency etc 
#########################################
fn = 'job_stats.tab'
u = read.table(fn)
colnames(u) = c('native_id','status','exitcode','time','time_req','mem_req','mem')
u$time = as.numeric(as.difftime(u$time, format="%H:%M:%S",'mins'))
u$time_req = as.numeric(as.difftime(u$time_req, format="%H:%M:%S",'mins'))

slurm_mem_conv = function(x) as.numeric(gsub("M$","",x))
u$mem = slurm_mem_conv(u$mem)
u$mem_req = slurm_mem_conv(u$mem_req)
u$mem_f = u$mem/u$mem_req
u$time_f = u$time/u$time_req
a = left_join(d[,c('job_name','native_id')],u)
a$job_name = factor(a$job_name,levels = names(pal))
pl1 = ggplot(a,aes(x = job_name, y = mem_f))+
#geom_violin(aes(fill = job_name))+
stat_boxplot(geom = 'errorbar',width = 0.25)+
geom_boxplot(outlier.shape = NA,aes(fill = job_name))+
scale_fill_manual(values = pal,breaks = names(pal))+
ylab("Fraction of memory used")

pl2 = ggplot(a,aes(x = job_name, y = time_f))+
#geom_violin(aes(fill = job_name))+
stat_boxplot(geom = 'errorbar',width = 0.25)+
geom_boxplot(outlier.shape = NA,aes(fill = job_name))+
scale_fill_manual(values = pal,breaks = names(pal))+
ylab("Fraction of time used")
pl1+pl2

