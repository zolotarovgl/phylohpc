library(stringr)
open = function(x){system(sprintf('open %s',x))}
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


##########################
# HG info
##########################
trace_fn = 'reports/trace.step2.txt'
ids = read.table(trace_fn,sep = '\t')[[4]]
ids = gsub("\\(|\\)","",str_split(ids,' ',simplify = T)[,2])
f = list.files('results/clusters',full = TRUE,pattern = 'fasta')
names(f) = gsub('.fasta','',basename(f))
count_seq = function(x){as.integer(system(sprintf('grep -c ">" %s',x),intern = TRUE))}
f = f[intersect(ids,names(f))]
count_seq_fast <- function(files) {
  as.integer(system(
    sprintf("grep -c '^>' %s | cut -f 2 -d :", paste(shQuote(files), collapse = " ")),
    intern = TRUE
  ))
}
n_seq <- count_seq_fast(f)
names(n_seq) <- names(f)
library(Biostrings)
library(stringr)
library(ggplot2)

mlen = sapply(f,FUN = function(x) median(width(readAAStringSet(x))))
length(mlen)
##########################
# Runtimes
##########################
# parse generax logs
f = list.files('results/generax',full = TRUE,pattern = '*log')
length(f)
p = do.call(rbind,lapply(f,parse_generax_log))
options(max.print = 30)
dim(p)

##########################
# Runtimes
##########################
library(ggplot2)
library(cowplot)
library(dplyr)
theme_set(
theme_cowplot()+
	theme(
    plot.title = element_text(
      face = "plain",  
      size = 20        
    ))
)
# count the number of sequences / the size of the tree? 

df = p%>%group_by(id)%>%slice_max(radius)
dim(df)
df$time = df$time_total
range(df$time/60)
summary(df$n_genes)
fit <- lm(log10(time/60) ~ log10(n_genes), data = df)
r2  <- round(summary(fit)$r.squared, 3)
ggplot(df, aes(n_genes, time/60)) +
	geom_point() +
	scale_x_log10() +
	scale_y_log10() +
	geom_smooth(method = "lm",
				formula = y ~ x,
				se = TRUE) +
	annotate("text",
			x = Inf, y = Inf,
			label = paste0("R² = ", r2),
			hjust = 1.1, vjust = 1.5)+
	xlab('# Genes')+
	ylab('Time, min')+
	ggtitle('GeneRax')
# that's very few genes yet, of course
df%>%mutate(min = time/60)%>%arrange(-min)%>%as.data.frame()
library(quantreg)
library(dplyr)
library(ggplot2)
library(quantreg)
df=p%>%group_by(id)%>%slice(which.max(radius))
fit=lm(log10(time/60)~log10(n_genes),data=df)
r2=round(summary(fit)$r.squared,3)
tau = 0.95
fit_q=rq(log10(time/60)~log10(n_genes),data=df,tau=tau)
newdata=data.frame(n_genes=seq(min(df$n_genes),max(df$n_genes),length.out=200))
pred_log=predict(fit_q,newdata=newdata)
newdata$time_pred=10^pred_log
ggplot(df,aes(n_genes,time/60))+
geom_point()+
scale_x_log10()+
scale_y_log10()+
geom_smooth(method="lm",formula=y~x,se=TRUE,color="blue")+
geom_line(data=newdata,aes(x=n_genes,y=time_pred),color="red",linewidth=1)+
annotate("text",x=Inf,y=Inf,label=paste0("R² = ",r2),hjust=1.1,vjust=1.5)+
xlab("# Genes")+
ylab("Time, min")+
ggtitle("GeneRax")
# what about memory? 

##########################
# Memory
##########################
.trace_time_convert=function(x){
  sapply(x,function(z){
    z=gsub(" ","",z)
    m=regmatches(z,regexec("^(?:(\\d+(?:\\.\\d+)?)m)?(?:(\\d+(?:\\.\\d+)?)s)?$",z))[[1]]
    if(length(m)==0) return(NA_real_)
    mins=ifelse(m[2]=="",0,as.numeric(m[2]))
    secs=ifelse(m[3]=="",0,as.numeric(m[3]))
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
library(RColorBrewer)
pal = setNames(brewer.pal(length(levs),'Spectral'),levs)
d$job_name = factor(d$job_name,levels = levs)

res_fn = 'resources.tsv'
res = read.table(res_fn,header = TRUE)
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

# resources contain the predictions 
.nf_convert_time = function(x){
	tmap = c('s' = 1,'min' = 60,'h' = 60^2)
	as.numeric(str_split(x,'\\.',simplify = T)[,1])*tmap[str_split(x,'\\.',simplify = T)[,2]]
}
.nf_convert_mem = function(x){
	as.integer(str_split(d$mem_res,'\\.',simplify = T)[,1])
}
d$time_pred = .nf_convert_time(d$time_res)
d$mem_pred = .nf_convert_mem(d$mem_res)

pl1=ggplot(d,aes(x = nseq,y = time/60))+
geom_point(aes(col = job_name))+
scale_x_log10()+
scale_y_log10()+
xlab('# nseq')+
ylab('Time, min')+
scale_color_manual(values = pal,breaks = names(pal))


pl2=ggplot(d[d$mem>=100,],aes(x = nseq,y = mem))+
geom_point(aes(col = job_name))+
scale_x_log10()+
scale_y_log10()+
xlab('# nseq')+
ylab('mem, MB')+
scale_color_manual(values = pal,breaks = names(pal))
library(patchwork)
pl1+pl2

#

ggplot(d,aes(x = nseq,y = time/60))+
geom_point(aes(col = job_name))+
scale_x_log10()+
xlab('# nseq')+
ylab('Time, min')+
scale_color_manual(values = pal,breaks = names(pal))

ggplot(d[,],aes(x = time/60, y = time_pred/60))+
geom_point(aes(col = job_name))+
geom_abline()+
scale_color_manual(values = pal,breaks = names(pal))


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
########

ggplot(d,aes(x = time, y = mem))+
geom_point(aes(col = job_name))+
scale_x_log10()+
scale_y_log10()+
xlab('Time, sec')

b = tapply(n_seq,str_split(names(n_seq),'\\.',simplify = T)[,1],sum)
b = sort(b,decreasing = TRUE)
ids = names(b[b>=2000])
ggplot(d[d$job_name == "ALN"&d$Class %in% ids,],aes(x = time, y = mem))+
geom_point(aes(col = Class))+
scale_x_log10()+
scale_y_log10()+
xlab('Time, sec')+
ylab("Mem, MB")


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
#########################
# Likelihoods
##########################
# Likelihood gains
options(max.print = 30) 
library(dplyr)
ggplot(p%>%group_by(id)%>%summarize(time = max(time),n_genes = mean(n_genes),n_species = mean(n_species),maxJointLL = max(JointLL),minJointLL = min(JointLL),dLL = max(JointLL)-min(JointLL)),
	aes(x = n_genes,y = dLL))+
	geom_point()+
	scale_y_log10()+
	ylab('Delta LL')+
	xlab('# Genes')+
	scale_x_log10()
p
df_ll = p%>%group_by(id)%>%summarize(time = max(time),n_genes = mean(n_genes),n_species = mean(n_species),maxJointLL = max(JointLL),minJointLL = min(JointLL),dJointLL = max(JointLL)-min(JointLL),dRecLL = max(RecLL)-min(RecLL),dLibpllLL = max(LibpllLL) - min(LibpllLL))
ggplot(df_ll,aes(x = time/60,y = dJointLL))+
geom_point()+
scale_y_log10()+
ylab('Delta LL')+
xlab('Time, min')+
scale_x_log10()

ggplot(df_ll,aes(x = n_genes,y = dJointLL))+
geom_point()+
scale_y_log10()+
ylab('Delta LL')+
xlab('# Genes')+
scale_x_log10()





#########################################
# Count the number of orthogroups before and after GeneRax 
# Does it at least improve the gene tree topology? 
# 
#########################################
# tmp directory for the test 
run_possvm = function(tree_fn,verbose = FALSE){
	log_fn = paste0(tree_fn,'.log')
	cmd = sprintf("python phylogeny/main.py possvm -t %s -l %s",tree_fn,log_fn)
	if(verbose){message(cmd)}
	system(cmd)
	return(log_fn)
}
get_nog_log = function(log_fn){
	o = readLines(log_fn)[grepl('n OGs',readLines(log_fn))]
	as.integer(gsub('n OGs = ','',str_extract(o[length(o)],'n OGs = [0-9]+')))
}

.compare = function(id,verbose = FALSE){
	system(sprintf('cp results/gene_trees/%s.treefile %s',id,dir))
	system(sprintf('cp results/generax/%s.generax.tree %s',id,dir))
	fn1 = sprintf('%s/%s.treefile',dir,id)
	log_fn = run_possvm(fn1)
	n1 = get_nog_log(log_fn)
	fn2 = sprintf('%s/%s.generax.tree',dir,id)
	log_fn = run_possvm(fn2)
	n2 = get_nog_log(log_fn)
	message(id)
	data.frame(id = id,n1 = n1, n2 = n2)
}

# more orthogroups? what? 
ids = unique(p$id)
library(parallel)
ncpu = 8
out = mclapply(ids,FUN = function(id){
	.compare(id)
},mc.cores = ncpu)
out = do.call(rbind,out)
out = as.data.frame(out)
ggplot(out,aes(x = n1,y = n2))+
geom_abline(a = 0,b = 1, col = 'grey')+
geom_point()+
coord_fixed()+
xlim(1,max(c(out$n1,out$n2)))+
ylim(1,max(c(out$n1,out$n2)))+
xlab('# OGs pre-generax')+
ylab('# OGs post-GeneRax')

#########################################
# Compare the tree topologies? 
#########################################
f = list.files(dir,pattern = 'tree$|treefile$',full = TRUE)
tr = lapply(f,ape::read.tree)
names(tr) = basename(f)
names(tr) = gsub('.treefile$','.input.tree',names(tr))
starness_bl <- function(tree) {
  if (is.null(tree$edge.length))
    stop("Tree has no branch lengths")

  # identify terminal edges
  tip_edges <- tree$edge[,2] <= Ntip(tree)

  terminal_sum <- sum(tree$edge.length[tip_edges])
  total_sum    <- sum(tree$edge.length)

  terminal_sum / total_sum
}
d = data.frame(fn = names(tr),Ntip = sapply(tr,Ntip),starness = sapply(tr,starness_bl))
d$type = str_split(d$fn,'\\.',simplify = T)[,4]
d$id = gsub('.generax.tree|.input.tree','',d$fn)
d$type = factor(d$type,levels = c('input','generax'))
ggplot(d,aes(x = type, y = starness))+
geom_line(aes(group = id),alpha = 0.5)+
geom_point()
library(reshape2)
u = dcast(d,id ~ type,value.var = 'starness')
u$dstar = u$input-u$generax
plot(density(u$dstar));abline(v = 0)
summary(u$dstar)

# so what? does generax INCREASE the starness?????
# wtf ????
# what's going on? 
# Clearly, the generax seesm to increase the starness of the trees
# high and low starness tree

id1 = d%>%arrange(-starness)%>%dplyr::pull(fn)
id2 = d%>%arrange(starness)%>%pull(fn)
id1 = id1[1]
id2 = id2[3]
par(mfrow = c(1,2))
plot(tr[[id1]], show.tip.label = FALSE)
plot(tr[[id2]], show.tip.label = FALSE)

#########################################
unique(df_ll$n_species)
df_ll%>%arrange(-dLL)
# Can we decide which trees makes sense to run the generax for? 
# also would be cool to run possvm for these cases to see if there are less orthogroups called there. 

