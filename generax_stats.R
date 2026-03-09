##########################
# Generax Runtimes
##########################
options(max.print = 30)
source('helper.R')
library(dplyr)
library(quantreg)
library(ggplot2)
set_theme()
# parse generax logs
f = list.files('results/generax',full = TRUE,pattern = '*log')
message(sprintf("%s generax logs found",length(f)))
p = do.call(rbind,lapply(f,parse_generax_log))
# Runtimes
df = p%>%group_by(id)%>%slice_max(radius)
dim(df)
df$time = df$time_total
range(df$time/60)
summary(df$n_genes)
fit <- lm(log10(time/60) ~ log10(n_genes), data = df)
r2  <- round(summary(fit)$r.squared, 3)
# that's very few genes yet, of course
df%>%mutate(min = time/60)%>%arrange(-min)%>%as.data.frame()
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

#########################
# Likelihoods
##########################
id = 'tfs.Forkhead.HG1'
# relative change - relative to the radius one
pp = p[,] %>%
  group_by(id) %>%
  select(-ncpu,-n_species) %>%
  arrange(radius, .by_group = TRUE) %>%
  mutate(rel_change = (JointLL - first(JointLL)) / abs(first(JointLL)))
pl = ggplot(pp,aes(x = radius, y = rel_change,group = id))+
geom_line(alpha = 0.5)+
scale_y_log10()+
ylab("Relative JointLL change")
plotname = 'downstream/generax.jointll_radius.pdf'
ggsave(plotname,height = 5,width = 5)
plotname = gsub(".pdf",".png",plotname)
ggsave(plotname,height = 5,width = 5,bg = 'white')
open(plotname)


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


########################################
# Count the number of orthogroups before and after GeneRax 
# Does it at least improve the gene tree topology? 
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
# high and low starness tre