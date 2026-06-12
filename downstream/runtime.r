# Do the stored models store the results properly
model_fn = 'workflow/models/models.rds'
trace_fn = 'reports/trace.step2.txt'
models = readRDS(model_fn)
names(models)
trace = read.table(trace_fn,sep = '\t',header = TRUE)
trace$job_name = str_split(trace$name,' ',simplify = T)[,1]
trace$hg_id = gsub('\\(|\\)','',str_split(trace$name,' ',simplify = T)[,2])

options(max.print = 100)
table(trace$job_name)
res = read.table('resources.tsv',header = TRUE)
#####################################
# 
#####################################
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
	message(cmd)
	
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
################################
# Explore the job stats - are they predicted correctly? 
################################
workdir = '~/no_backup/asebe/gzolotarov/nextflow/phylohpc/work_step2'
job_name = 'ALN'
d = trace[trace$job_name == job_name,]
d = d[d$status == 'COMPLETED',]
hash = 'c4/82bf46'
l = get_job_workdir(workdir = workdir,hash = hash)
l
# aha, can not access 
d$time = .trace_time_convert(d$duration)
summary(d$time)
# so prediction vs requested? 
