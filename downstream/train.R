
#!/usr/bin/env Rscript

library(argparse)
library(stringr)
library(data.table)
library(quantreg)
library(jsonlite)

.trace_conv_time <- function(x){
  sapply(x,function(z){
    z=gsub(" ","",z)
    m=regmatches(z,regexec("^(?:(\\d+(?:\\.\\d+)?)m)?(?:(\\d+(?:\\.\\d+)?)s)?$",z))[[1]]
    if(length(m)==0) return(NA_real_)
    mins=ifelse(m[2]=="",0,as.numeric(m[2]))
    secs=ifelse(m[3]=="",0,as.numeric(m[3]))
    mins*60+secs
  })
}

.trace_conv_mem <- function(x){
  sapply(x,function(z){
    z=gsub(" ","",z)
    m=regmatches(z,regexec("^(\\d+(?:\\.\\d+)?)(MB|GB)$",z))[[1]]
    if(length(m)==0) return(NA_real_)
    value=as.numeric(m[2])
    unit=m[3]
    if(unit=="GB") value*1024 else value
  })
}

.trace_parse_name <- function(x){
  data.frame(
    job_name = str_split(x,' ',simplify=TRUE)[,1],
    id = gsub('\\(|\\)','',str_split(x,' ',simplify=TRUE)[,2])
  )
}

# -----------------------------
# CLI arguments
# -----------------------------

parser <- ArgumentParser(description="Train resource prediction models from Nextflow trace")

parser$add_argument("--trace", required=TRUE)
parser$add_argument("--seq_stats", required=TRUE)
parser$add_argument("--outfile", required=TRUE)
parser$add_argument("--plotfile", required=FALSE, default = 'plot.pdf',
                    help="Output PDF with model diagnostics")
parser$add_argument("--tau", type="double", default=0.95)
parser$add_argument("--min_n", type="double", default=50)
parser$add_argument("--jobs", default=NULL)

args <- parser$parse_args()

min_n = as.integer(args$min_n) # minimum number of sequences in fasta (to avoid fitting on small files)


jobs = args$jobs
if(!is.null(jobs)){jobs <- strsplit(jobs,",")[[1]]}

cat("Quantile tau:", args$tau, "\n")

# -----------------------------
# Load sequence stats
# -----------------------------

input <- fread(args$seq_stats, header=FALSE)
colnames(input) <- c("id","nseq","mlen")

# -----------------------------
# Load trace
# -----------------------------

trace <- fread(args$trace, sep="\t", header=TRUE)

trace <- trace[trace$status=="COMPLETED",]

trace <- cbind(trace, .trace_parse_name(trace$name))

trace$time <- .trace_conv_time(trace$realtime)
trace$mem  <- .trace_conv_mem(trace$peak_rss)

d <- trace[,c("job_name","id","time","mem")]

if(is.null(jobs)){
  jobs = unique(trace$job_name)
}
cat("Training models for:", paste(jobs, collapse=", "), "\n")

d <- merge(d, input, by="id")

d = d[d$nseq>=min_n,]

# -----------------------------
# Ensure coefficients exist
# -----------------------------

terms <- c("(Intercept)","log(nseq)","log(mlen)","log(nseq):log(mlen)")

get_coefs <- function(model){
  cfs <- coef(model)
  out <- setNames(rep(0,length(terms)),terms)
  out[names(cfs)] <- cfs
  as.list(out)
}

# -----------------------------
# Train models
# -----------------------------

models <- list()

pdf(args$plotfile,width = 10,height = 5)

for(job in jobs){

	dd <- d[d$job_name==job]

	dd <- dd[
	is.finite(mem) &
	is.finite(time) &
	nseq > 0 &
	mlen > 0
	]

  if(nrow(dd)==0){
    warning(paste("No data for job",job))
    next
  }

  cat("Fitting model for",job,"(",nrow(dd),"observations)\n")

  m1 <- rq(log(mem) ~ log(nseq) * log(mlen), data=dd, tau=args$tau)
  m2 <- rq(log(time) ~ log(nseq) * log(mlen), data=dd, tau=args$tau)

  models[[job]] <- list(
    mem  = get_coefs(m1),
    time = get_coefs(m2)
  )

  # -----------------------------
  # Diagnostic plots
  # -----------------------------

  pred_mem <- exp(predict(m1))
  pred_time <- exp(predict(m2))

  par(mfrow=c(1,2))

  plot(dd$mem, pred_mem,
       log="xy",
       xlab="Observed memory (MB)",
       ylab="Predicted memory (MB)",
       main=paste(job,"memory"),
       pch=16)

  abline(0,1,col="red",lwd=2)

  plot(dd$time, pred_time,
       log="xy",
       xlab="Observed time (s)",
       ylab="Predicted time (s)",
       main=paste(job,"runtime"),
       pch=16)

  abline(0,1,col="red",lwd=2)

}

dev.off()

# -----------------------------
# Save JSON
# -----------------------------

write_json(models, args$outfile, pretty=TRUE, auto_unbox=TRUE)

cat("Saved models to:", args$outfile, "\n")
cat("Saved plots to:", args$plotfile, "\n")