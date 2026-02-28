#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(argparse)
  library(stringr)
  library(dplyr)
})

parser <- ArgumentParser(description = "Predict resources and write resources TSV")

parser$add_argument("--ids_fn", required = TRUE, help = "Path to ids.txt")
parser$add_argument("--cluster_dir", required = TRUE, help = "Directory with cluster FASTA files")
parser$add_argument("--models_rds", required = TRUE, help = "Path to models.rds")
parser$add_argument("--outfile", required = TRUE, help = "Output TSV path")

parser$add_argument("--min_mem", type = "integer", required = TRUE, help = "Minimum memory in MB")
parser$add_argument("--min_time", type = "integer", required = TRUE, help = "Minimum time in minutes")
parser$add_argument("--max_mem", type = "integer", required = TRUE, help = "Maximum memory in MB")
parser$add_argument("--max_time", type = "integer", required = TRUE, help = "Maximum time in minutes")
parser$add_argument("--increase", type = "double", required = TRUE, help = "Increase factor, e.g. 0.1 for +10%")

parser$add_argument("--default_aln_mem", type = "integer", default = 300, help = "Default ALN memory in MB")
parser$add_argument("--default_aln_time", type = "integer", default = 10, help = "Default ALN time in minutes")
parser$add_argument("--default_phy_mem", type = "integer", default = 300, help = "Default PHY memory in MB")
parser$add_argument("--default_phy_time", type = "integer", default = 10, help = "Default PHY time in minutes")
parser$add_argument("--default_pvm_mem", type = "integer", default = 500, help = "Default PVM memory in MB")
parser$add_argument("--default_pvm_time", type = "integer", default = 2, help = "Default PVM time in minutes")

args <- parser$parse_args()

ids_fn <- args$ids_fn
cluster_dir <- args$cluster_dir
models_rds <- args$models_rds
outfile <- args$outfile

defaults <- setNames(
  c(args$default_aln_mem, args$default_aln_time,
    args$default_phy_mem, args$default_phy_time,
    args$default_pvm_mem, args$default_pvm_time),
  c("aln_mem", "aln_time", "phy_mem", "phy_time", "pvm_mem", "pvm_time")
)

min_mem <- args$min_mem
min_time <- args$min_time
max_mem <- args$max_mem
max_time <- args$max_time
increase <- args$increase

ids <- readLines(ids_fn)
f <- list.files(cluster_dir, full.names = TRUE, pattern = "fasta")
names(f) <- gsub("\\.fasta$", "", basename(f))
f <- f[intersect(ids, names(f))]
message(sprintf("%s ids. %s .fastas found", length(ids), length(f)))

predict_res <- function(model_m, model_t, input) {
  pred_m <- exp(predict(model_m, newdata = input))
  pred_t <- exp(predict(model_t, newdata = input))
  data.frame(mem = pred_m, time = pred_t)
}

.get_nseq <- function(x) setNames(length(grep(">", readLines(x))), "nseq")

.get_mlen <- function(x) {
  lines <- readLines(x)
  seq_lines <- lines[!grepl("^>", lines)]
  setNames(median(nchar(seq_lines)), "mlen")
}

input <- t(sapply(f, FUN = function(x) c(.get_nseq(x), .get_mlen(x))))
input <- as.data.frame(input)

models <- readRDS(models_rds)

o <- lapply(names(models), FUN = function(job_name) {
  model_m <- models[[job_name]][["mem"]]
  model_t <- models[[job_name]][["time"]]
  u <- predict_res(model_m, model_t, input)
  u$mem[u$mem < min_mem] <- min_mem
  u$time[u$time < min_time] <- min_time
  u
})

names(o) <- tolower(names(models))

for (i in seq_along(o)) {
  colnames(o[[i]]) <- paste0(names(o)[i], "_", colnames(o[[i]]))
}

o <- do.call(cbind, setNames(o, NULL))

pred <- o
for (i in seq_along(pred)) {
  pred[, i] <- ceiling(pred[, i] * (1 + increase))
}

for (x in names(defaults)) {
  if (!x %in% colnames(pred)) {
    pred[, x] <- defaults[x]
  }
}

timcols <- grep("_time", colnames(pred), value = TRUE)
memcols <- grep("_mem", colnames(pred), value = TRUE)

pred <- as.data.frame(apply(pred, 2, ceiling))

rng <- t(apply(pred, 2, range))
msg <- paste0(
  sapply(
    1:nrow(rng),
    FUN = function(i) sprintf(
      "%s: %s - %s,%s",
      rownames(rng)[i], rng[i, 1], rng[i, 2],
      ifelse(grepl("time", rownames(rng)[i]), "min", "MB")
    )
  ),
  collapse = "\n"
)
message(msg)

b1 <- t(apply(pred[, memcols, drop = FALSE] >= max_mem, 1, FUN = function(x) c(sum(x), length(x))))
b2 <- t(apply(pred[, timcols, drop = FALSE] >= max_time, 1, FUN = function(x) c(sum(x), length(x))))

if (colSums(b1)[1] > 0) {
  message(sprintf("%s / %s jobs with MEM >= max_mem (%s)", colSums(b1)[1], sum(b1), max_mem))
}
if (colSums(b2)[1] > 0) {
  message(sprintf("%s / %s jobs with TIME>= max_time (%s)", colSums(b2)[1], sum(b2), max_time))
}

for (i in memcols) {
  pred[, i][pred[, i] >= max_mem] <- max_mem
}
for (i in timcols) {
  pred[, i][pred[, i] >= max_time] <- max_time
}

round_base <- function(x, base = 60) {
  ifelse(x >= base, ceiling(x / base) * base, x)
}

for (i in timcols) {
  pred[, i] <- round_base(pred[, i], base = 60)
}
for (i in memcols) {
  pred[, i] <- round_base(pred[, i], base = 1024)
}

convert_mem <- function(x) paste0(ceiling(x / 100) * 100, ".MB")
convert_time <- function(x, min_time = 5) paste0(ceiling(pmax(x, min_time)), ".min")

for (i in grep("_mem", colnames(pred), value = TRUE)) {
  pred[, i] <- convert_mem(pred[, i])
}
for (i in grep("_time", colnames(pred), value = TRUE)) {
  pred[, i] <- convert_time(pred[, i])
}

pred <- cbind(data.frame(id = rownames(pred)), pred)
rownames(pred) <- NULL

write.table(pred, outfile, sep = "\t", quote = FALSE, row.names = FALSE)
message(sprintf("Created: %s", outfile))