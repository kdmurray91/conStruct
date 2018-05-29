#!/usr/bin/env Rscript

usage <- "
Usage:
  xval_conStruct.R (data directory) (output directory) n.iter
"

args <- commandArgs(TRUE)
if (length(args) < 3) {
    cat(sprintf("Only %d arguments:\n", length(args)))
    cat(paste(c(args, "\n"), collapse=" : "))
    stop(usage)
}

basedir <- args[1]
outdir <- args[2]
n.iter <- as.numeric(args[3])
stopifnot(n.iter > 0)
stopifnot(file.exists(basedir)
          && file.exists(outdir))

library(conStruct)

# input
freqfile <- file.path(basedir, "results.freq.tsv")
popfile <- file.path(basedir, "results.pop.tsv")

if (!file.exists(freqfile)
    || !file.exists(popfile)) {
    stop("The relevant files do not exist.")
}

freqs <- t(read.table(freqfile, sep='\t', header=TRUE))
pops <- read.table(popfile, sep='\t', header=TRUE)

geog.coords <- as.matrix(pops[,c("x","y")])
geog.dist <- as.matrix(dist(geog.coords))

x.validation(n.reps=10,
          K=1:7,
          freqs=freqs,
          geoDist=geog.dist,
          coords=geog.coords,
          prefix=file.path(outdir, "xval"),
          n.iter=n.iter,
          make.figs=TRUE,
          save.files=TRUE)


