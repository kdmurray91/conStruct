#!/usr/bin/env Rscript

usage <- "
Usage:
  run_conStruct.R (data directory) (output directory) K n.iter n.chains
"

args <- commandArgs(TRUE)
if (length(args) < 5) {
    cat(sprintf("Only %d arguments:\n", length(args)))
    cat(paste(c(args, "\n"), collapse=" : "))
    stop(usage)
}

basedir <- args[1]
outdir <- args[2]
Kval <- as.numeric(args[3])
n.iter <- as.numeric(args[4])
n.chains <- as.numeric(args[5])
stopifnot(Kval > 0 && n.iter > 0 && n.chains > 0)
stopifnot(file.exists(basedir)
          && file.exists(outdir))

library(conStruct)

if (n.chains > 1) { options(mc.cores = n.chains) }

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

conStruct(spatial=TRUE,
          K=Kval,
          freqs=freqs,
          geoDist=geog.dist,
          coords=geog.coords,
          prefix=file.path(outdir, "cs"),
          n.chains=n.chains,
          n.iter=n.iter,
          make.figs=TRUE,
          save.files=TRUE)


