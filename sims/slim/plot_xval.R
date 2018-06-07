#!/usr/bin/env Rscript

usage <- "
Usage:
  ./plot_xval.R (directory with crossvalidation results in) (K value to zoom in above)
"

args <- commandArgs(TRUE)
if (length(args) < 2) {
    stop(usage)
}

source("construct_plot_fns.R")

dir <- args[1]
K <- as.numeric(args[2])
outfile <- file.path(dirname(dir), paste0(basename(dir), "_results.pdf"))

pdf(file=outfile, width=6, height=3, pointsize=10)
plot.sim.xvals(dir, Kzoom=K)
dev.off()
