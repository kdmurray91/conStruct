#!/usr/bin/env Rscript

usage <- "
Usage:
  ./plot_xval.R (directory with crossvalidation results in)
"

args <- commandArgs(TRUE)
if (length(args) < 1) {
    stop(usage)
}

source("xval_plot_fns.R")

dir <- args[1]
outfile <- file.path(dirname(dir), paste0(basename(dir), "_results.pdf"))

pdf(file=outfile, width=6, height=3, pointsize=10)
plot.sim.xvals(dir)
dev.off()
