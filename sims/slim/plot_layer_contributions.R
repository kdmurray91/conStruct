#!/usr/bin/env Rscript

usage <- "
Usage:
  ./plot_layer_contributions.R (directory with construct results in)
"

args <- commandArgs(TRUE)
if (length(args) < 1) {
    stop(usage)
}

dir <- args[1]

library(conStruct)
source("construct_plot_fns.R")

outfile <- file.path(dir, "layer_contributions.pdf")

pdf(file=outfile, width=6, height=3, pointsize=10)
plot.layer.contribs(dir)
dev.off()
