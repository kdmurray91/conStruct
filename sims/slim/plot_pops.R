#!/usr/bin/env Rscript

usage <- "
    ./plot_pops.R (directory with .indiv.tsvs and .pop.tsv file in)
"

args <- if (interactive()) { scan(what='') } else { commandArgs(TRUE) }
if (length(args) != 1) {
    stop(usage)
}
basedir <- args[1]
outfile <- file.path(basedir, "pop_locations.pdf")


# use fact that pops are written out in order to individual file
indivs <- read.table(file.path(basedir, "results.indiv.tsv"), header=TRUE)
pops <- read.table(file.path(basedir, "results.pop.tsv"), header=TRUE)
indivs$pop <- rep(pops$pop, pops$num_individuals)

pop_cols <- rainbow(1.5*nrow(pops))[1:nrow(pops)]
names(pop_cols) <- pops$pop


pdf(file=outfile, width=5, height=5, pointsize=10)
plot(indivs$location_x, indivs$location_y, col=pop_cols[indivs$pop], pch=20, asp=1,
     xlab='eastings', ylab='northings')
points(pops$x, pops$y, pch=21, cex=4,
       col=adjustcolor(pop_cols, 0.4),
       bg=adjustcolor(pop_cols, 0.4))
text(pops$x, pops$y, pops$pop)
dev.off()
