## For plotting layer contributions
layer.colors <- c("blue", "red", "goldenrod1", "forestgreen", "darkorchid1", "deepskyblue", "darkorange1", "seagreen2", "yellow1", "black")

cluster.2.layer <- function(conStruct.results){
	names(conStruct.results[[1]][["posterior"]])[which(names(conStruct.results[[1]][[1]])=="cluster.params")] <- "layer.params"
		names(conStruct.results[[1]]$posterior$layer.params) <- gsub("Cluster","Layer",names(conStruct.results[[1]]$posterior$layer.params))
			for(k in 1:length(conStruct.results[[1]]$posterior$layer.params)){
				names(conStruct.results[[1]]$posterior$layer.params[[k]])[which(names(conStruct.results[[1]]$posterior$layer.params[[k]])=="cluster.cov")] <- "layer.cov"
			}
	names(conStruct.results[[1]][["MAP"]])[which(names(conStruct.results[[1]]$MAP)=="cluster.params")] <- "layer.params"
		names(conStruct.results[[1]]$MAP$layer.params) <- gsub("Cluster","Layer",names(conStruct.results[[1]]$MAP$layer.params))
			for(k in 1:length(conStruct.results[[1]]$MAP$layer.params)){
				names(conStruct.results[[1]]$MAP$layer.params[[k]])[which(names(conStruct.results[[1]]$MAP$layer.params[[k]])=="cluster.cov")] <- "layer.cov"
			}
	return(conStruct.results)
}

plot.layer.contribs <- function (dir) {
    sp.dirs <- sort(list.files(dir, "spK._[0-9]*", full.names=TRUE))
    load(file.path(sp.dirs[1], "cs_data.block.Robj"))
    res.files <- unlist(lapply(sp.dirs, list.files, "cs_conStruct.results.Robj", full.names=TRUE))
    output.list.sp <- lapply(res.files, function (fn) { 
                                 a <- load(fn)
                                 conStruct.results <- cluster.2.layer(get(a))
                                 for (j in 1:length(conStruct.results[[1]]$MAP$layer.params)) {
                                     names(conStruct.results[[1]]$MAP$layer.params[[j]])[4] <- "phi"
                                 }
                                 return(conStruct.results)  })
    Kvals <- as.integer(gsub("_.*", "", gsub(".*spK", "", sp.dirs)))
    K <- max(Kvals)
    laycon.sp <- matrix(0, nrow=K, ncol=length(Kvals))
    colnames(laycon.sp) <- Kvals
    for (j in seq_along(Kvals)) {
        k <- Kvals[j]
        csr <- output.list.sp[[j]][[1]]
        laycon.sp[1:k,j] <- calculate.layer.contribution(csr, data.block, 1:k)
    }

    {
        par(mar=c(4,4,3,0.5))
        barplot(laycon.sp,	
                col=layer.colors,
                xlab="", ylab="layer contribution")
            mtext(side=1, text="number of layers", padj=4.25, adj=-1)
            mtext(side=3, text=paste("Layer contributions for", dir),
                  padj=-1.3, adj=12, font=2, cex=1.2)
    }
}

## For plotting the actual crossvalidation log-likelihoods

plot.sim.xvals <- function (dir, label, Kzoom=2) {
    if (missing(label)) { label <- dir }
    lnl.files <- list.files(dir, ".*.lnl.Robj", full.names=TRUE)
    n.reps <- length(lnl.files)
    x.vals <- lapply(lnl.files, function (x) { load(x); test.lnl })
    K <- length(x.vals[[1]][[1]])
	x.vals.std <- lapply(x.vals,conStruct:::standardize.xvals)
	xval.CIs <- conStruct:::get.xval.CIs(x.vals.std, K)

    layout(t(1:2))
	plot.xval.CIs(xval.CIs, K, jitter=0.1)
		legend(x="bottomright", pch=19, col=c("blue", "green"), legend=c("spatial", "nonspatial"))
	mtext("Predictive accuracy", side=2, padj=-5)
	plot.xval.CIs(xval.CIs, K, ylim=range(xval.CIs$sp.means[Kzoom:K]))
		legend(x="bottomleft", pch=c(19, NA), lty=c(NA, 1), lwd=c(NA, 2), col=c(1, adjustcolor(1, 0.8)), legend=c("mean", "95% CI"))
	mtext("number of layers", side=1, adj=-0.85, padj=4)
	mtext(paste("Cross-validation results", label), side=3, adj=9.25, padj=-2, font=2, cex=1.2)
}

plot.xval.CIs <- function(xval.CIs, K, k.range=c(1:K), ylim=NULL, cex=1.5, jitter=0, ...){
	if(is.null(ylim)){
		ylim <- range(c(unlist(lapply(k.range,function(k){xval.CIs$sp.CIs[[k]]})),
						 unlist(lapply(k.range,function(k){xval.CIs$nsp.CIs[[k]]}))))
	}
	plot(xval.CIs$sp.means,
			ylim=ylim,
			main= "",
			ylab="",
			xlab="",type='n',...)
		lapply(1:K,function(k){
				segments(k,xval.CIs$sp.CIs[[k]][1],
						 k,xval.CIs$sp.CIs[[k]][2],
						 col=adjustcolor(4,0.5),lwd=3)})
		lapply(1:K,function(k){
				segments(k+jitter,xval.CIs$nsp.CIs[[k]][1],
						 k+jitter,xval.CIs$nsp.CIs[[k]][2],
						 col=adjustcolor("green",0.5),lwd=3)})
		points(xval.CIs$sp.means,pch=19,col=4,cex=cex)			
		points(1:K + jitter, xval.CIs$nsp.means,col="green",pch=19,cex=cex)
	return(invisible("plotted"))
}
