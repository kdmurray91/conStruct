
## For plotting the actual crossvalidation log-likelihoods

plot.sim.xvals <- function (dir, label) {
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
	plot.xval.CIs(xval.CIs, K, ylim=range(xval.CIs$sp.means[2:K]))
		legend(x="bottomleft", pch=c(19, NA), lty=c(NA, 1), lwd=c(NA, 2), col=c(1, adjustcolor(1, 0.8)), legend=c("mean", "95% CI"))
	mtext("number of layers", side=1, adj=-0.85, padj=4)
	mtext(bquote(paste("Cross-validation results", label)), side=3, adj=9.25, padj=-2, font=2, cex=1.2)
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
