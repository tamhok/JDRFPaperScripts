data=read.csv(file="cust_nrmlz.csv", header=TRUE, stringsAsFactors=FALSE)
ndata=data[,-1]
#Scale to the mock data
ndata=as.data.frame(t(ndata))
colnames(ndata)=data[,1]
sdata=ndata
scaling=colMeans(ndata[1:3,])
sdata=t(log(apply(ndata,1,function(x) x/scaling)))
pdata=matrix(ncol=ncol(sdata), nrow=nrow(sdata)/3)
for (i in 1:6) {
	pdata[i,]=log(colMeans(exp(sdata[seq(3*(i-1)+1, 3*i),])))
}
colnames(pdata)=colnames(sdata)
rownames(pdata)=c("Mock", "Day 1", "Day 4", "Day 7", "Day 14", "Day 28")
rows=ncol(sdata)
sums=apply(ndata,2,sum)

source("pheatmap.R")

draw_map=function(data, filename=NA, sidelabel=NA, breaks=col_breaks, rbreaks=cbreaks, legend_breaks=c(-4,-2,0,2,4), legend_labels=c(-4,-2,0,2,4)) {
	if(!is.na(filename)) {
		width=nrow(data)*36/72+4
		tiff(filename, width=width, units="in", res=300, height=ncol(data)*20/72+3)
	}
	pheatmap(t(data.matrix(data)), col=greenred(64), cluster_rows=FALSE, cluster_cols=FALSE, legend_breaks=legend_breaks, legend_labels=legend_labels,
		family="Sans", fontsize=18, cellwidth=36, cellheight=20, main=sidelabel,
		breaks=breaks, rbreaks=rbreaks, scale="none", legend=TRUE)
	if(!is.na(filename)) {
		dev.off()
	}
}
draw_map(sdata, filename="cust_nrmlz_scaled.tiff")
draw_map(ndata, filename="cust_nrmlz_unscaled.tiff", breaks=NA, rbreaks=NA, legend_breaks=NA, legend_labels=NA)

data=read.csv(file="infl_nrmlz.csv", header=TRUE)
ndata=data[,-1]
rownames(ndata)=data[,1]
nsampl=4
scaling=rowMeans(ndata[,seq(1,nsampl)])
sdata=log(ndata / scaling)
ndata=as.data.frame(t(ndata))
sdata=as.data.frame(t(sdata))
pdata=matrix(ncol=ncol(ndata), nrow=nrow(ndata)/nsampl)
for (i in 1:6) {
	pdata[i,]=log(colMeans(exp(ndata[seq(nsampl*(i-1)+1, nsampl*i),])))
}
colnames(pdata)=colnames(ndata)
rownames(pdata)=c("Day 0", "Day 1", "Day 4", "Day 7", "Day 14", "Day 28")
rows=ncol(adata)

draw_map(sdata, filename="infl_nrmlz_scaled.tiff")
draw_map(ndata, filename="infl_nrmlz_unscaled.tiff", breaks=NA, rbreaks=NA, legend_breaks=NA, legend_labels=NA)
