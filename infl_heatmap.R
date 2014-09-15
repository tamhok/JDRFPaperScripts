library(gplots)
#library(randomForest)
#library(mRMRe)
data=read.csv(file="inflammation.csv", header=TRUE)
ndata=data[,-1]
nsampl=4
scaling=rowMeans(ndata[,seq(1,nsampl)])
ndata=log(ndata / scaling)
ndata=as.data.frame(t(ndata))
colnames(ndata)=data[,1]
outc=ordered(as.character(c(rep("Mok",nsampl), rep(1,nsampl), 
	rep(4,nsampl), rep(7,nsampl),rep(14,nsampl), rep(28,nsampl))))
pval2=apply(ndata, 2, function(x) summary(aov(x~outc))[[1]]$"Pr(>F)"[1])
pdata=matrix(ncol=ncol(ndata), nrow=nrow(ndata)/nsampl)
for (i in 1:6) {
	pdata[i,]=log(colMeans(exp(ndata[seq(nsampl*(i-1)+1, nsampl*i),])))
}
colnames(pdata)=colnames(ndata)
rownames(pdata)=c("0", "1", "4", "7", "14", "28")
adata=pdata[, pval2 < 0.01/ncol(ndata)]
rows=ncol(adata)
pdf("inflammation_pooled.pdf",width=7, height=7)
heatmap.2(t(data.matrix(adata)), scale="none", Colv=FALSE, dendrogram="row",
	col=greenred(32),density.info="none",trace="none", cexRow=0.7,
	cexCol=1, ylab="log Expression", symbreaks=TRUE, 
	keysize=1,margins=c(7,7),
	main="Custom Panel\nF-test p < 0.01")
dev.off()

sums=apply(ndata,2,sum)
rows=ncol(ndata)
pdf("infl_ord.pdf",width=7, height=rows*8/75+7)
heatmap.2(t(data.matrix(pdata[,order(sums,decreasing=TRUE)])), scale="none", Colv=FALSE, dendrogram="none",
	col=greenred(32),density.info="none",trace="none", cexRow=0.7,
	cexCol=1, ylab="log Expression", symbreaks=TRUE, 
	keysize=1, Rowv=FALSE,
	margins=c(7,7),  
	lhei=c(1, rows*8/75),
	main="Inflammation Panel")
dev.off()

infl_acc=read.csv("infl_acc.csv", header=TRUE)
accs=as.character(infl_acc[colnames(adata),2])
alist=sapply(accs, function(x) strsplit(x, '[.]')[[1]][1])
rownames(infl_acc)=infl_acc[,1]
alist=sapply(as.character(infl_acc[colnames(ndata),2]), function(x) strsplit(x, '[.]')[[1]][1])
outfile=file("infl_all.txt")
writeLines(alist, outfile)
close(outfile) 

alidest=sapply(accs, function(x) strsplit(x, '[.]')[[1]][1])
outfile=file("infl_out.txt")
writeLines(alist, outfile)
close(outfile)