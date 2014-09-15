library(gplots)
#library(randomForest)
#library(mRMRe)
data=read.csv(file="miRNA_comb.csv", header=TRUE)
ndata=data[,-1]
nsampl=2
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
rownames(pdata)=c("NT", "1", "4", "7", "14", "28")
adata=pdata[, pval2 < 0.01]
colnames(adata)
rows=ncol(adata)
pdf("miRNA_comb_pooled.pdf",width=7, height=7)
heatmap.2(t(data.matrix(adata)), scale="none", Colv=FALSE, dendrogram="row",
	col=greenred(32),density.info="none",trace="none", cexRow=0.7,
	cexCol=1, ylab="log Expression", symbreaks=TRUE, 
	keysize=1,margins=c(7,7),
	main="Custom Panel\nF-test p < 0.01")
dev.off()
#miRNA_acc=read.csv("miRNA_acc.csv", header=TRUE)
#alist=sapply(accs, function(x) strsplit(x, '[.]')[[1]][1])
#rownames(miRNA_acc)=miRNA_acc[,1]
#alist=sapply(as.character(miRNA_acc[colnames(ndata),2]), function(x) strsplit(x, '[.]')[[1]][1])
#oufile=file("miRNA_all.txt")
#writeLines(alist, outfile)
#close(outfile) 
#accs=as.character(miRNA_acc[colnames(adata),2])
#alist=sapply(accs, function(x) strsplit(x, '[.]')[[1]][1])
alist=colnames(adata)
outfile=file("miRNA_out.txt")
writeLines(alist, outfile)
close(outfile)