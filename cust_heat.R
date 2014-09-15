library(gplots)
#library(randomForest)
#library(mRMRe)
data=read.csv(file="cap_custom_only.csv", header=TRUE)
ndata=data[,-1]
scaling=rowMeans(ndata[,seq(1,3)])
ndata=log(ndata / scaling)
ndata=as.data.frame(t(ndata))
colnames(ndata)=data[,1]
outc=ordered(as.character(c(rep("Mok",3), rep(1,3), rep(4,3), rep(7,3),
	 rep(14,3), rep(28,3))))
#nndata=ndata
#nndata$outc=outc
#rftree=randomForest(x=ndata, y=outc, ntree=5000, mtry=120, importance=TRUE)
#varImpPlot(rftree)
#adata=ndata[,importance(rftree)[,7] > 12]
#pval=apply(ndata, 2, function(x) kruskal.test(x, outc)$p.value)
pval2=apply(ndata, 2, function(x) summary(aov(x~outc))[[1]]$"Pr(>F)"[1])
adata=ndata[, pval2 < 0.01/ncol(ndata)]
pdata=matrix(ncol=ncol(ndata), nrow=nrow(ndata)/3)
for (i in 1:6) {
	pdata[i,]=log(colMeans(exp(ndata[seq(3*(i-1)+1, 3*i),])))
}
colnames(pdata)=colnames(ndata)
rownames(pdata)=c("Mok", "1", "4", "7", "14", "28")
adata=pdata[, pval2 < 0.01/ncol(ndata)]
rows=ncol(adata)
pdf("cust_pooled.pdf",width=7, height=7)
heatmap.2(t(data.matrix(adata)), scale="none", Colv=FALSE, dendrogram="row",
	col=greenred(32),density.info="none",trace="none", cexRow=0.7,
	cexCol=1, ylab="log Expression", symbreaks=TRUE, 
	keysize=1,margins=c(7,7),
	main="Custom Panel\nF-test p < 0.01")
dev.off()