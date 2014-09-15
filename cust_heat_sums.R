library(gplots)
data=read.csv(file="cap_custom_only.csv", header=TRUE)
ndata=data[,-1]
scaling=rowMeans(ndata[,seq(1,3)])
ndata=log(ndata / scaling)
ndata=as.data.frame(t(ndata))
colnames(ndata)=data[,1]
pdata=matrix(ncol=ncol(ndata), nrow=nrow(ndata)/3)
for (i in 1:6) {
	pdata[i,]=log(colMeans(exp(ndata[seq(3*(i-1)+1, 3*i),])))
}
colnames(pdata)=colnames(ndata)
rownames(pdata)=c("Mok", "1", "4", "7", "14", "28")
rows=ncol(ndata)
sums=apply(ndata,2,sum)

pdf("cust_ord.pdf",width=4, height=rows*8/75+2)
heatmap.2(t(data.matrix(ndata[,order(sums,decreasing=TRUE)])), scale="none", Colv=FALSE, dendrogram="none",
	col=greenred(32),density.info="none",trace="none", cexRow=0.7,
	cexCol=1, symbreaks=TRUE, family="sans",
	keysize=1,margins=c(7,7), Rowv=FALSE, lhei=c(1, rows*8/75))
dev.off()

pdf("cust_pooled_ord.pdf",width=4, height=rows*8/75+2)
heatmap.2(t(data.matrix(pdata[,order(sums,decreasing=TRUE)])), scale="none", Colv=FALSE, dendrogram="none",
	col=greenred(32),density.info="none",trace="none", cexRow=0.7,
	cexCol=1, symbreaks=TRUE, family="sans",
	keysize=1,margins=c(7,7), Rowv=FALSE, lhei=c(1, rows*8/75))
dev.off()

type_file=scan("cust_pops.csv", what="", sep="\n")
types=strsplit(type_file, ',+')
names(types)=sapply(types, '[[',1)
types=lapply(types, '[', -1)
col_breaks=seq(-4.5,4.5,length.out=33)
for(i in 1:length(types)) {
	tdata=ndata[,types[[i]]]
	sums=apply(tdata,2,sum)
	pdf(file=paste("cust_ord_", names(types)[i], ".pdf", sep=""),width=7, height=rows*8/75+2)
	heatmap.2(t(data.matrix(tdata[,order(sums,decreasing=TRUE)])), scale="none", Colv=FALSE, dendrogram="none",
		col=greenred(32),density.info="none",trace="none", cexRow=0.7,
		cexCol=1, ylab="log Expression", breaks=col_breaks, family="sans",
		keysize=1,margins=c(7,7), Rowv=FALSE, lhei=c(1, rows*8/75),
		main=paste("Custom Panel: ", names(types)[i], sep=""))
	dev.off()
}

for(i in 1:length(types)) {
	tdata=pdata[,types[[i]]]
	sums=apply(tdata,2,sum)
	pdf(file=paste("cust_ord_pooled", names(types)[i], ".pdf", sep=""),width=7, height=rows*8/75+2)
	heatmap.2(t(data.matrix(tdata[,order(sums,decreasing=TRUE)])), scale="none", Colv=FALSE, dendrogram="none",
		col=greenred(32),density.info="none",trace="none", cexRow=0.7,
		cexCol=1, breaks=col_breaks, family="sans",
		keysize=1,margins=c(7,7), Rowv=FALSE, lhei=c(1, rows*8/75))
	dev.off()
}


data=read.csv("cap_fat_mn.csv", header=TRUE)	
ndata=data.matrix(data[,seq(-3,-1)])
scaling=data[,3]
ndata=log(ndata / scaling)
ndata=as.data.frame(t(ndata))
colnames(ndata)=data[,1]
rows=ncol(ndata)
caps=ndata[1:5,];
fats=ndata[6:10,];
ndata=rbind(fats,caps)
sums=apply(ndata,2,sum)
pdf(file="cap_fat_cust_ord.pdf",width=7, height=rows*8/75+2) 
heatmap.2(t(data.matrix(ndata[,order(sums,decreasing=TRUE)])), scale="none", Colv=FALSE, dendrogram="none",
	col=greenred(32),density.info="none",trace="none", cexRow=0.7,
	cexCol=1, ylab="log Expression", symbreaks=TRUE, family="sans",
	keysize=1,margins=c(7,7), Rowv=FALSE, lhei=c(1, rows*8/75),
	main="Custom Panel\nCap vs Fat")
dev.off()

for(i in 1:length(types)) {
	tdata=ndata[,types[[i]]]
	sums=apply(tdata,2,sum)
	pdf(file=paste("cap_fat_ord_pooled", names(types)[i], ".pdf", sep=""),width=7, height=rows*8/75+2)
	heatmap.2(t(data.matrix(tdata[,order(sums,decreasing=TRUE)])), scale="none", Colv=FALSE, dendrogram="none",
		col=greenred(32),density.info="none",trace="none", cexRow=0.7,
		cexCol=1, ylab="log Expression", breaks=col_breaks, family="sans",
		keysize=1,margins=c(7,7), Rowv=FALSE, lhei=c(1, rows*8/75),
		main=paste("Custom Panel: ", names(types)[i], sep=""))
	dev.off()
}


