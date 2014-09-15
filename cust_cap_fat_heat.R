library(gplots)
data=read.csv(file="cap_fat_custom.csv", header=TRUE)
#Get actual data from frame, scale by average of NT mice
#and then convert to fold changes
ndata=data[,-1]
scaling=rowMeans(ndata[,seq(2,3)])
ndata=log(ndata / scaling)
ndata=as.data.frame(t(ndata))
#Take averages of each group
colnames(ndata)=data[,1]
pdata=matrix(ncol=ncol(ndata), nrow=nrow(ndata)/3)
for (i in 1:nrow(pdata)) {
	pdata[i,]=log(colMeans(exp(ndata[seq(3*(i-1)+1, 3*i),])))
}

colnames(pdata)=colnames(ndata)
rownames(pdata)=c("Non", "Mock D1", "Mock D14", "Mock D28", "Fat D1", "Fat D4",
	"Fat D7", "Fat D14", "Fat D28", "Cap D1", "Cap D4", "Cap D7", "Cap D14", 
	"Cap D28")
rows=ncol(ndata)

#Use sums over capsules to sort.
sums=apply(pdata[10:14,] ,2,sum)

#Make a chart with all of the groups, ordered
pdf("cust_pooled_ord.pdf",width=7, height=rows*8/75+2)
heatmap.2(t(data.matrix(pdata[,order(sums,decreasing=TRUE)])), scale="none", Colv=FALSE, dendrogram="none",
	col=greenred(32),density.info="none",trace="none", cexRow=0.7,
	cexCol=1, ylab="log Expression", symbreaks=TRUE, 
	keysize=1,margins=c(7,7), Rowv=FALSE, lhei=c(1, rows*8/75),
	main="Custom Panel")
dev.off()

#Load subpopulation data
type_file=scan("cust_pops.csv", what="", sep="\n")
types=strsplit(type_file, ',+')
names(types)=sapply(types, '[[',1)
types=lapply(types, '[', -1)

#Make plots for each type
col_breaks=seq(-4.5,4.5,length.out=33)
for(i in 1:length(types)) {
	tdata=pdata[,types[[i]]]
	sums=apply(tdata,2,sum)
	pdf(file=paste("cust_ord_", names(types)[i], ".pdf", sep=""),width=7, height=rows*8/75+2)
	heatmap.2(t(data.matrix(tdata[,order(sums,decreasing=TRUE)])), scale="none", Colv=FALSE, dendrogram="none",
		col=greenred(32),density.info="none",trace="none", cexRow=0.7,
		cexCol=1, ylab="log Expression", breaks=col_breaks,
		keysize=1,margins=c(7,7), Rowv=FALSE, lhei=c(1, rows*8/75),
		main=paste("Custom Panel: ", names(types)[i], sep=""))
	dev.off()
}



