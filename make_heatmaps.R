#R script to process data from the nanoString CSV file and
#a) do hierarchical clustering
#b) make a nice heatmap
#Function takes a CSV file filename, a plot title, and an array of
#numbers of columns to be used for scaling (excluding the gene names).
library(gplots)
setwd("C:/Users/tamhok/Documents/MIT/Anderson/Fibrosis")

#function to make the plots
make_heatmap_plot<-function(filename, title, scale) {
	data=read.csv(file=paste(filename, ".csv", sep=""), header=TRUE)	
	ndata=data.matrix(data[,-1])
	scaling=rowMeans(ndata[,scale])
	ndata=log(ndata / scaling)
	rows=nrow(ndata)
	pdf(file=paste(filename, ".pdf",sep=""),width=7, height=rows*8/75+2) 
	heatmap.2(ndata, scale="none", Colv=FALSE, dendrogram="row",
		col=greenred, labRow=data[,"Gene.Name"],density.info="none",
		trace="none", cexRow=0.7, cexCol=1, ylab="log Expression",
		symbreaks=TRUE,
		main=title, keysize=1,margins=c(7,7),lhei=c(2, rows*8/75));
	dev.off()
}

make_miRNA_heatmap_plot<-function(filename, title, scale) {
	data=read.csv(file=paste(filename, ".csv", sep=""), header=TRUE)	
	ndata=data.matrix(data[,seq(-3,-1)])
	scaling=data[,scale]
	ndata=log(ndata / scaling)
	rows=nrow(ndata)
	pdf(file=paste(filename, scale, ".pdf",sep=""),width=7, height=rows*8/75+2) 
	heatmap.2(ndata, scale="none", Colv=FALSE, dendrogram="none",
		col=greenred, labRow=data[,"Gene.Name"],density.info="none",
		trace="none", cexRow=0.7, cexCol=1, ylab="log Expression",
		symbreaks=TRUE, Rowv=FALSE
		main=title, keysize=1,margins=c(7,7),lhei=c(2, rows*8/75));
	dev.off()
}

make_miRNA_heatmap_plot_2<-function(filename, title) {
	data=read.csv(file=paste(filename, ".csv", sep=""), header=TRUE)	
	ndata=data.matrix(data[,seq(-3,-1)])
	ndata=log(ndata[,seq(1,5)] / ndata[,seq(6,10)])
	rows=nrow(ndata)
	pdf(file=paste(filename, ".pdf",sep=""),width=7, height=rows*8/75+2) 
	heatmap.2(ndata, scale="none", Colv=FALSE, dendrogram="row",
		col=greenred, labRow=data[,"Gene.Name"],density.info="none",
		trace="none", cexRow=0.7, cexCol=1, ylab="log Expression",
		symbreaks=TRUE,
		main=title, keysize=1,margins=c(7,7),lhei=c(2, rows*8/75));
	dev.off()
}

make_miRNA_heatmap_plot("fat_cap_miRNA", "miRNA Panel, Non control", 2)
make_miRNA_heatmap_plot("fat_cap_miRNA", "miRNA Panel, Mok control", 3)
make_miRNA_heatmap_plot_2("fat_cap_miRNA", "miRNA Panel, Fat control")
make_heatmap_plot("cap_custom_only", "Custom Panel", seq(1,3))
make_heatmap_plot("inflammation", "Inflammation Panel", seq(1,4))
