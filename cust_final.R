library(gplots)
source("pheatmap.R")
minVal=10
#Load data
data=read.csv(file="cap_custom_neg.csv", header=TRUE, stringsAsFactors=FALSE)
ndata=data[,-1]

gene_names=read.csv(file="gene_names.csv", header=FALSE,stringsAsFactors=FALSE)
rownames(gene_names)=gene_names[,1]
cbreaks=seq(-4,4,length.out=63)
col_breaks=c(-10, cbreaks, 10)

draw_map=function(data, filename=NA, sidelabel=NA) {
	cbreaks=seq(-4,4,length.out=63)
	col_breaks=c(-100, cbreaks, 100)
	if(!is.na(filename)) {
		width=nrow(data)*36/72+4
		tiff(filename, width=width, units="in", res=300, height=ncol(data)*20/72+1)
	}
	pheatmap(t(data.matrix(data)), col=greenred(64), cluster_rows=FALSE, cluster_cols=FALSE, legend_breaks=c(-4,-2,0,2,4), legend_labels=c(-4,-2,0,2,4),
		family="Sans", fontsize=18, cellwidth=36, cellheight=20, main=sidelabel,
		breaks=col_breaks, rbreaks=cbreaks, scale="none", legend=TRUE)
	if(!is.na(filename)) {
		dev.off()
	}
}

get_size=function(data, sidelabel=NA, filename, showcolns) {
	dim=pheatmap(data, col=greenred(64), cluster_rows=FALSE, 
		cluster_cols=FALSE, family="Sans", fontsize=18, cellwidth=36,
		cellheight=20, main=sidelabel, filename=filename, show_colnames=showcolns,
		breaks=col_breaks, scale="none", legend=FALSE)
	return(dim)
}

#Scale to the mock data
scaling=rowMeans(ndata[,seq(1,3)])
zero_data=apply(ndata,1,function(x) min(x) < minVal)
#zero_data=scaling < minVal & rowMeans(ndata) < minVal;
ndata=log(ndata / scaling)
ndata=as.data.frame(t(ndata))
colnames(ndata)=gene_names[data[,1],2]
ndata[t(data[,-1]) < 1]=-1000
pdata=matrix(ncol=ncol(ndata), nrow=nrow(ndata)/3)
for (i in 1:6) {
	pdata[i,]=log(colMeans(exp(ndata[seq(3*(i-1)+1, 3*i),])))
}
colnames(pdata)=colnames(ndata)
rownames(pdata)=c("Mock", "Day 1", "Day 4", "Day 7", "Day 14", "Day 28")
rows=ncol(ndata)
pdata=pdata[,!zero_data]
zn=colnames(ndata)[zero_data]
sums=apply(pdata,2,sum)
#pdata[,zero_data]=matrix(data=-1000, nrow=nrow(pdata), ncol=ncol(pdata[,zero_data]))
# ddata=pdata;
# for (i in 2:6) {
	# ddata[i,]=pdata[i,]-pdata[(i-1),]
# }
# draw_map(ddata[,order(sums,decreasing=TRUE)], filename="cust_pooled_diffs.tiff")
# rownames(ddata)=rownames(pdata)

updata=apply(pdata,2,max)
downdata=apply(pdata,2,min)
updata_all=updata
downdata_all=downdata
print(length(downdata_all))
updata=updata[updata > sd(pdata)]
updata=updata[order(-updata)]
downdata=downdata[downdata < -1 * sd(pdata)]
downdata=downdata[order(downdata)]
write.table(updata, file="up_out_cust.txt", append=FALSE, quote=FALSE, sep="\t", col.names=FALSE)
write.table(downdata, file="down_out_cust.txt", append=FALSE, quote=FALSE, sep="\t", col.names=FALSE)
# for (i in 2:6) {
	
	# wdata=ddata[i, ddata[i,] > 1.5 | ddata[i,] < -1.5]
	# wdata=wdata[order(wdata,  decreasing=TRUE)]
	# write.table(wdata, file="up_down_out.txt", append=TRUE, col.names=rownames(pdata)[i], quote=FALSE, sep="\t")
# }
#Make plots and the panel
draw_map(pdata[,order(sums,decreasing=TRUE)], filename="cust_pooled_ord.tiff")

type_file=scan("cust_pops.csv", what="", sep="\n")
types=strsplit(type_file, ',+')
names(types)=sapply(types, '[[',1)
types=lapply(types, '[', -1)

panel_sizes=vector();
tdatas=list();
panels=c("MACRO", "NEUTRO", "B", "ANGIO")
pan_names=c("Macrophage", "Neutrophil", "B Cell", "Angiogenesis")

showcolns=c(T,T,F,F);
for(i in 1:length(panels)) {
	elts=setdiff(gene_names[types[[panels[i]]],2],zn)
	tdata=pdata[,elts]
	sums=apply(tdata,2,sum)
	tdatas[[i]]=t(data.matrix(tdata[,order(sums, decreasing=TRUE)]))
	panel_size=get_size(tdatas[[i]], pan_names[i], paste("cust_ord_pooled", panels[i], ".pdf", sep=""), showcolns[i])
	panel_sizes=c(panel_sizes, panel_size)
}
names(panel_sizes)=panels;
names(tdatas)=panels;

min_width=panel_sizes[1] + max(panel_sizes[c(3,5,7)])
min_height=max(panel_sizes[2], sum(panel_sizes[c(4,6,8)]))

tiff("cust_panel.tiff", units="in", width=min_width+2, height=min_height+2, res=300)

p1_width=unit(panel_sizes[1], "inches")
break_width=unit(0.5, "inches")
p2_width=unit(max(panel_sizes[c(3,5,7)]), "inches")
print(panel_sizes)
p1_height=unit(panel_sizes[4], "inches")
p2_height=unit(panel_sizes[6], "inches")
p3_height=unit(panel_sizes[8], "inches")
vbrk_ht=unit((panel_sizes[2] - sum(panel_sizes[c(4,6,8)]))/2, "inches")
widths=unit.c(p1_width, p2_width, break_width)
heights=unit.c(p1_height, vbrk_ht, p2_height, vbrk_ht, p3_height)

pushViewport(viewport(layout = grid.layout(nrow = 5, ncol = 3, widths=widths, heights=heights), gp = gpar(fontsize=18)))

pushViewport(vplayout(1:5, 1,just=c("left", "top")))
pheatmap(tdatas[[1]], col=greenred(64), cluster_rows=FALSE, 
		cluster_cols=FALSE, family="Sans", fontsize=18, cellwidth=36,
		cellheight=20, main=pan_names[1], show_colnames=TRUE,
		breaks=col_breaks, scale="none", legend=FALSE)
upViewport()
for(i in 2:4) {
pushViewport(vplayout(((i-1)-1)*2+1, 2, just=c("left", "top")))
pheatmap(tdatas[[i]], col=greenred(64), cluster_rows=FALSE, 
		cluster_cols=FALSE, family="Sans", fontsize=18, cellwidth=36,
		cellheight=20, main=pan_names[i], show_colnames=showcolns[i],
		breaks=col_breaks, scale="none", legend=FALSE)
upViewport()
}

pushViewport(vplayout(1:2,3, just="bottom"))
legend = c(-4,-2,0 ,2,4)
names(legend)=legend
draw_legend(color=greenred(64), ht = 0.8, cbreaks, legend, fontsize = 18, family="Sans")
upViewport()

dev.off()
dev.off()

data=read.csv(file="infl_neg.csv", header=TRUE)
ndata=data[,-1]
nsampl=4
scaling=rowMeans(ndata[,seq(1,nsampl)])
zero_data=apply(ndata,1,function(x) min(x) < minVal)
#zero_data=scaling < minVal & rowMeans(ndata) < minVal;
ndata=log(ndata / scaling)
ndata=as.data.frame(t(ndata))
ndata[t(data[,-1]) < minVal]=-1000
colnames(ndata)=data[,1]
outc=ordered(as.character(c(rep("Mok",nsampl), rep(1,nsampl), 
	rep(4,nsampl), rep(7,nsampl),rep(14,nsampl), rep(28,nsampl))))
pval2=apply(ndata, 2, function(x) summary(aov(x~outc))[[1]]$"Pr(>F)"[1])
pdata=matrix(ncol=ncol(ndata), nrow=nrow(ndata)/nsampl)
for (i in 1:6) {
	pdata[i,]=log(colMeans(exp(ndata[seq(nsampl*(i-1)+1, nsampl*i),])))
}

colnames(pdata)=colnames(ndata)
rownames(pdata)=c("Day 0", "Day 1", "Day 4", "Day 7", "Day 14", "Day 28")
adata=pdata[, pval2 < 0.01/ncol(ndata)]
rows=ncol(adata)
#Remove data that is too low
pdata=pdata[,!zero_data]
#pdata[,zero_data]=matrix(data=-1000, nrow=nrow(pdata), ncol=ncol(pdata[,zero_data]))
sums=apply(pdata,2,sum)
adata=pdata[,sums > 2 | sums < -2]
asums=sums[sums > 2 | sums < -2]

updata=apply(pdata,2,max)
downdata=apply(pdata,2,min)
updata_all=c(updata_all, updata)
downdata_all=c(downdata_all, downdata)
print(length(downdata_all))
updata=updata[updata > sd(pdata)]
updata=updata[order(-updata)]
downdata=downdata[downdata < -1 * sd(pdata)]
downdata=downdata[order(downdata)]
write.table(updata, file="up_out_infl.txt", append=FALSE, quote=FALSE, sep="\t", col.names=FALSE)
write.table(downdata, file="down_out_infl.txt", append=FALSE, quote=FALSE, sep="\t", col.names=FALSE)

draw_map(adata[,order(asums, decreasing=TRUE)], "infl_min_max.tiff")

brk_len=ceiling(ncol(pdata)/3)
infl_breaks=c(0,brk_len, 2*brk_len, ncol(pdata))

pdata=pdata[,order(sums, decreasing=TRUE)]
panel_sizes=vector()
for(i in 1:3) {
	tdata=pdata[, (infl_breaks[i]+1):infl_breaks[i+1]]
	panel_size=get_size(t(data.matrix(tdata)), filename=paste("cust_ord_pooled", panels[i], ".pdf", sep=""), showcolns=TRUE)
	panel_sizes=c(panel_sizes, panel_size)
}

min_width=sum(panel_sizes[c(1,3,5)])
min_height=max(panel_sizes[c(2,4,6)])
print(c(min_width, min_height))
tiff("infl_panel.tiff", units="in", width=min_width+5, height=min_height+2, res=600)

grid.newpage()

p1_width=unit(panel_sizes[1], "inches")
p2_width=unit(panel_sizes[3], "inches")
p3_width=unit(panel_sizes[5], "inches")
break_width=unit(0.4, "inches")
p1_height=unit(min_height, "inches")
widths=unit.c(p1_width, break_width, p2_width, break_width, p3_width, break_width)
heights=p1_height

pushViewport(viewport(layout = grid.layout(nrow = 1, ncol = 6, widths=widths, heights=heights), gp = gpar(fontsize=18)))

for(i in 1:3) {
tdata=pdata[, (infl_breaks[i]+1):infl_breaks[i+1]]
pushViewport(vplayout(1, (i-1)*2+1, just=c("left", "top")))
pheatmap(t(data.matrix(tdata)), col=greenred(64), cluster_rows=FALSE, 
		cluster_cols=FALSE, family="Sans", fontsize=18, cellwidth=36,
		cellheight=20, show_colnames=TRUE,
		breaks=col_breaks, scale="none", legend=FALSE)
upViewport()
}

pushViewport(vplayout(1,6, just="bottom"))
legend = c(-4,-2,0 ,2,4)
names(legend)=legend
draw_legend(color=greenred(64), ht=0.95, cbreaks, legend, fontsize = 18, family="Sans")
upViewport()

dev.off()
upViewport()

data=read.csv(file="miRNA_comb.csv", header=TRUE)
ndata=data[,-1]
nsampl=2
scaling=rowMeans(ndata[,seq(1,nsampl)])
zero_data=apply(ndata,1,function(x) min(x) < minVal)
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
pdata=pdata[,!zero_data]

updata=apply(pdata,2,max)
downdata=apply(pdata,2,min)
updata_all=c(updata_all, updata)
downdata_all=c(downdata_all, downdata)
print(length(downdata_all))
updata=updata[updata > sd(pdata)]
updata=updata[order(-updata)]
downdata=downdata[downdata < -1 * sd(pdata)]
downdata=downdata[order(downdata)]
write.table(updata, file="up_out_miRNA.txt", append=FALSE, quote=FALSE, sep="\t", col.names=FALSE)
write.table(downdata, file="down_out_miRNA.txt", append=FALSE, quote=FALSE, sep="\t", col.names=FALSE)

updata=updata_all[updata_all>2*sd(updata_all)]
updata=updata[order(-updata)]
write.table(updata, file="up_out.txt", append=FALSE, quote=FALSE, sep="\t", row.names=names(updata), col.names=FALSE)
print(length(downdata_all))
print(sd(downdata_all))
downdata=downdata_all[downdata_all< -2*sd(downdata_all)]
print(length(downdata_all))
downdata=downdata[order(downdata)]
alldata=c(updata, downdata)
alldata=alldata[order(-alldata)]
write.table(downdata, file="down_out.txt", append=FALSE, quote=FALSE, sep="\t", row.names=names(downdata), col.names=FALSE)
write.table(alldata, file="all_out.txt", append=FALSE, quote=FALSE, sep="\t", row.names=names(alldata), col.names=FALSE)