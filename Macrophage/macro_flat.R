setwd("C:/Users/tamhok/Documents/MIT/Anderson/Fibrosis/Macrophage");
library('NanoStringNorm')
library('gplots')
# data=read.csv("macro_raw.csv", header=TRUE, stringsAsFactors=FALSE);

# raw=NanoStringNorm(data, 
			# anno=NA, 
			# CodeCount = 'geo.mean',
			# SampleContent='housekeeping.sum', 
			# Background='mean.2sd',
			# return.matrix.of.endogenous.probes=TRUE,
			# );

raw=read.csv("macro_comp.csv", header=TRUE, row.names=1)			
#To not divide by zero
low_data=t(raw < 0.0001)
ndata=raw+0.0001;

source("../pheatmap.R")
minVal=0;
scaling=rowMeans(ndata[,4:6])

#zero_data=apply(ndata,1,function(x) min(x) <= minVal)

ndata=log(ndata / scaling)
ndata=as.data.frame(t(ndata))
outc=ordered(as.character(c(rep("NT",3), rep("Mok",3), rep(1,3),
	rep(4,3), rep(7,3),rep(14,3), rep(28,3))))
pval2=apply(ndata, 2, function(x) summary(aov(x~outc))[[1]]$"Pr(>F)"[1])
pdata=matrix(ncol=ncol(ndata), nrow=7)
lowdata=matrix(ncol=ncol(ndata), nrow=7)
thesizes=c(3,3,3,3,3,3,3);
curSize=0;
for (i in 1:7) {
	pdata[i,]=log(colMeans(exp(ndata[seq(curSize+1, curSize+thesizes[i]),])))
	lowdata[i,]=colSums(low_data[seq(curSize+1, curSize+thesizes[i]),])
	curSize=curSize+thesizes[i];
}


colnames(pdata)=colnames(ndata)
rownames(pdata)=c("NT", "MOCK", "1", "4", "7", "14", "28")
adata=pdata[, pval2 < 0.01/ncol(ndata)]
rows=ncol(adata)
#Remove data that is too low
#pdata=pdata[,!zero_data]

sizes=c(1,4,7,14,28);
diff=apply(pdata,2,function (x) lm(x[c(3:7)] ~ sizes)$coefficients[2])
sums=apply(pdata,2,sum)
adata=pdata[,sums > 2 | sums < -2]
asums=sums[sums > 2 | sums < -2]

pdata=pdata[,order(sums, decreasing=TRUE)]
lowdata=lowdata[,order(sums, decreasing=TRUE)]
zerodata=colSums(lowdata==0)<2;
pdata=pdata[,!zerodata]
pdata[lowdata > 0.1] = -4

draw_map(pdata[-1,], "macro_flat_Mock.tiff")

brk_len=ceiling(ncol(pdata)/3)
infl_breaks=c(0,brk_len, 2*brk_len, ncol(pdata))

panel_sizes=vector()
for(i in 1:3) {
	tdata=pdata[, (infl_breaks[i]+1):infl_breaks[i+1]]
	panel_size=get_size(t(data.matrix(tdata)), filename="test.pdf", showcolns=TRUE)
	panel_sizes=c(panel_sizes, panel_size)
}

min_width=sum(panel_sizes[c(1,3,5)])
min_height=max(panel_sizes[c(2,4,6)])
print(c(min_width, min_height))
tiff("infl_panel.tiff", units="in", width=min_width+5, height=min_height+2, res=300)

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



			
