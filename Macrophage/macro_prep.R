setwd("C:/Users/tamhok/Documents/MIT/Anderson/Fibrosis/Macrophage");
library('NanoStringNorm')
library('ggplot2')
library('plotrix')
library('gplots')
library('gridBase')
source("macro_pie.R")
source("C:/Users/tamhok/Documents/MIT/Anderson/Fibrosis/pheatmap.R")

# data=read.csv("macro_raw.csv", header=TRUE, stringsAsFactors=FALSE);
# data[data$Name %in% 
	# c("Cltc", "Gapdh", "Gusb", "Hprt1", "Pgk1", "Tubb5"),
	 # 'Code.Class'] = 'Housekeeping';

# classes=read.csv("macro_raw_groups.csv", header=TRUE, stringsAsFactors=FALSE);
	 
# raw=NanoStringNorm(data, 
			# anno=NA, 
			# CodeCount = 'geo.mean',
			# SampleContent='housekeeping.sum', 
			# Background='mean.2sd',
			# return.matrix.of.endogenous.probes=TRUE,
			# );

# write.csv(raw, "macro_proc.csv")
gdata=read.csv("macro_groups.csv", header=TRUE, stringsAsFactors=FALSE);
grps=gdata[which(gdata[,2]!="None"),]
raw=read.csv("macro_proc.csv", header=TRUE, row.names=1)

#To not divide by zero
low_data=t(raw < 0.0001)
ndata=raw+0.0001;

minVal=0;
scaling=rowMeans(ndata[,4:6])

#make groups
ndata=log(ndata / scaling)
ndata=as.data.frame(t(ndata))
outc=ordered(as.character(c(rep("NT",3), rep("Mock",3), rep(1,3),
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

sizes=c(1,4,7,14,28);
diffs=apply(pdata,2,function (x) lm(x[c(3:7)] ~ sizes)$coefficients[2])
sums=apply(pdata,2,sum)
adata=pdata[,sums > 2 | sums < -2]
asums=sums[sums > 2 | sums < -2]
zerodata=colSums(lowdata==0)<2;
pdata[lowdata > 0.1] = -3.95;
pdata=pdata[,!zerodata]
sums=sums[!zerodata]
diffs=diffs[!zerodata]
grp.in=which(grps[,1] %in% colnames(pdata))
grps=grps[grp.in,];
udata=rbind(pdata, sums)
udata=t(udata[,grps[,1]])
grps=cbind(grps, udata);


postscript("macro_pie.eps", horizontal = FALSE, onefile = FALSE, paper = "special", width=25, height=5)
par(mfrow=c(1, 7))
for (i in 2:7) {
	make_pie(grps[,c(1:3, i+3, 11)], title=outc[3*i-2])
}
plot.new()
vps=baseViewports()
pushViewport(vps$figure)
cbreaks=seq(-4,4,length.out=63)
legend=c(-4,-2,0,2,4)
names(legend)=legend
draw_legend(color=greenred(64), ht=0.95, cbreaks, legend, fontsize = 18, family="Sans")

dev.off()
