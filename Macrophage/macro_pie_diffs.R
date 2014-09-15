source("macro_pie.R")

gdata=read.csv("macro_groups.csv", header=TRUE, stringsAsFactors=FALSE);
grps=gdata[which(gdata[,2]!="None"),]
raw=read.csv("macro_proc.csv", header=TRUE, row.names=1)

to_binary=function(val, high=4, low=-4, thresh=0.25) {
	im=val;
	im[im > thresh] = high;
	im[im > -thresh & im < thresh] = 0
	im[im < -thresh] = low;
	return(im);
}

#To not divide by zero
low_data=t(raw < 0.0001)
ndata=raw+0.0001;

minVal=0;
scaling=rowMeans(ndata[,7:9])

#make groups
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

sizes=c(1,4,7,14,28);
diffs=apply(pdata,2,function (x) lm(x[c(3:7)] ~ sizes)$coefficients[2])
ddata=diff(pdata); 
pdata=rbind(ddata[2,], pdata[3,], ddata[3:6,])
#pdata=pdata[2:7,]
pdata=to_binary(pdata);
zerodata=colSums(lowdata==0)<2;
#pdata[lowdata[2,7] > 0.1] = -3.95;
pdata=pdata[,!zerodata]
sums=apply(pdata[2:6,],2,function(x) 2*sum(x==4)+sum(x==0))
diffs=diffs[!zerodata]
grp.in=which(grps[,1] %in% colnames(pdata))
grps=grps[grp.in,];
udata=rbind(pdata,sums)
udata=t(udata[,grps[,1]])
grps=cbind(grps, udata);

tiff("macro_pie_diffs.tiff", units="in", width=25, height=5, res=300)
titles=c("MOK", "Day 1", "Day 4", "Day 7", "Day 14", "Day 28");
par(mfrow=c(1, 6))
for (i in 1:6) {
	make_pie(grps[,c(1:3, i+3, ncol(grps))], title=titles[i])
}
dev.off()
