setwd("C:/Users/tamhok/Documents/MIT/Anderson/Fibrosis/Macrophage/macro_new");
source("macro_functions.R")
nano_data=load_data("macro_raw_all", 3);
proc_data=average_data(nano_data$raw, nano_data$class, nano_data$limit); 
#Determine which scalings to use
scaling=vector(length=nrow(proc_data$pdata))
labs=rownames(proc_data$pdata);
ids=sub(".*_", "", labs)
scale_inds=c(which(labs=="C_C_IP"), which(labs=="C_C_Fat"), which(labs=="C_C_Cap"))
scaled_data=proc_data$pdata[scale_inds,]
#scaled_data[scaled_data < 1]=1;
#scaled_data=rbind(scaled_data, rep(1, ncol(scaled_data)));
scaling[which(ids=="Cap")]=3
scaling[which(ids=="IP")]=1
scaling[which(ids=="Fat")]=2
kdata=proc_data$pdata[scale_inds,]
sdata=scale_data(proc_data$pdata, scaling, scaled_data)
#sdata=scale_data(proc_data$pdata, scaling, mean(nano_data$limit), scaling_mat=scaled_data)
kdata[kdata<0.001] = -2000
kdata[kdata > 0.001] = 0
sdata[scale_inds,]=kdata
idx=which(ids=="Cap" & labs != "C_C_Cap")
#sdata[idx,]=log((exp(sdata[idx,]) + proc_data$negs[idx])/proc_data$negs[idx])
sdata=sdata[order(scaling),]
rdata=remove_zero_data(sdata, proc_data$zdata)
ip_d=rdata[which(ids=="IP"),];
cap_d=rdata[which(ids=="Cap"),];
fat_d=rdata[which(ids=="Fat"),];

#draw_map(ip_d, "macro_ip.tiff")
#draw_map(cap_d, "macro_cap.tiff")
#draw_map(fat_d, "macro_fat.tiff")

#Import groups data
gdata=read.csv("../macro_groups.csv", header=TRUE, stringsAsFactors=FALSE);
grps=gdata[which(gdata[,2]!="None"),]

#Segment data further
u_cls=t(unique(t(nano_data$class)))
classes=rbind(u_cls, create_class_names(u_cls))
colnames(classes)=classes[4,]
classes=classes[, rownames(rdata)]
source("../macro_pie.R")

#obtain series and make pie charts
spaces=as.factor(classes[3,]);
sizes=as.factor(classes[1,]);
spacesl=levels(spaces);
sizesl=levels(sizes);

plotctr=0;
sums=apply(rdata,2,sum)
for (i in 1:length(spacesl)) {
	postscript(paste("macro_pie", spacesl[4-i], ".eps", sep=""), width=15, height=10)
	par(mfrow=c(2, 3))
	for (j in 1:2) {
		indices=which(spaces==spacesl[4-i] & sizes==sizesl[j])
		if(length(indices) > 0) {
			cdata=rdata[indices,]
			grp_data=obtain_group_data(cdata, grps, sums)
			for (k in 1:length(indices)) {
				make_pie(grp_data[,c(1:3, k+3, ncol(grp_data))], title=paste(spacesl[4-i], sizesl[j], "mm Day", classes[2, indices[k]]))
			}
			if( k < 3) {
				plot.new()
			}
		}
	}
	dev.off()
	#tiff(paste("macro_flat_all", spacesl[4-i], ".tiff", sep=""), units="in", width=15, height=5.5, res=600)
	#make_flats(rdata[which(spaces==spacesl[4-i] & sizes %in% sizesl[1:2]),], grps)
	#dev.off()
}
postscript("macro_pie_legend.eps", width=6, height=6)
tdata=cdata[1,]
tdata[1:length(tdata)]=-10
tgrp_data=obtain_group_data(tdata, grps, sums)
make_pie(tgrp_data, title="Legend")
dev.off()

cobj=hclust(dist(t(rdata)))
rdata=rdata[,cobj$order]

cobq=hclust(dist(t(sdata)))
sdata=sdata[,cobj$order]

sums=apply(rdata,2,sum)
gdata=obtain_group_data(rdata,grps,sums)
types=as.factor(gdata[,2])
lvls=levels(types)
for (i in 1:3) {
	nom=gdata[types==lvls[i],1]
	td=data.matrix(gdata[types==lvls[i],4:(ncol(gdata)-1)])
	#windows()
	#heatmap.2(td, labRow=nom, col=greenred(62), dendrogram="row", Colv=NULL, breaks=cbreaks, trace="none")
}

grps=load_groups("../macro_groups")
new_names=do_stats(nano_data, grps, filt_word="None")
do_stat_alt(nano_data, grps, filt_word = "None")
grps=load_groups("../macro_groups_new")

comp_order=c("IP", "Fat", "Cap")
size_order=c("500", "1500")
clsnames=classes[4,]
index=0;

cdata=rbind(nano_data$class, nano_data$raw)

take_anova=function(data) {
  dats=log(as.numeric(data))
  newdf=data.frame(days,dats)
  return(anova(aov(dats~days, newdf))$'Pr(>F)'[[1]])
}

for(i in 1:length(comp_order)) {
  for (j in 1:length(size_order)) {
    temp_data=cdata[c(-1,-3),which(cdata[1,]==size_order[j] & cdata[3,]==comp_order[i])]
    days=as.vector(t(temp_data[1,]))
    
  }
}

#build class data
classdata=matrix(nrow=0,ncol=3);
for(i in 1:length(comp_order)) {
  cur_set=classes[,which(classes[3,]==comp_order[i])]
  for (j in 1:length(size_order)) {
    sub_set=cur_set[, which(cur_set[1,]==size_order[j])]
    sub_set=sub_set[c(2,4), order(sub_set[2,])]
    C_set=cur_set[c(1,4), which(cur_set[1,]=="C"), drop=FALSE]
    if(length(C_set) > 0) {
      sub_set=cbind(C_set, sub_set);
    }
    index=index+1;
    sub_set=t(rbind(sub_set, rep(index, ncol(sub_set))));
    classdata=rbind(classdata, sub_set)
  }
}
cbreaks=seq(-4,4,length.out=63)
col_breaks=c(-15, cbreaks, 15)
colscale1=list(cbreaks=cbreaks, col_breaks=col_breaks, legend=c(-4,-2,0,2,4), scale=greenred(64))
cbreaks=seq(0,10,length.out=63)
col_breaks=c(-15, cbreaks, 30)
colscale2=list(cbreaks=cbreaks, col_breaks=col_breaks, legend=c(0,2,4,6,8,10), scale=colorpanel(64, low="#000000", high="#FF0000"))
colscale=list(colscale1, colscale1, colscale2)

splitdata=split_by_class(sdata, classdata)

playout=rbind(c(2,4,6), c(1,3,5))
rownames(playout)=c("1500", "500")
colnames(playout)=c("IP", "Fat", "Cap")
#playout=t(playout)
new_names=do_stats(nano_data, grps)
make_grouped_plot("cell_types.eps", splitdata, grps, playout, colscale, byrow=2, labrow=3, opt_names=new_names, cluster=FALSE)
#allgrp=rbind(grps[1,], rep("All", nrow(grps)), grps[3,])
#make_grouped_plot("everything.eps", splitdata, allgrp, playout, colscale, byrow=2, labrow=3)


