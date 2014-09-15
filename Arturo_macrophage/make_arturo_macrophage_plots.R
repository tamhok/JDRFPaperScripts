source("C:/Users/tamhok/Documents/MIT/Anderson/Fibrosis/Macrophage/macro_new/macro_functions.R");
setwd("C:/Users/tamhok/Documents/MIT/Anderson/Fibrosis/Arturo_macrophage/")
nano_data=load_data("art_macro_raw_all", 2);
#nano_data=new_nano_data;
proc_data=average_data(nano_data$raw, nano_data$class, nano_data$limit);

#Determine which scalings to use
scaling=vector(length=nrow(proc_data$pdata))
labs=rownames(proc_data$pdata);
ids=sub("_.*", "", labs)
scale_inds=c(which(labs=="IP_Mock"), which(labs=="Fat_Mock"), which(labs=="Cap_SLG20"))
scaled_data=proc_data$pdata[scale_inds,]
#scaled_data[scaled_data < 1]=1;
scaled_data=rbind(scaled_data, rep(1, ncol(scaled_data)));
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

u_cls=t(unique(t(nano_data$class)))
classes=rbind(u_cls, create_class_names(u_cls))
colnames(classes)=classes[3,]
classes=classes[, rownames(rdata)]

comp_order=c("IP", "Fat", "Cap")
mat_order=c("Mock","E9", "VLVG/SLG100", "SLG20")
cap_mat_order = c("E9", "VLVG/SLG100", "SLG20")
#mat_order=c("Mock","E9", "VLVG/SLG100", "SLG20", "7")

cobj=hclust(dist(t(rdata)))
rdata=rdata[,cobj$order]

cobq=hclust(dist(t(sdata)))
sdata=sdata[,cobj$order]

#build class data
classdata=matrix(nrow=0,ncol=3);
for(i in 1:length(comp_order)) {
  cur_set=classes[,which(classes[1,]==comp_order[i])]
  colnames(cur_set)=cur_set[2,]
  if(comp_order[i] == "Cap") {
    cur_set=cur_set[,cap_mat_order]
  } else {
    cur_set=cur_set[,mat_order]
  }
  cur_set=cbind(t(cur_set)[,2:3], rep(i, dim(cur_set)[2]))
  classdata=rbind(classdata, cur_set)
}

splitdata=split_by_class(rdata, classdata)
splitdata2=split_by_class(sdata, classdata)
e9data=rdata[which(classes[2,]=="E9"),]
slg20data=rdata[which(classes[2,]=="SLG20"),]
diffdata=slg20data-e9data
ranks2=apply(diffdata, 2, sum)

#Compute the ones that change the most, rearrange to order by that metric
for (i in 1:length(splitdata)) {
  cur_set=splitdata[[i]]$data
  diffs=cur_set["E9",]-cur_set["SLG20",]
  cur_set=cur_set[,order(ranks2)]
  splitdata[[i]]$data=cur_set
}
#Import groups data for pies
grps=load_groups("art_macro_groups_wheels")
grps=grps[,-which(grps[2,]=="")]
grps=t(grps)
source("../Macrophage/macro_pie.R")

#obtain series and make pie charts
sums=apply(rdata,2,sum)
for (i in 1:length(comp_order)) {
  postscript(paste("macro_pie", comp_order[i], ".eps", sep=""), width=20, height=10, paper="special")
  par(mfrow=c(1, 4))
  cdata=splitdata[[i]]$data
  grp_data=obtain_group_data(data.matrix(cdata), grps, sums)
  for (k in 1:length(mat_order)) {
    grp_txt=grp_data[,c(1:3)]
    grp_num=grp_data[,c(k+3,ncol(grp_data))]
    class(grp_num)="numeric"
    grp_combo=data.frame(grp_txt, grp_num, stringsAsFactors=FALSE)
    make_pie(grp_combo, title=mat_order[k])
  }
  
  dev.off()
  #tiff(paste("macro_flat_all", spacesl[4-i], ".tiff", sep=""), units="in", width=15, height=5.5, res=600)
  #make_flats(rdata[which(spaces==spacesl[4-i] & sizes %in% sizesl[1:2]),], grps)
  #dev.off()
}
postscript("macro_pie_legend.eps", width=6, height=6)
tdata=cdata[1,]
tdata[1:length(tdata)]=-10
grp_data=obtain_group_data(tdata, grps, sums)
grp_txt=grp_data[,c(1:3)]
grp_num=grp_data[,c(4,ncol(grp_data))]
class(grp_num)="numeric"
grp_combo=data.frame(grp_txt, grp_num, stringsAsFactors=FALSE)
make_pie(grp_combo, title="Legend")
dev.off()

#do stats
raw_mat=nano_data$raw
cnames=create_class_names(nano_data$class)
colnames(raw_mat)=cnames
raw_data=melt(raw_mat)
raw_grps=t(rbind(apply(raw_data, 1, function(x) unlist(strsplit(x[2], "_")))))
vals=raw_data[,3]
is_zero = vals==0
vals=log(raw_data[,3]+min(vals[vals > 0]) / 2)

fdata=data.frame(vals=vals, genes=raw_data[,1], mat=raw_grps[,2], comp=raw_grps[,1], is_zero=is_zero)
#f_filt=fdata[-which(fdata$comp=="Cap"),]
f_filt=fdata
model=aov(vals ~ mat * comp + genes + is_zero, data=f_filt)
comp_order = c("Fat", "IP", "Cap")
mat_order = c("E9", "SLG20","VLVG/SLG100")
gene_order = levels(as.factor(as.character((f_filt$genes))))

apply_tukey = function(i,j, comparator="Mock") {
  cur_set=f_filt[which(f_filt$comp==comp_order[i] & f_filt$genes==j),]
  model=aov(vals ~ mat, data=cur_set)
  res = TukeyHSD(model)$mat
  labels = c("SLG20-E9", "VLVG/SLG100-E9")
  if(res[labels[1],1] < 0) res[labels[1],4]=res[labels[1],4]+1
  if(res[labels[2],1] < 0) res[labels[2],4]=res[labels[2],4]+1
  return(res[labels, 4])
}

rows_nans = function(x) {
  unlist(apply(x, 1, function(x) {sum(is.na(x)) > 0}))
}
comp_order = c("Fat", "IP", "Cap")
compare_order=c("Mock", "Mock", "E9")
for (i in 1:length(comp_order)) {
  print(comp_order[i])
  combdvals=lapply(gene_order, function(j) apply_tukey(i, j, compare_order[i]))
  combdvals = Reduce(rbind, combdvals)
  adjvals = apply(combdvals[which(!rows_nans(combdvals)),], 2, function(x) p.adjust(x, "fdr"))
  diffvals = apply(adjvals < 0.05, 2, sum)
  print(diffvals)
}

#Load all groups
grps=load_groups("art_macro_groups")
allgrp=load_groups("eveything_groups")
grpname=allgrp[2,]
grpname[grpname != ""] = " ";
allgrp[2,]=grpname
#Plotting information
cbreaks=seq(-4,4,length.out=63)
col_breaks=c(-15, cbreaks, 15)
colscale1=list(cbreaks=cbreaks, col_breaks=col_breaks, legend=c(-4,-2,0,2,4), scale=greenred(64))
cbreaks=seq(0,4,length.out=63)
col_breaks=c(-15, cbreaks, 30)
colscale2=list(cbreaks=cbreaks, col_breaks=col_breaks, legend=c(0,2,4,6,8,10), scale=colorpanel(64, low="#000000", high="#FF0000"))
colscale=list(colscale1, colscale1, colscale1)

playout=t(as.matrix(c(1,2,3)))
rownames(playout)=c("")
colnames(playout)=c("IP", "Fat", "Cap")
make_grouped_plot("cell_types_arturo.eps", splitdata, grps, playout, colscale, byrow=2, labrow=3)#, opt_names=new_names)
make_grouped_plot("cell_types_all_arturo.eps", splitdata2, grps, playout, colscale, byrow=2, labrow=3)#, opt_names=new_names)
make_grouped_plot("everything_arturo.eps", splitdata2, allgrp, playout, colscale, byrow=2, labrow=3)
