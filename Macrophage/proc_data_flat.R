source("C:/Users/tamhok/Documents/MIT/Anderson/Fibrosis/R/nanostring_prep.R")
setwd("C:/Users/tamhok/Documents/MIT/Anderson/Fibrosis/macrophage/macro_new")

nano_data=load_data("macro_raw_all", 3);
proc_data=average_data(nano_data$raw, nano_data$class, nano_data$limit); 
#Determine which scalings to use
scaling=vector(length=nrow(proc_data$pdata))
classes=proc_data$class_lst;
scaled_data=proc_data$pdata[c(which(classes[4,]=="Mock_14_IP"), which(classes[4,]=="Mock_14_Fat"), which(classes[4,]=="1500_1_Cap")),]
scaled_data[scaled_data < 1]=1;
scaled_data=rbind(scaled_data, rep(1, ncol(scaled_data)));
scaling[which(classes[3,]=="Cap")]=4
scaling[which(classes[3,]=="IP")]=1
scaling[which(classes[3,]=="Fat")]=2
sdata=scale_data(proc_data$pdata, scaling, mean(nano_data$limit), scaling_mat=scaled_data)

sdata=remove_zero_data(sdata, proc_data$zdata)
sums=apply(sdata, 2, sum)
sdata=sdata[order(classes[3,], classes[1,], classes[2,]),order(sums, decreasing=TRUE)]

#load groups
grps=load_groups("../macro_groups_all")

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
		mock_set=cur_set[c(1,4), which(cur_set[1,]=="Mock"), drop=FALSE]
		if(length(mock_set) > 0) {
			sub_set=cbind(mock_set, sub_set);
		}
		index=index+1;
		sub_set=t(rbind(sub_set, rep(index, ncol(sub_set))));
		classdata=rbind(classdata, sub_set)
	}
}
cbreaks=seq(-4,4,length.out=63)
col_breaks=c(-15, cbreaks, 15)
colscale1=list(cbreaks=cbreaks, col_breaks=col_breaks, legend=c(-4,-2,0,2,4), scale=greenred(64))
cbreaks=seq(0,4,length.out=63)
col_breaks=c(-15, cbreaks, 15)
colscale2=list(cbreaks=cbreaks, col_breaks=col_breaks, legend=c(0,2,4), scale=colorpanel(64, low="#000000", high="#FF0000"))
colscale=list(colscale1, colscale1, colscale2)

splitdata=split_by_class(sdata, classdata)

playout=rbind(c(2,4,6), c(1,3,5))
rownames(playout)=c("1500", "500")
colnames(playout)=c("IP", "Fat", "Cap")
#playout=t(playout)
make_grouped_plot("cell_types.eps", splitdata, grps, playout, colscale, byrow=2, labrow=3)
allgrp=rbind(grps[1,], rep("All", nrow(grps)), grps[3,])
make_grouped_plot("everything.eps", splitdata, allgrp, playout, colscale, byrow=2, labrow=3)
#make_break_plots(sdata, 1, "cyto_data")

