source("C:/Users/tamhok/Documents/MIT/Anderson/Fibrosis/R/nanostring_prep.R")
setwd("C:/Users/tamhok/Documents/MIT/Anderson/Fibrosis/jd_nano")

nano_data=load_data("jd_nano", 3);
proc_data=average_data(nano_data$raw, nano_data$class, nano_data$limit); 
#Determine which scalings to use
scaling=vector(length=nrow(proc_data$pdata))
classes=proc_data$class_lst;
scaled_data=proc_data$pdata[c(which(classes[4,]=="Mock_WT_Fat")),]
scaled_data[scaled_data < 1]=1;

scaled_data=rbind(scaled_data, rep(1, ncol(proc_data$pdata)));
scaling = rep(1, dim(classes)[2])
sdata=scale_data(proc_data$pdata, scaling, mean(nano_data$limit), scaling_mat=scaled_data)
sdata=remove_zero_data(sdata, proc_data$zdata)
sums=apply(sdata, 2, sum)
sdata=sdata[order(classes[3,], classes[1,], classes[2,]),order(sums, decreasing=TRUE)]

#Create scaled data for each type of mouse
mocks = classes[,which(classes[1,]=="Mock")]
scaled_data2 = proc_data$pdata[which(classes[1,]=="Mock"),]
scaled_data[scaled_data < 1] = 1

assign_elt = function(x) {
  elt = which(mocks[2,] == x[2])
  if(any(elt)) return(elt)
  return(1)
}

scaling2 = apply(classes, 2, assign_elt)
sdata2=scale_data(proc_data$pdata, scaling2, mean(nano_data$limit), scaling_mat=scaled_data2)
sdata2=remove_zero_data(sdata2, proc_data$zdata)
classdata2 = cbind(classes[4,], classes[4,], scaling2)
# 
# #load groups
# grps=load_groups("fibrosis_groups")
# 
# comp_order=c("IP", "Fat", "Cap")
# size_order=c("500", "1500")
# clsnames=classes[4,]
# index=0;
# 
# cdata=rbind(nano_data$class, nano_data$raw)
# 
# take_anova=function(data) {
#   dats=log(as.numeric(data))
#   newdf=data.frame(days,dats)
#   return(anova(aov(dats~days, newdf))$'Pr(>F)'[[1]])
# }
# 
# for(i in 1:length(comp_order)) {
#   for (j in 1:length(size_order)) {
#     temp_data=cdata[c(-1,-3),which(cdata[1,]==size_order[j] & cdata[3,]==comp_order[i])]
#     days=as.vector(t(temp_data[1,]))
#     
#   }
# }

classdata=cbind(classes[4,], classes[4,], rep(1, ncol(classes)))
  
# #build class data
# classdata=matrix(nrow=0,ncol=3);
# for(i in 1:length(comp_order)) {
# 	cur_set=classes[,which(classes[3,]==comp_order[i])]
# 	for (j in 1:length(size_order)) {
# 		sub_set=cur_set[, which(cur_set[1,]==size_order[j])]
# 		sub_set=sub_set[c(2,4), order(sub_set[2,])]
# 		mock_set=cur_set[c(1,4), which(cur_set[1,]=="Mock"), drop=FALSE]
# 		if(length(mock_set) > 0) {
# 			sub_set=cbind(mock_set, sub_set);
# 		}
# 		index=index+1;
# 		sub_set=t(rbind(sub_set, rep(index, ncol(sub_set))));
# 		classdata=rbind(classdata, sub_set)
# 	}
# }
cbreaks=seq(-4,4,length.out=63)
col_breaks=c(-15, cbreaks, 15)
colscale1=list(cbreaks=cbreaks, col_breaks=col_breaks, legend=c(-4,-2,0,2,4), scale=greenred(64))
cbreaks=seq(0,8,length.out=63)
col_breaks=c(-15, cbreaks, 15)
colscale2=list(cbreaks=cbreaks, col_breaks=col_breaks, legend=c(0,2,4,6,8), scale=colorpanel(64, low="#000000", high="#FF0000"))
#colscale=list(colscale1, colscale1, "autosym")
colscale=list(colscale1, colscale1)
splitdata=split_by_class(sdata, classdata)

playout=as.matrix(1);
rownames(playout)="";
colnames(playout)="";

# playout=rbind(c(2,4), c(1,3))
# rownames(playout)=c("1500", "500")
# colnames(playout)=c("IP", "Fat")
#playout=t(playout)
#make_grouped_plot("markers.eps", splitdata, grps, playout, colscale, byrow=2)
#make_grouped_plot("activation.eps", splitdata, grps, playout, colscale, byrow=3)
allgrp=rbind(colnames(sdata), rep("All", length(colnames(sdata))), colnames(sdata))
make_grouped_plot("most.eps", splitdata, allgrp, playout, colscale, byrow=2, labrow=3)
#make_break_plots(sdata, 1, "cyto_data")

splitdata2 = split_by_class(sdata2, classdata2)
playout2 = as.matrix(t(1:ncol(mocks)))
rownames(playout2)= c("")
colnames(playout2) = mocks[2,]
colscale2 = rep(list(colscale1), ncol(mocks))
allgrp=rbind(colnames(sdata2), rep("All", length(colnames(sdata2))), colnames(sdata2))
make_grouped_plot("most_grpd.eps", splitdata2, allgrp, playout2, colscale2, byrow=2, labrow=3)
