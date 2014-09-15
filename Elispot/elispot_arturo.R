setwd("C:/Users/tamhok/Documents/MIT/Anderson/Fibrosis/Elispot");
library('NanoStringNorm')
library('ggplot2')
library('plotrix')
library('gplots')
library('reshape')
source("C:/Users/tamhok/Documents/MIT/Anderson/Fibrosis/R/nanostring_prep.R")

#Loads the data from nSolver Export and normalizes it
#writes files
#identifiers gives the number of attributes of each set (size, days, etc.)
#Returns list of raw data data.frame and classes
load_data=function(filename, identifiers) {
  data=read.csv(paste(filename, ".csv", sep=""), header=TRUE, stringsAsFactors=FALSE);
  data=data.frame(data)
  classes=data[1:identifiers,-1]
  l_data=data[-c(1:identifiers),-1]
  rownames(l_data)=data[-c(1:identifiers),1]
  return(list(data=data.matrix(l_data), classes=classes))
}

#Remove background data & scale to positive control
scale_data=function(data, pos_id, neg_id) {
  scale_fun=function(x) {
    pos_mean=mean(x[pos_id])
    neg_mean=mean(x[neg_id])
    tehdata=x[c(-pos_id, -neg_id)]
    return(unlist(lapply(tehdata, function(x) max((x-neg_mean)/(pos_mean-neg_mean),0))))
  }
  return(apply(data, 2, scale_fun))
}

data=load_data("elispot_arturo", 1)
data$data=log2(data$data)
pos_id=which(rownames(data$data) %in% c("Pos1", "Pos2", "Pos3"))
neg_id=which(rownames(data$data) %in% c("Negative Ctrl"))
scaled_data=scale_data(data$data, pos_id, neg_id)
p_data=average_data(scaled_data, data$classes, 0)

index=0;
classes=p_data$class_lst
#build class data
classdata=cbind(t(classes), rep(1, ncol(classes)))

#cbreaks=seq(0,1,length.out=63)
cbreaks=c(seq(0,0.05, length.out=32), seq(0.06,1, length.out=31))
col_breaks=c(-15, cbreaks, 15)
colscale1=list(cbreaks=cbreaks, col_breaks=col_breaks, legend=c(0,0.25,0.5,0.75,1), scale=colorpanel(64, low="#000000", high="#FF0000"))
colscale=list(colscale1, colscale1)

splitdata=split_by_class(p_data$pdata, classdata)

playout=as.matrix(c(1))
rownames(playout)=c("")
colnames(playout)=c("")
#playout=t(playout)
#new_names=do_stats(nano_data, grps)
#make_grouped_plot("cell_types.eps", splitdata, grps, playout, colscale, byrow=2, labrow=3, opt_names=new_names)
allgrp=rbind(colnames(p_data$pdata), rep(" ", ncol(p_data$pdata)), colnames(p_data$pdata))
make_grouped_plot("everything_arturo.eps", splitdata, allgrp, playout, colscale, byrow=2, labrow=3, cluster=TRUE)

# ratio_data=splitdata[[2]]$data - splitdata[[1]]$data
# rat_data=list(list(data=ratio_data, name="1"))
# playout=as.matrix(1)
# rownames(playout)=c("1500-500")
# colnames(playout)=c("")
# playout=t(playout)
# cbreaks=seq(-0.7,0.7,length.out=63)
# col_breaks=c(-15, cbreaks, 15)
# colscale1=list(cbreaks=cbreaks, col_breaks=col_breaks, legend=c(-0.7,-0.35,0,0.35,0.7), scale=greenred(64))
# colscale=list(colscale1)
# make_grouped_plot("ratio_arturo.eps", rat_data, allgrp, playout, colscale, byrow=2, labrow=3)
