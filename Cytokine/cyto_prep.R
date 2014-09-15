source("C:/Users/tamhok/Documents/MIT/Anderson/Fibrosis/R/nanostring_prep.R")
setwd("C:/Users/tamhok/Documents/MIT/Anderson/Fibrosis/Cytokine")

nano_data=load_data("raw_cytomur", 3);
proc_data=average_data(nano_data$raw, nano_data$class, nano_data$limit); 
#Determine which scalings to use
scaling=vector(length=nrow(proc_data$pdata))
labs=rownames(proc_data$pdata);
ids=sub("_*", "", labs)
proc_data$pdata=proc_data$pdata[c(3,1,2),]
scaling=rep(1, nrow(proc_data$pdata))
sdata=scale_data(proc_data$pdata, scaling, mean(nano_data$limit))
#sdata=remove_zero_data(sdata, proc_data$zdata)
sums=sdata[2,]+sdata[1,]
sdata=sdata[,order(sums, decreasing=TRUE)]
make_break_plots(sdata, 4, "cyto_data")