setwd("C:/Users/tamhok/Documents/MIT/Anderson/Fibrosis/Smallmolec")
library("gplots")
library("RCurl")
source("spec_pheatmap.R")
library("grid")
toHTML = function (x) {
	x=gsub("^\\s+|\\s+$", "", x)
	x=gsub("\\s", "%20", x)
	x=gsub("\\(", "%28",x)
	x=gsub("\\)", "%29",x)
	return(x)
}

getSMILES=function(name) {
	name=toHTML(name)
	url=paste("http://www.emolecules.com/lookup?q=", name, sep="")
	output=getURL(url);
	if(output=="__END__\n") {
		return(NA)
	}
	rv=strsplit(output, "\t")[[1]][1:2]
	return(rv)
}



proc_bb_struct=function(x) {
	if(x[1] == "N/A") {
		return (x[2:3])
	} else if(x[3] == "N/A") {
		return (x[1:2])
	} else {
		return (x[c(1,3)])
	}
}

#gets ids of duplicated elements
get_duplicates=function(x, do_t=FALSE) {
  if(do_t) x=t(x)
  return(which(duplicated(x) | duplicated(x, fromLast=TRUE)))
}

#Combines duplicates in COLUMNS based on indexing columns index_cols
#Does comb_fn on each set of duplicates on a per-column basis using
#comdb_cols
#Use do_t argument to take transpose
combine_duplicates=function(data, index_cols, comb_cols, comb_fn, do_t=FALSE) {
  if(do_t) data=t(data)
  
  #Form a single factor
  index=data[,index_cols]
  new_ind=do.call(paste, c(index, sep=":"))
  
  #Combine everything
  cdat=Reduce(cbind, lapply(comb_cols, function(x) {tapply(data[,x], new_ind, comb_fn)}))
  
  #make new names
  new_ind=Reduce(rbind, strsplit(rownames(cdat), ":"))

  #stick together
  rdata=cbind(new_ind, cdat)
  colnames(rdata)=c(colnames(index), colnames(data)[comb_cols])
  if(do_t) rdata=t(rdata) 
  return(rdata)
}

make_heatmap=function(data, filename, raw_data) {
  codes=read.csv("ht_codes.csv", stringsAsFactors=FALSE, header=TRUE, strip.white=TRUE)
  rownames(codes)=codes[,1]
  b1=as.factor(data[1,])
  b2=as.factor(data[2,])
  levels(b1)=codes[levels(b1),2]
  levels(b2)=codes[levels(b2),2]
  print(length(levels(b1)))
  print(length(levels(b2)))
  molec_mat=matrix(data=-50, nrow=length(levels(b1)), ncol=length(levels(b2)))
  rownames(molec_mat)=levels(b1)
  colnames(molec_mat)=levels(b2)
  for(i in 1:nrow(raw_data)) {
    molec_mat[b1[i], b2[i]]=raw_data[i,5]
  }
  print(max(molec_mat))
  pheatmap(t(molec_mat), color=c("#000000", colorpanel(74, "#00FF00", "#FFFF00", "#FF0000")), 
           breaks=c(-0.0001, seq(0.001,0.8,length.out=40), seq(0.81,1,length.out=20), seq(1.1, 4, length.out=14)),
           legend_breaks=c(0,0.25,0.5,0.75,1, 2, 3, 4), 
           legend_labels=c("0", "0.25", "0.5", "0.75", "1 (VLVG)", "2", "3", "4"),
           cellwidth = 15, cellheight = 12, fontsize = 8, filename = filename)
#   pdf(filename, width=6, height=12)
#   heatmap.2(t(molec_mat), trace="none", col=c("#0000FF", "#FFFFFF", greenred(62)),
#             breaks=c(-0.001,-0.0001,seq(0.00001,1,length.out=63)), margins=c(10,0), 
#             cexRow=0.5, cexCol=0.4, lhei=c(1,7), lwid=c(1,5,2), 
#             lmat=rbind(c(0,3,4), c(2,1,0)), density.info="none", 
#             distfun=dist, xlab="Base")
#   dev.off()  
}

create_map_from_file=function(file, outfile) {
  raw_data=read.csv(file, header=TRUE, stringsAsFactors=FALSE, strip.white=TRUE)
  new_bb=apply(raw_data[,1:3], 1, proc_bb_struct)
  make_heatmap(new_bb, paste(outfile, ".eps", sep=""), raw_data)  
}

create_map_from_file("ht_library.csv","ht_clust")
create_map_from_file("ht_library_1.csv","ht_clust1")
create_map_from_file("ht_library_2.csv","ht_clust2")
create_map_from_file("azido.csv","azido")
create_map_from_file("iodo.csv","iodo")
# newer_bb=new_bb[,-get_duplicates(new_bb, do_t=TRUE)]
# make_heatmap(newer_bb, "ht_clust2.pdf")
# dup_bb=new_bb[, get_duplicates(new_bb, do_t=TRUE)]
# make_heatmap(dup_bb, "ht_clust_dup.pdf")
# cbreaks=seq(0.000001,4,length.out=63)
# col_breaks=c(0, cbreaks, 10)
# pheatmap(molec_mat, col=greenred(64), treeheight_row=100, treeheight_col=100,
	# legend_breaks	c(0,1,2,3,4), legend_labels=c(0,1,2,3,4), filename="ht_clust.pdf",
		# family="Sans", fontsize=18, cellwidth=36, cellheight=20, 
		# breaks=col_breaks, rbreaks=cbreaks, scale="none", legend=TRUE)