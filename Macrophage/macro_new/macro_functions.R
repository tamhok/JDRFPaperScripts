setwd("C:/Users/tamhok/Documents/MIT/Anderson/Fibrosis/Macrophage/macro_new");
library('NanoStringNorm')
library('ggplot2')
library('plotrix')
library('gplots')
library('reshape')
source("C:/Users/tamhok/Documents/MIT/Anderson/Fibrosis/R/nanostring_prep.R")
source("macro_pheatmap.R")

#Loads the data from nSolver Export and normalizes it
#writes files
#identifiers gives the number of attributes of each set (size, days, etc.)
#Returns list of raw data data.frame and classes
load_data=function(filename, identifiers) {
	data=read.csv(paste(filename, ".csv", sep=""), header=TRUE, stringsAsFactors=FALSE);
	data[data$Name %in% 
		c("Cltc", "Gapdh", "Gusb", "Hprt1", "Pgk1", "Tubb5", "CLTC", "HPRT1","Tub5"),
		 'Code.Class'] = 'Housekeeping';
	data=data.frame(data)
	classes=data[1:identifiers,-c(1:3)]
	nano_data=cbind(data[-c(1:identifiers),1:3], data.matrix(data[-c(1:identifiers),-c(1:3)]))
	raw_data=NanoStringNorm(nano_data, 
				anno=NA, 
				CodeCount = 'geo.mean',
				SampleContent='housekeeping.sum', 
				Background='mean',
				#return.matrix.of.endogenous.probes=TRUE,
				);
	raw=data.frame(raw_data$normalized.data);
	negs=raw[raw$Code.Class=="Negative",4:ncol(raw)]
	negs=apply(negs, 2, mean)
	raw=raw[raw$Code.Class=="Endogenous",c(2,4:ncol(raw))]
	rownames(raw)=raw[,1]
	raw=as.matrix(raw[,-1])
	raw=t(t(raw)-negs)
	raw[raw<0]=0
	raw2=rbind(classes, raw)
	write.csv(raw2, paste(filename, "_proc.csv", sep=""))
	write.csv(classes, paste(filename, "_classes.csv", sep=""))
	return(list("raw"=raw, "class"=classes, "limit"=negs));
}

#Create a vector of names for each "type"
create_class_names=function(classes)
	return(apply(classes, 2, function(x) paste(x, collapse="_")))
end

#Averages the data into groups
#Returns a matrix with columns of genes and rows of avgd samples
#and matrix with same but of percent zeros.
average_data=function(raw_data, classes, negs) {
	class_names=as.factor(create_class_names(classes));
	lvls=levels(class_names);
	pdata=matrix(ncol=nrow(raw_data), nrow=length(lvls))
	zdata=matrix(ncol=nrow(raw_data), nrow=length(lvls))
	neg_pooled=vector(length=length(levels))
	for (i in 1:length(lvls)) {
		indices=which(class_names==lvls[i])
		sel_data=raw_data[,indices]
		neg_pooled[i]=mean(negs[indices])
		pdata[i,]=apply(sel_data,1,mean);
		zdata[i,]=apply(sel_data,1,function(x) sum(x < 0.001)/length(x));
	}
	rownames(pdata)=lvls;
	colnames(pdata)=rownames(raw_data)
	rownames(zdata)=lvls;
	colnames(zdata)=rownames(raw_data)
	return(list("pdata"=pdata,"zdata"=zdata,"class"=class_names, "negs"=neg_pooled))
}

#Does the scaling. Requires pdata and a vector of which things
#to use as scales. Also puts it in log scale.
scale_data=function(pdata, scaling, scaling_mat=pdata) {
	scale_fun=function(x) {
	  ids=which(scaling_mat[x[ncol(pdata)+1],]<0.001)
    if(length(ids) > 0) {
      under=x[ids]
      x=x[1:ncol(pdata)]/ scaling_mat[x[ncol(pdata)+1],]
      under[under < 1]=1
      x[ids]=under
      return(x[1:ncol(pdata)])
    } else {
      return(x[1:ncol(pdata)]/ scaling_mat[x[ncol(pdata)+1],])
    }
	}
  sdata=apply(cbind(pdata, scaling), 1, scale_fun )
	sdata=t(sdata)
	sdata[sdata<0.001]=0.001
	return(log(sdata));
}

#Removes rows with all > thresh_indiv (fraction) zeros
remove_zero_data=function(sdata, zdata, thresh_indiv=0.5, thresh_overall=0.6) {
	elim_data=apply(zdata, 2, function (x) sum(x > thresh_indiv)) < thresh_overall * nrow(zdata) 
	rdata=sdata[,elim_data]
	return(rdata);
}

obtain_group_data=function(rdata, grps, ranker) {
	if(is.vector(rdata)) {
		rdata=t(as.matrix(rdata))
	}
	grp.in=which(grps[,1] %in% colnames(rdata))
	grps=grps[grp.in,];
	udata=rbind(rdata, ranker)
	udata=t(udata[,grps[,1]])
	grps=cbind(grps, udata)
	return(grps)
}

make_flats=function(rdata, grps) {
	grp.in=which(grps[,1] %in% colnames(rdata))
	grps=grps[grp.in,];
	ranker=apply(rdata,2, sum);
	udata=rbind(rdata, ranker)
	udata=t(udata[,grps[,1]])
	grps=cbind(grps, udata);
	
	ngrps=levels(as.factor(grps[,2]));
	ngens=length(grps[,1]);
	
	grid.newpage();
	break_width=unit(0.02, "npc");
	panel_width=unit(0.3, "npc");
	widths=unit.c(panel_width, break_width, panel_width, break_width, panel_width, break_width)
	heights=unit(1, "npc")
	pushViewport(viewport(layout = grid.layout(nrow = 1, ncol = 6, widths=widths, heights=heights), gp = gpar(fontsize=18)))
	
	for (i in 1:length(ngrps)) {		
		sg=which(grps[,2]==ngrps[i]);
		this.data=grps[sg,4:(ncol(grps)-1)];
		rownames(this.data)=grps[sg,3]
		d.order=order(grps[sg,ncol(grps)]);		
		this.data=this.data[d.order,];
		
		pushViewport(vplayout(1, (i-1)*2+1, just=c("left", "top")))
		
		pheatmap(data.matrix(this.data), col=greenred(64), cluster_rows=FALSE, 
		cluster_cols=FALSE, family="Sans", fontsize=18, cellwidth=36,
		cellheight=20, show_colnames=TRUE, main=ngrps[i],
		breaks=col_breaks, scale="none", legend=FALSE)
		upViewport()
	}

	pushViewport(vplayout(1,6, just="bottom"))
	legend = c(-4,-2,0 ,2,4)
	names(legend)=legend
	draw_legend(color=greenred(64), ht=0.95, cbreaks, legend, fontsize = 18, family="Sans")
	upViewport()
}

do_stats=function(nano_data, grps, filt_word = "") {
  #Do statistics
  raw_mat=nano_data$raw
  cnames=create_class_names(nano_data$class)
  colnames(raw_mat)=cnames
  raw_data=melt(raw_mat)
  raw_grps=t(rbind(apply(raw_data, 1, function(x) unlist(strsplit(x[2], "_")))))
  
  grp_factors=as.factor(grps[2,])
  grp_lvls=levels(grp_factors)
  if(grp_lvls[1]=="") {
    grp_lvls=grp_lvls[-1]
  }
  
  gene_class=vector(mode="character", length=length(raw_data[,1]))
  for(i in 1:length(grp_lvls)) {
    gene_list=grps[1,which(grps[2,]==grp_lvls[i])]
    index=which(raw_data[,1] %in% gene_list)
    gene_class[index]=grp_lvls[i]
  }
  vals=log(raw_data[,3]+0.00001)
  fdata=data.frame(vals=vals, genes=raw_data[,1], days=raw_grps[,2],
                   size=raw_grps[,1], comp=raw_grps[,3], class=gene_class)
  
  f_filt=fdata[-which(fdata$class==filt_word),]
  
  model=aov(vals ~ size*days*comp*class + genes, data=f_filt)
  class_order=c("M[Classic]", "M[Reg]", "M[Wound]")
  p_vals=vector("list", length(comp_order)*length(size_order))
  ctr=1
  for(i in 1:length(comp_order)) {
    for (j in 1:length(size_order)) {
      temp_pvals=vector("list", length(class_order))
      for (k in 1:length(class_order)) {
        cur_set=f_filt[f_filt$comp==comp_order[i] & f_filt$size==size_order[j] & f_filt$class==class_order[k],]
        C_set=f_filt[f_filt$comp==comp_order[i] & f_filt$days=="C" & f_filt$class==class_order[k],]
        test_set=rbind(C_set, cur_set)
        #print(levels(as.factor(test_set$class)))
        model=aov(vals ~ days + genes, data=test_set)
        print(paste(comp_order[i], size_order[j], class_order[k]))
        tukey=TukeyHSD(model)
        stuk=tukey$days[unlist(lapply(rownames(tukey$days), function(x) substr(x, 1, 1)=="C")),]
        print(tukey$days)
        pval=min(stuk[,"p adj"])*18
        mod_name=class_order[k]
        if(pval< 0.000005) {
          mod_name=paste(mod_name, "***")
        } else if(pval < 0.0005) {
          mod_name=paste(mod_name, "**")
        } else if(pval < 0.05) {
          mod_name=paste(mod_name, "*")
        }
        temp_pvals[[k]]=list(name=class_order[k], pval=pval, mod_name=mod_name)
      }
      p_vals[[ctr]]=temp_pvals
      ctr=ctr+1
    }
  }
  return(p_vals)
}

do_stats=function(nano_data, grps, filt_word = "") {
  #Do statistics
  raw_mat=nano_data$raw
  cnames=create_class_names(nano_data$class)
  colnames(raw_mat)=cnames
  raw_data=melt(raw_mat)
  raw_grps=t(rbind(apply(raw_data, 1, function(x) unlist(strsplit(x[2], "_")))))
  
  grp_factors=as.factor(grps[2,])
  grp_lvls=levels(grp_factors)
  if(grp_lvls[1]=="") {
    grp_lvls=grp_lvls[-1]
  }
  
  gene_class=vector(mode="character", length=length(raw_data[,1]))
  for(i in 1:length(grp_lvls)) {
    gene_list=grps[1,which(grps[2,]==grp_lvls[i])]
    index=which(raw_data[,1] %in% gene_list)
    gene_class[index]=grp_lvls[i]
  }
  vals=log(raw_data[,3]+0.001)
  fdata=data.frame(vals=vals, genes=raw_data[,1], days=raw_grps[,2],
                   size=raw_grps[,1], comp=raw_grps[,3], class=gene_class)
  
  f_filt=fdata[-which(fdata$class==filt_word),]
  
  model=aov(vals ~ size+days+comp+class + genes, data=f_filt)
  class_order=c("M[Classic]", "M[Reg]", "M[Wound]")
  p_vals=vector("list", length(comp_order)*length(size_order))
  ctr=1
  for(i in 1:length(comp_order)) {
    for (j in 1:length(size_order)) {
      temp_pvals=vector("list", length(class_order))
      for (k in 1:length(class_order)) {
        cur_set=f_filt[f_filt$comp==comp_order[i] & f_filt$size==size_order[j] & f_filt$class==class_order[k],]
        C_set=f_filt[f_filt$comp==comp_order[i] & f_filt$days=="C" & f_filt$class==class_order[k],]
        test_set=rbind(C_set, cur_set)
        #print(levels(as.factor(test_set$class)))
        model=aov(vals ~ days + genes, data=test_set)
        print(paste(comp_order[i], size_order[j], class_order[k]))
        tukey=TukeyHSD(model)
        stuk=tukey$days[unlist(lapply(rownames(tukey$days), function(x) substr(x, 1, 1)=="C")),]
        print(tukey$days * 18)
        pval=min(stuk[,"p adj"])*18
        mod_name=class_order[k]
        if(pval< 0.000005) {
          mod_name=paste(mod_name, "***")
        } else if(pval < 0.0005) {
          mod_name=paste(mod_name, "**")
        } else if(pval < 0.05) {
          mod_name=paste(mod_name, "*")
        }
        temp_pvals[[k]]=list(name=class_order[k], pval=pval, mod_name=mod_name)
      }
      p_vals[[ctr]]=temp_pvals
      ctr=ctr+1
    }
  }
  return(p_vals)
}

do_stat_alt = function(nano_data, grps, filt_word = "") {
  #Do statistics
  raw_mat=nano_data$raw
  cnames=create_class_names(nano_data$class)
  colnames(raw_mat)=cnames
  raw_data=melt(raw_mat)
  raw_grps=t(rbind(apply(raw_data, 1, function(x) unlist(strsplit(x[2], "_")))))
  
  grp_factors=as.factor(grps[2,])
  grp_lvls=levels(grp_factors)
  if(grp_lvls[1]=="") {
    grp_lvls=grp_lvls[-1]
  }
  
  gene_class=vector(mode="character", length=length(raw_data[,1]))
  for(i in 1:length(grp_lvls)) {
    gene_list=grps[1,which(grps[2,]==grp_lvls[i])]
    index=which(raw_data[,1] %in% gene_list)
    gene_class[index]=grp_lvls[i]
  }
  vals=log(raw_data[,3]+0.001)
  fdata=data.frame(vals=vals, genes=raw_data[,1], days=raw_grps[,2],
                   size=raw_grps[,1], comp=raw_grps[,3], class=gene_class)
  
  f_filt=fdata[-which(fdata$class==filt_word),]
  
  model=aov(vals ~ size*days*comp*class + genes, data=f_filt)
  class_order=c("M[Classic]", "M[Reg]", "M[Wound]")
  day_order = c(1,4,7)
  p_vals=vector("list", length(comp_order)*length(size_order))
  ctr=1
  for(i in 1:length(comp_order)) {
    for (j in 1:length(day_order)) {
        cur_set=f_filt[f_filt$comp==comp_order[i] & f_filt$days == day_order[j] &
                         (f_filt$size == 500 | f_filt $ size == 1500),]
        C_set=f_filt[f_filt$comp==comp_order[i] & f_filt$size=="C",]
        cur_set=rbind(C_set, cur_set)
        #print(levels(as.factor(test_set$class)))
        model=aov(vals ~ size + genes, data=cur_set)
        print(paste(comp_order[i], day_order[j]))
        tukey=TukeyHSD(model)
        print(tukey$size)
    }
  }
}
