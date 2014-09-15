library('NanoStringNorm')
library('ggplot2')
library('plotrix')
library('gplots')
library('grid')
source("C:/Users/tamhok/Documents/MIT/Anderson/Fibrosis/pheatmap.R")

#Loads the data from nSolver Export and normalizes it
#writes files
#identifiers gives the number of attributes of each set (size, days, etc.)
#Returns list of raw data data.frame and classes
load_data=function(filename, identifiers) {
  data=read.csv(paste(filename, ".csv", sep=""), header=TRUE, stringsAsFactors=FALSE);
  data[data$Name %in% 
         c("Bact", "Cltc", "Gapdh", "Gusb", "Hprt1", "Pgk1", "Tubb5", "CLTC", "HPRT1","Tub5", "Tubb5"),
       'Code.Class'] = 'Housekeeping';
#  data[data$Name %in% 
#         c("Gapdh", "Gusb", "Hprt1", "Pgk1", "Tubb5", "CLTC", "HPRT1","Tub5", "Tubb5"),
#       'Code.Class'] = 'Endogenous';
  data=data.frame(data)
  classes=data[1:identifiers,-c(1:3)]
  nano_data=cbind(data[-c(1:identifiers),1:3], data.matrix(data[-c(1:identifiers),-c(1:3)]))
  raw_data=NanoStringNorm(nano_data, 
                          anno=NA, 
                          CodeCount = 'geo.mean',
                          SampleContent='housekeeping.geo.mean', 
                          Background='none',
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
create_class_names=function(classes) {
	classes=as.matrix(classes)
  if(dim(classes)[1] > 1 && dim(classes)[2] > 1) {
    return(apply(classes, 2, function(x) paste(x, collapse="_")))
	} else {
    return(classes)
	}
}

#Averages the data into groups
#Returns a matrix with columns of genes and rows of avgd samples
#and matrix with same but of percent zeros.
average_data=function(raw_data, classes, negs) {
  rm_dup_classes=classes[!duplicated(t(classes))]
  rm_dup_classes=as.matrix(rm_dup_classes)
  if(dim(rm_dup_classes)[2]==1) {
    rm_dup_classes=t(rm_dup_classes)
  }    
  rm_dup_classes=as.matrix(rbind(rm_dup_classes, create_class_names(rm_dup_classes)))	
  lvls=rm_dup_classes[dim(rm_dup_classes)[1],]                          
  class_names=create_class_names(classes);
  pdata=matrix(ncol=nrow(raw_data), nrow=length(lvls))
	zdata=matrix(ncol=nrow(raw_data), nrow=length(lvls))
	neg_pooled=vector(length=length(lvls))
	for (i in 1:length(lvls)) {
		indices=which(class_names==lvls[i])
		sel_data=as.matrix(raw_data[,indices])
		neg_pooled[i]=mean(negs[indices])
		pdata[i,]=apply(sel_data,1,mean);
		zdata[i,]=apply(sel_data,1,function(x) sum(x < 0.001)/length(x));
	}
	rownames(pdata)=lvls;
	colnames(pdata)=rownames(raw_data)
	rownames(zdata)=lvls;
	colnames(zdata)=rownames(raw_data)
	return(list("pdata"=pdata,"zdata"=zdata,"class"=lvls, "class_lst"=rm_dup_classes, "negs"=neg_pooled))
}

#Does the scaling. Requires pdata and a vector of which things
#to use as scales. Also puts it in log scale.
scale_data=function(pdata, scaling, limit, scaling_mat=pdata) {
	scaling_mat[scaling_mat==0]==limit;
	scale_fun=function(x) {
		ids=which(scaling_mat[x[ncol(pdata)+1],] > 0.001) 
		x[ids]=x[ids] / scaling_mat[x[ncol(pdata)+1],ids]
		x[-ids]=(x[-ids] + limit) / limit
		j=x[-ids]
		j[j==1]=NaN
		x[-ids]=j
		return(x[1:ncol(pdata)])
	}
	sdata=apply(cbind(pdata, scaling), 1, scale_fun)
	sdata=t(sdata)
	sdata[sdata<0.01]=0.001
	return(log2(sdata));
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
	class(udata)="numeric"
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

#load group data
load_groups=function(filename) {
	groups=read.csv(paste(filename, ".csv", sep=""), header=TRUE, stringsAsFactors=FALSE);
	return(t(groups));
}



split_by_class=function(pdata, classdata) {
	#identify all possible groups
	lvls=levels(as.factor(classdata[,3]))
	
	#remove blank elements
	if(lvls[1]=="") {
		lvls=lvls[-1]
	}
	
	#iterate through groups and pick out relevant ones from pdata
	rv=vector("list", length(lvls))
	for(i in 1:length(lvls)) {
		if(lvls[i]=="") {
			next;
		}
		members=which(classdata[,3]==lvls[i])
		ids=classdata[members,2]
		data=as.matrix(pdata[ids,,drop=FALSE])
		rownames(data)=classdata[members,1]
		rv[[i]]=list(data=data, name=lvls[i])
	}
	return(rv)
}

#form subgroups
make_subgroups=function(pdata, groups, byrow=2, labrow=4) {
	#Check to ensure groups in pdata
  groups=groups[,which(groups[1,] %in% colnames(pdata))]
  name_order=1:dim(pdata)[2]
  names(name_order)=colnames(pdata)
  #identify all possible groups
	lvls=levels(as.factor(groups[byrow,]))
	
	#remove blank elements
	if(lvls[1]=="") {
		lvls=lvls[-1]
	}
	
	#iterate through groups and pick out relevant ones from pdata
	rv=vector("list", length(lvls))
	for(i in 1:length(lvls)) {
		if(lvls[i]=="") {
			next;
		}
		members=which(groups[byrow,]==lvls[i])
		ids=groups[1,members]
    temp_order=name_order[ids]
		data=as.matrix(pdata[,ids, drop=FALSE])
		colnames(data)=groups[labrow, members]
    data=data[,order(temp_order),drop=FALSE]
		rv[[i]]=list(data=data, name=lvls[i])
	}

	return(rv)
}

make_grouped_plot=function(filename, grp_data, group_data, plot_layout, colscales, byrow=2, labrow=4, opt_names=NULL, cluster=FALSE, show_dend_col=FALSE, show_dend_row=FALSE, flip_data=FALSE) {
	n_subp=length(grp_data);
		
	datalist=vector("list", n_subp);
	panel_sizes=matrix(ncol=2, nrow=n_subp)
	max_val_per_col=vector(length=ncol(plot_layout))
	
	cbreaks=seq(-3,3,length.out=63)
	col_breaks=c(-15, cbreaks, 15)
	scale=greenred(64)
	
	for(i in 1:ncol(plot_layout)) {
		colmaxval=0;
		for(j in 1:nrow(plot_layout)) {
			num=plot_layout[j,i]
			if(num > 0) {				
				new_data=make_subgroups(grp_data[[num]]$data, group_data, byrow, labrow)
				if(!is.null(opt_names)) {
				  for(k in 1:length(new_data)) {
            new_data[[k]]$name=opt_names[[num]][[k]]$mod_name
				  }
				}
        tmaxval=max(unlist(lapply(new_data, function(x) {max(x$data)})))
				if(colmaxval < tmaxval) {
					colmaxval=tmaxval;
				}
				sizes=make_single_plot(new_data, 0, showcolns=(j==TRUE),
                               col_breaks, cbreaks, scale_cols, cluster, flip_data=flip_data,
				                       show_dend_col=show_dend_col, show_dend_row=show_dend_row)
				panel_sizes[num,]=c(sizes$min_width, sizes$min_height)
				
        datalist[[num]]=list(data=new_data, size=sizes);

			}
		}
		max_val_per_col[i]=colmaxval;
	}
	
	max_val_per_col=ceiling(max_val_per_col*2)/2
	
	p_widths=vector();
	p_heights=vector();
	
	for (i in 1:ncol(plot_layout)) {
		p_widths=c(p_widths, max(panel_sizes[plot_layout[,i],1]))
	}
	
	for (i in 1:nrow(plot_layout)) {
		p_heights=c(p_heights, max(panel_sizes[plot_layout[i,],2]))
	}

	min_width=sum(p_widths)+ncol(plot_layout)*3
	min_height=sum(p_heights)+1*nrow(plot_layout)
	postscript(filename, horizontal = FALSE, onefile = FALSE, paper = "special", width=min_width, height=min_height)
	grid.newpage();
	
	title_ht=unit(48, "bigpts")
	title_w=unit(48, "bigpts")
	
	break_ht=unit(1, "inches");
	
	hts=unit.c(title_ht,unit(p_heights[1], "inches"))
	if(nrow(plot_layout) > 1) {
		for(i in 2:nrow(plot_layout)) {
			p_ht=unit(p_heights[i], "inches")
			hts=unit.c(hts, break_ht, p_ht)
		}
	}
	
	break_w=unit(1.5, "inches");
	p_w=unit(p_widths[1], "inches")
	wds=unit.c(title_w, p_w)
	wds=unit.c(wds, break_w)
	if(ncol(plot_layout) > 1) {
		for(i in 2:ncol(plot_layout)) {
			p_w=unit(p_widths[i], "inches")
			wds=unit.c(wds, p_w, break_w)
		}
	}
	pushViewport(viewport(layout = grid.layout(nrow = length(hts), ncol = length(wds), widths=wds, heights=hts), gp = gpar(fontsize=18)))

	#plot row titles
	for (i in 1:nrow(plot_layout)) {
		pushViewport(vplayout(i*2, 1))
		grid.text(rownames(plot_layout)[i], x=unit(0.5, "npc"), vjust=0.5,hjust=0.5, gp = gpar(fontface = "bold", fontsize=24), rot=90)
		upViewport()
	}
	
	for(i in 1:ncol(plot_layout)) {
		#plot titles
		pushViewport(vplayout(1, i*2))
		grid.text(colnames(plot_layout)[i], x=unit(0.5, "npc"), vjust=0.5,hjust=0.5, gp = gpar(fontface = "bold", fontsize=24), rot=0)
		upViewport()
		
		if(is.list(colscales[[i]])) {
			cbreaks=colscales[[i]]$cbreaks
			col_breaks=colscales[[i]]$col_breaks
			scale_cols=colscales[[i]]$scale
			legend = colscales[[i]]$legend
		} else {
			if(colscales[[i]]=="autosym") {
				cbreaks=seq(-max_val_per_col[i],max_val_per_col[i],length.out=63)
				col_breaks=c(-(max_val_per_col[i]+5), cbreaks, max_val_per_col[i]+5)
				scale_cols=greenred(64)
				legend=seq(-max_val_per_col[i], max_val_per_col[i], length.out=5)
			} else if(colscales[[i]]=="autopos") {
			  cbreaks=seq(0,max_val_per_col[i],length.out=63)
			  col_breaks=c(-5, cbreaks, max_val_per_col[i]+5)
			  scale_cols=greenred(64)
			  legend=seq(-max_val_per_col[i], max_val_per_col[i], length.out=5)
      } else {
				cbreaks=seq(0,max_val_per_col[i],length.out=63)
				col_breaks=c(-15, cbreaks, max_val_per_col[i]+5)
				scale_cols=colorpanel(64, low="#000000", high="#FF0000")
				legend=seq(0, max_val_per_col[i], length.out=4)
			}
		}
				
		for(j in 1:nrow(plot_layout)) {
			num=plot_layout[j,i]
			if(num > 0) {
				tdata=datalist[[num]]$data
				pushViewport(vplayout(j*2, i*2, just=c("left", "top")))
				make_single_plot(tdata, datalist[[num]]$size, showcolns=(j==TRUE), 
                         col_breaks, cbreaks, scale_cols, cluster, flip_data=flip_data,
				                 show_dend_col=show_dend_col, show_dend_row=show_dend_row)
				upViewport()
			}
		}
		pushViewport(vplayout(2,i*2+1, just=c("right", "bottom")))
		names(legend)=legend
		draw_legend(color=scale_cols, ht=0.95, cbreaks, legend, fontsize = 18, family="Sans")
		upViewport()
	}
	
	upViewport()
	dev.off()
}

make_single_plot=function(grp_data, size_array=0, showcolns=TRUE, col_breaks, cbreaks, scale_cols, cluster=FALSE, show_dend_col=FALSE, show_dend_row=FALSE, flip_data=FALSE) {
	n_subp=length(grp_data)
	
	if(!is.list(size_array)) {
		sizes=matrix(nrow=0, ncol=3)
		for(i in 1:n_subp) {
			tdata=grp_data[[i]]$data
			if(flip_data) {
			  tdata=t(tdata)
			}
			panel_size=get_size(t(data.matrix(tdata)), sidelabel=grp_data[[i]]$name, filename="test_lvl2.pdf", showcolns=(i==TRUE & showcolns))
			sizes=rbind(sizes, c(panel_size[[1]], panel_size[[2]], panel_size[[3]]))
		}
	
		min_width=max(sizes[,1])+0.5
		min_height=sum(sizes[,2])+1
		min_main=max(sizes[,3])
		return(list(panel=sizes, min_width=min_width, min_height=min_height -0.05* n_subp, main_width=min_main))
	}
	panel_sizes=size_array$panel;
	preset_main_width=size_array$main_width
	p_ht=vector();
	break_ht=unit(-0.001, "inches")
	p_ht=unit(panel_sizes[1,2], "inches")
	hts=unit.c(p_ht, break_ht)
	if(n_subp > 1) {
		for(i in 2:n_subp) {
			p_ht=unit(panel_sizes[i,2], "inches")
			hts=unit.c(hts, p_ht, break_ht)
		}
	}
	
	p_w=unit(size_array$min_width, "inches")
	
	pushViewport(viewport(layout = grid.layout(nrow = n_subp*2 - 1, ncol = 1, widths=p_w, heights=hts), gp = gpar(fontsize=18)))

	for(i in 1:n_subp) {
		tdata=grp_data[[i]]$data
		if(flip_data) {
		  tdata=t(tdata)
		}
    pushViewport(vplayout((i-1)*2+1, 1, just=c("left", "top")))
		pheatmap(t(data.matrix(tdata)), col=scale_cols, cluster_rows=cluster, 
			cluster_cols=cluster, family="Sans", fontsize=18, cellwidth=36,
			cellheight=20, show_colnames=(i==TRUE & showcolns), preset_main_width=preset_main_width,
			breaks=col_breaks, scale="none", legend=FALSE, main=grp_data[[i]]$name,
			show_dend_col=show_dend_col, show_dend_row=show_dend_row)
		upViewport()
	}

	upViewport()
}

#plots all things with broken up plots to make large reference plots
#	pdata-ordered pdata
#	num_breaks-number of breaks to make
make_break_plots=function(pdata, num_breaks, filename) {
	brk_len=ceiling(ncol(pdata)/num_breaks)
	break_pos=seq(0,ncol(pdata), brk_len)
	if(length(break_pos) < num_breaks) {
		break_pos=c(break_pos, ncol(pdata))
	}

	panel_sizes=vector()
	for(i in 1:num_breaks) {
		tdata=pdata[, (break_pos[i]+1):break_pos[i+1]]
		panel_size=get_size(t(data.matrix(tdata)), filename="test.pdf", showcolns=TRUE)
		panel_sizes=c(panel_sizes, panel_size[[1]], panel_size[[2]])
	}

	min_width=sum(panel_sizes[(1:num_breaks)*2-1])
	min_height=max(panel_sizes[(1:num_breaks)*2])
	pdf(paste(filename, ".pdf", sep=""), width=min_width+1*num_breaks, height=min_height+2)

	grid.newpage()
	p_width=vector();
	break_width=unit(0.4, "inches")
	p_width=unit(panel_sizes[1], "inches")
	widths=unit.c(p_width, break_width)
	for(i in 2:num_breaks) {
		p_width=unit(panel_sizes[i*2-1], "inches")
		widths=unit.c(widths, p_width, break_width)
	}

	p_height=unit(min_height, "inches")
	heights=p_height

	pushViewport(viewport(layout = grid.layout(nrow = 1, ncol = num_breaks*2, widths=widths, heights=heights), gp = gpar(fontsize=18)))

	for(i in 1:num_breaks) {
	tdata=pdata[, (break_pos[i]+1):break_pos[i+1]]
	pushViewport(vplayout(1, (i-1)*2+1, just=c("left", "top")))
	pheatmap(t(data.matrix(tdata)), col=greenred(64), cluster_rows=FALSE, 
			cluster_cols=FALSE, family="Sans", fontsize=18, cellwidth=36,
			cellheight=20, show_colnames=TRUE,
			breaks=col_breaks, scale="none", legend=FALSE)
	upViewport()
	}

	pushViewport(vplayout(1,num_breaks*2, just="bottom"))
	legend = c(-4,-2,0 ,2,4)
	names(legend)=legend
	draw_legend(color=greenred(64), ht=0.95, cbreaks, legend, fontsize = 18, family="Sans")
	upViewport()

	dev.off()
	upViewport()
}

