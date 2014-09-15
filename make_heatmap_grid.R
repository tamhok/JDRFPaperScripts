lo = function(rown, coln, nrow, ncol, cellheight = NA, cellwidth = NA, treeheight_col, treeheight_row, legend, annotation, annotation_colors, annotation_legend, main, fontsize, fontsize_row, fontsize_col, ...){
	# Get height of colnames and length of rownames
	if(!is.null(coln[1])){
		longest_coln = which.max(nchar(coln))
		gp = list(fontsize = fontsize_col, ...)
		coln_height = unit(1.1, "grobheight", textGrob(coln[longest_coln], rot = 90, gp = do.call(gpar, gp)))
	}
	else{
		coln_height = unit(5, "bigpts")
	}
	
	if(!is.null(rown[1])){
		longest_rown = which.max(nchar(rown))
		gp = list(fontsize = fontsize_row, ...)
		rown_width = unit(1.2, "grobwidth", textGrob(rown[longest_rown], gp = do.call(gpar, gp)))
	}
	else{
		rown_width = unit(5, "bigpts")
	}
	
	gp = list(fontsize = fontsize, ...)
	# Legend position
	if(!is.na(legend[1])){
		longest_break = which.max(nchar(names(legend)))
		longest_break = unit(1.1, "grobwidth", textGrob(as.character(names(legend))[longest_break], gp = do.call(gpar, gp)))
		title_length = unit(1.1, "grobwidth", textGrob("Scale", gp = gpar(fontface = "bold", ...)))
		legend_width = unit(12, "bigpts") + longest_break * 1.2
		legend_width = max(title_length, legend_width)
	}
	else{
		legend_width = unit(0, "bigpts")
	}
	
	# Set main title height
	if(is.na(main)){
		main_height = unit(0, "npc")
	}
	else{
		main_height = unit(1.5, "grobheight", textGrob(main, gp = gpar(fontsize = 1.3 * fontsize, ...)))
	}
	
	# Column annotations
	if(!is.na(annotation[[1]][1])){
		# Column annotation height 
		annot_height = unit(ncol(annotation) * (8 + 2) + 2, "bigpts")
		# Width of the correponding legend
		longest_ann = which.max(nchar(as.matrix(annotation)))
		annot_legend_width = unit(1.2, "grobwidth", textGrob(as.matrix(annotation)[longest_ann], gp = gpar(...))) + unit(12, "bigpts")
		if(!annotation_legend){
			annot_legend_width = unit(0, "npc")
		}
	}
	else{
		annot_height = unit(0, "bigpts")
		annot_legend_width = unit(0, "bigpts")
	}
	
	# Tree height
	treeheight_col = unit(treeheight_col, "bigpts") + unit(5, "bigpts")
	treeheight_row = unit(treeheight_row, "bigpts") + unit(5, "bigpts") 
	
	# Set cell sizes
	if(is.na(cellwidth)){
		matwidth = unit(1, "npc") - rown_width - legend_width - treeheight_row - annot_legend_width
	}
	else{
		matwidth = unit(cellwidth * ncol, "bigpts")
	}
	
	if(is.na(cellheight)){
		matheight = unit(1, "npc") - main_height - coln_height - treeheight_col - annot_height
	}
	else{
		matheight = unit(cellheight * nrow, "bigpts")
	}	
	
	
	# Produce layout()
	pushViewport(viewport(layout = grid.layout(nrow = 5, ncol = 5, widths = unit.c(treeheight_row, matwidth, rown_width, legend_width, annot_legend_width), heights = unit.c(main_height, treeheight_col, annot_height, matheight, coln_height)), gp = do.call(gpar, gp)))
	
	# Get cell dimensions
	pushViewport(vplayout(4, 2))
	cellwidth = convertWidth(unit(0:1, "npc"), "bigpts", valueOnly = T)[2] / ncol
	cellheight = convertHeight(unit(0:1, "npc"), "bigpts", valueOnly = T)[2] / nrow
	upViewport()
	
	# Return minimal cell dimension in bigpts to decide if borders are drawn
	mindim = min(cellwidth, cellheight) 
	return(mindim)
}
