library(gplots)

expr=function(val) {
	b=parse(text=paste("expression(", val, ")"))
	return(eval(b));
}

plot_pie=function (x, labels = names(x), draw_borders, outer_labels, out_lbl_x, edges = 200, radius = 0.8, clockwise = FALSE, 
    init.angle = if (clockwise) 90 else 10, density = NULL, angle = 45, 
    col = NULL, border = NULL, lty = NULL, main = NULL, ...) 
{
   
   if (!is.numeric(x) || any(is.na(x) | x < 0)) 
        stop("'x' values must be positive.")
    if (is.null(labels)) 
        labels <- as.character(seq_along(x))
	if(is.null(draw_borders))
		draw_borders=rep(TRUE, length(x))
    else labels <- as.graphicsAnnot(labels)
    x <- c(0, cumsum(x)/sum(x))
    dx <- diff(x)
    nx <- length(dx)
    plot.new()
    pin <- par("pin")
    xlim <- ylim <- c(-1, 1)
    if (pin[1L] > pin[2L]) 
        xlim <- (pin[1L]/pin[2L]) * xlim
    else ylim <- (pin[2L]/pin[1L]) * ylim
    dev.hold()
    on.exit(dev.flush())
    plot.window(xlim, ylim, "", asp = 1)
    if (is.null(col)) 
        col <- if (is.null(density)) 
            c("white", "lightblue", "mistyrose", "lightcyan", 
                "lavender", "cornsilk")
        else par("fg")
    if (!is.null(col)) 
        col <- rep_len(col, nx)
    if (!is.null(border)) 
        border <- rep_len(border, nx)
    if (!is.null(lty)) 
        lty <- rep_len(lty, nx)
    angle <- rep(angle, nx)
    if (!is.null(density)) 
        density <- rep_len(density, nx)
    twopi <- if (clockwise) 
        -2 * pi
    else 2 * pi
    t2xy <- function(t, rad=radius) {
        t2p <- twopi * t + init.angle * pi/180
        list(x = rad * cos(t2p), y = rad * sin(t2p))
    }
	for (i in 1L:nx) {
        n <- max(2, floor(edges * dx[i]))
        P <- t2xy(seq.int(x[i], x[i + 1], length.out = n))
        if(draw_borders[i]) {
		polygon(c(P$x, 0), c(P$y, 0), density = density[i], angle = angle[i], 
            border = border[i], col = col[i], lty = lty[i])
		}
        P <- t2xy(mean(x[i + 0:1]))
        lab <- as.character(labels[i])
        if (!is.na(lab) && nzchar(lab)) {
            if(draw_borders[i]) {
			lines(c(1, 1.05) * P$x, c(1, 1.05) * P$y)
            text(1.1 * P$x, 1.15 * P$y, lab, xpd = TRUE, 
                adj = ifelse(P$x < 0, 1, 0), cex=1.0, ...)
			}
        }
    }
	if(length(outer_labels)>1) {
		o = c(0, cumsum(out_lbl_x)/sum(out_lbl_x));
		do=diff(o)
		no=length(do)
		for (i in 1:length(outer_labels)) {
			if(!is.null(outer_labels[i]) & !is.na(outer_labels[i])) {
				P <- t2xy(mean(o[i + 0:1]), radius*1.5)
				lab <- as.character(outer_labels[i])
				if (!is.na(lab) && nzchar(lab)) {
					text(1.1 * P$x, 1.1 * P$y, expr(outer_labels[i]), xpd = TRUE, 
						adj = ifelse(P$x < 0, 1, 0), cex=2,...)
					lines(c(1, 1.05) * P$x, c(1, 1.05) * P$y)
				}
				n <- max(2, floor(edges * do[i]))
				P <- t2xy(seq.int(o[i], o[i+1], length.out = n), radius*1.5)
				segments(x0=P$x[-n], y0=P$y[-n], x1=P$x[-1], y1=P$y[-1], density = density[i],
				border = border[i], lty = lty[i])
				
			}
		}
	}
	
    title(main = main, cex.main=3,...)
    invisible(NULL)
}

to_colors=function(val, vmax, vmin, lb) {
	res=100;
	cscl=c("#FFFFFF", greenred(res+1));
	im=round(val*res/(vmax-vmin)+res/2)+2;
	im[im > res+2] = res+1;
	im[im < 2] = 2;
	im[val <= lb | is.na(val)] = 1;
	
	return(cscl[im]);
}

#Calls the draw pie chart function. Sorts data into groups.
#Need to pass it a matrix with 
make_pie=function(groups, title=NA) {
	ngrps=levels(as.factor(groups[,2]));
	ngens=length(groups[,1]);

	cols=vector();
	labl=vector();
	rvals=vector();
	draw_b=vector();
	out_lbls=vector();
	out_lbl_pos=vector();
	for (i in 1:length(ngrps)) {		
		sg=which(groups[,2]==ngrps[i]);
		this.data=groups[sg,4];
		genes=groups[sg,3];
		d.order=order(groups[sg,5]);		
		out_lbls=c(out_lbls, as.character(ngrps[i]), NA);
		this.data=this.data[d.order];
		
		cols=c(cols, to_colors(this.data, 4,-4,-8), "#FFFFFF");
		labl=c(labl, genes[d.order], " "); 
		rvals=c(rvals, rep(1, length(d.order)), 0.7);
		draw_b=c(draw_b, rep(TRUE, length(d.order)));
		out_lbl_pos=c(out_lbl_pos, length(d.order), 0.7)
		draw_b=c(draw_b, FALSE);
	}
	plot_pie(rvals, col=cols, draw_border=draw_b, lab=labl, 
	outer_labels=out_lbls, out_lbl_x=out_lbl_pos, radius=0.7, main=title);
	
}