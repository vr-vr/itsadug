#' Visualization of EEG topo maps.
#'
#' @export
#'
#' @param model A gam object, produced by \code{\link[mgcv]{gam}} or 
#' \code{\link[mgcv]{bam}}.
#' @param view A two-value vector containing the names of the two main effect 
#' terms to be displayed on the x and y dimensions of the plot. Note that 
#' variables coerced to factors in the model formula won't work as view 
#' variables.
#' @param el.pos A list with X and Y positions and Electrodes, which are used
#' for fitting the model.
#' @param fun Text string, "fvisgam", "pvisgam", or "plot_diff2" signalling 
#' which function to use for plotting.
#' @param add.color.legend Logical: whether or not to add a color legend. 
#' Default is TRUE. If FALSE (omitted), one could use the function
#' \code{\link{gradientLegend}} to add a legend manually at any position.
#' @param size Size in inch of plot window.
#' @param n.grid  The number of grid nodes in each direction used for 
#' calculating the plotted surface. 
#' @param col The colors for the background of the plot.
#' @param color The color scheme to use for plots. One of "topo", "heat", 
#' "cm", "terrain", "gray" or "bw". 
#' @param pch The type of points as indications for the electrode positions. 
#' The value NA will suppress the plotting of electrode positions.
#' @param bg The background color of the points.
#' @param xlab Label x-axis. Default excluded.
#' @param ylab Label y-axis. Default excluded.
#' @param ... other options to pass on to \code{\link{fvisgam}},  
#' \code{\link{pvisgam}}, or \code{\link{plot_diff2}}. 
#' @section Notes:
#' X-positions of electrodes should have lower values for electrodes on the
#' left hemisphere (e.g. T7) than for electrodes on the right 
#' hemisphere.
#' Y-positions of electrodes should have lower values for electrodes at the 
#' back of the head than for the frontal electrodes.
#' @author Jacolien van Rij
#' @examples
#'
#' data(eeg)
#'
#' \dontrun{
#' # simple GAMM model:
#' m1 <- gam(Ampl ~ te(Time, X, Y, k=c(10,5,5), 
#'     d=c(1,2)), data=eeg)
#'
#' # topo plot, by default uses fvisgam 
#' # and automatically selects a timestamp (270ms):
#' plot_topo(m1, view=c("X", "Y"))
#' 
#' # add electrodes:
#' electrodes <- eeg[,c('X','Y','Electrode')]
#' electrodes <- as.list( electrodes[!duplicated(electrodes),] )
#' plot_topo(m1, view=c("X", "Y"), el.pos=electrodes)
#' 
#' # some formatting options:
#' plot_topo(m1, view=c("X", "Y"), el.pos=electrodes,
#'     main="Topo plot", zlim=c(-.5,.5), 
#'     pch=15, col='red', color='terrain')
#' 
#' # plotting more than one panel only works if 
#' # each figure region is a square:
#' dev.new(width=12, height=4) 
#' par(mfrow=c(1,3))
#'
#' for(i in c(100, 200, 300)){
#'     # make sure to keep zlim constant:
#'	   plot_topo(m1, view=c('X', 'Y'), zlim=c(-.5, .5), 
#'     cond=list(Time=i), el.pos=electrodes,
#'     main=i)
#' }
#' 
#' dev.new(width=12, height=4) 
#' par(mfrow=c(1,3), cex=1.1)
#' # The three different functions for plotting:
#' plot_topo(m1, view=c('X', 'Y'), zlim=c(-.5, .5), 
#'     el.pos=electrodes,
#'     fun='fvisgam', main='fvisgam', 
#'     cond=list(Time=200), rm.ranef=TRUE)
#' plot_topo(m1, view=c('X', 'Y'), zlim=c(-.5, .5), 
#'     el.pos=electrodes, select=1,
#'     fun='pvisgam', main='pvisgam', 
#'     cond=list(Time=200))
#' plot_topo(m1, view=c('X', 'Y'), zlim=c(-.5, .5), 
#'     el.pos=electrodes, comp=list(Time=c(300,100)),
#'     fun='plot_diff2', main='plot_diff2', 
#'     plotCI=TRUE)
#' 
#' # Add labels:
#' plot_topo(m1, view=c('X', 'Y'), zlim=c(-.5, .5), 
#'     fun='fvisgam', main='', 
#'     cond=list(Time=200), add.color.legend=FALSE)
#' text(electrodes[['X']], electrodes[['Y']], 
#'     labels=electrodes[['Electrode']], cex=.75, 
#'     xpd=TRUE)
#' }
#' @family functions for interpreting nonlinear effects

plot_topo <- function(model, view, el.pos=NULL, fun='fvisgam', 
	add.color.legend=TRUE, 
	size=5, n.grid=100, col=1, pch=21, bg=alpha(1), 
	color='topo', xlab="", ylab="", ...){

	# save old settings
	oldpar =  list(mai=par()$mai, pin=par()$pin,
		xaxt = par()$xaxt, yaxt = par()$yaxt,
		bty = par()$bty)

	if( !is.null( names(dev.list()) )){
		if(round( par()$fin[1], 4) != round( par()$fin[2], 4)){
			dev.new(width=size, height=size, noRStudioGD = TRUE)
		}
	}

	par(pin=c(size,size), mai=rep(.5,4), xaxt='n', yaxt='n', bty='n')

	# range X
	r.x <- range(model$model[,view[1]])
	# range Y
	r.y <- range(model$model[,view[2]])
	# center
	center <- c(r.x[1]+(r.x[2]-r.x[1])/2, r.y[1]+(r.y[2]-r.y[1])/2)
	# radius
	radius <- max(c((r.x[2]-r.x[1])/2, (r.y[2]-r.y[1])/2))

	if(!is.null(el.pos)){
		# range X
		r.x <- range(c(model$model[,view[1]], el.pos[['X']]))
		# range Y
		r.y <- range(c(model$model[,view[2]], el.pos[['Y']]))
		# center
		center <- c(r.x[1]+(r.x[2]-r.x[1])/2, r.y[1]+(r.y[2]-r.y[1])/2)
		# radius
		radius <- max(c((r.x[2]-r.x[1])/2, (r.y[2]-r.y[1])/2))
	}

	xlim = c(center[1]-radius-0.1*radius, center[1]+radius+0.1*radius)
	ylim = c(center[2]-radius-0.1*radius, center[2]+radius+0.1*radius)	

	# plot effects
	if(fun=='fvisgam'){
		pp <- fvisgam(model, view, add.color.legend=FALSE, n.grid=n.grid, 
			xlab=xlab, ylab=ylab, 
			xlim=xlim, ylim=ylim, color=color, ...)
	}else if (fun=='pvisgam'){
		pp <- pvisgam(model, view, add.color.legend=FALSE, n.grid=n.grid, 
			xlab=xlab, ylab=ylab, 
			xlim=xlim, ylim=ylim, color=color, ...)
	}else if (fun=='plot_diff2'){
		pp <- plot_diff2(model, view, add.color.legend=FALSE, n.grid=n.grid, 
			xlab=xlab, ylab=ylab, 
			xlim=xlim, ylim=ylim, color=color, ...)
	}else{
		stop("Function unknown.")
	}
	par(bty='o')
	box(col='white', lwd=1)

	# add circle
	gfc <- getFigCoords()

	im <- expand.grid(x=seq(gfc[1], gfc[2], length=n.grid),
		y=seq(gfc[3], gfc[4], length=n.grid))
	
	im$val <- NA
	im$val <- ifelse( sqrt(im$x^2+im$y^2 ) >= (radius+.1*radius),
		alpha('white', f=1), alpha('white', f=0))

	mat <- matrix(im$val, byrow=TRUE, ncol=n.grid)
	rasterImage(mat, xleft=gfc[1], ybottom=gfc[3],
		xright=gfc[2], ytop=gfc[4])
	points(el.pos[['X']], el.pos[['Y']], col=col, pch=pch, bg=bg, xpd=TRUE)
	text(center[1], center[2]+radius, labels="^", pos=3, xpd=TRUE)
	if(add.color.legend==TRUE){
		gradientLegend(round(pp$zlim,3), color=color, pos=0.125, side=1, inside=FALSE)
	}
	
	for(i in names(oldpar)){
		eval(parse(text=sprintf("par(%s=%s)", i, 
			ifelse(is.character(oldpar[[i]]), 
				sprintf("'%s'", oldpar[[i]]), 
				sprintf("c(%s)", paste(oldpar[[i]], collapse=","))))))
	}
}