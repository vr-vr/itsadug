#' Plot difference curve based on model predictions.
#' 
#' @export
#' @aliases plotDiff
#' @param model A GAMM model, resulting from the functions
#' \code{\link[mgcv]{gam}} or \code{\link[mgcv]{bam}}.
#' @param  view Name of continuous predictor that should be plotted on the x-
#' axis.
#' @param comp Named list with the grouping predictor (categorical variable)
#' and the 2 levels to calculate the difference for.
#' @param cond A named list of the values to use for the predictor terms. 
#' Variables omitted from this list will have the closest observed value to 
#' the median for continuous variables, or the reference level for factors. 
#' @param plotCI Logical: whether or not to plot confidence intervals. 
#' Default is TRUE.
#' @param f A number to scale the standard error. Defaults to 1.96, resulting 
#' in 95\% confidence intervals. For 99\% confidence intervals use a value of 
#' 2.58.
#' @param rm.ranef Logical: whether or not to remove random effects. 
#' Default is FALSE. Alternatively a string (or vector of strings) with the 
#' name of the random effect(s) to remove.
#' @param eegAxis Logical: whether or not to reverse the y-axis, plotting 
#' negative values upwards. Default is FALSE.
#' @param col Line color. Shading color is derived from line color.
#' @param shade Logical: plot shaded confidence interval (TRUE) 
#' or dashed lines that indicate confidence region (FALSE).
#' @param ylim Range of y-axis. If not specified, the function automatically 
#' generates an appropriate y-axis.
#' @param main Text string, alternative title for plot.
#' @param xlab Text string, alternative label for x-axis. 
#' @param ylab Text string, alternative label for y-axis. 
#' @param n.grid Number of data points sampled as predictions. Defaults to 100.
#' @param add Logical: whether or not to add the line to an existing plot. 
#' Default is FALSE.
#' When no plot window is available and \code{add=TRUE}, 
#' the function will generate an error.
#' @param plot Logical: whether or not to plot the difference. If FALSE, then 
#' the output is returned as a list, with the estimated difference 
#' (\code{est}) and the standard error over the estimate (\code{se.est}) and 
#' the x-values (\code{x}). Default is TRUE.
#' @param hide.label Logical: whether or not to hide the label 
#' (i.e., "difference"). Default is FALSE.
#' @param print.summary Logical: whether or not to print the summary. 
#' Default set to the print info messages option 
#' (see \code{\link{infoMessages}}).
#' @param ... Optional arguments for plot.
#' @return If the result is not being plotted, a list is 
#' returned with the estimated difference (\code{est}) and the standard error 
#' over the estimate (\code{se}) and the x-values (\code{x}) is returned.
#' @author Martijn Wieling, Jacolien van Rij
#' @examples
#' data(simdat)
#' \dontrun{
#' m1 <- bam(Y ~ Group + te(Time, Trial, by=Group),
#'     data=simdat)
#' plot_diff(m1, view='Time', comp=list(Group=c("Children", "Adults")))
#' # Reversed y-axis (for EEG data):
#' plot_diff(m1, view='Time', comp=list(Group=c("Children", "Adults")), 
#'     eegAxis=TRUE)
#' # retrieving plot values...
#' out <- plot_diff(m1, view='Time', comp=list(Group=c("Children", "Adults")), 
#'    plot=FALSE)
#' #... which might be used for indicating differences:
#' x <- find_difference(out$est, out$se, f=1.96, xVals=out$xVals)
#' # add lines:
#' arrows(x0=x$start, x1=x$end, y0=0, y1=0,code=3, length=.1, col='red')
#' }
#'
#' @family functions for interpreting nonlinear effects

plot_diff <- function(model, view, comp, cond=NULL, plotCI=TRUE, f=1.96, 
	eegAxis=FALSE, col="black", shade=TRUE, n.grid=100, add=FALSE,
	print.summary=getOption('itsadug_print'), plot=TRUE, rm.ranef=NULL,
	main=NULL, ylab=NULL, xlab=NULL, ylim=NULL, 
	hide.label=FALSE, ...) { 

	dat = model$model

	xvar <- NULL
	by_predictor <- NULL

	# check view
	if(length(view) > 1){
		warning("Only first element of 'view' is being used. Use plot_diff2 for plotting difference surfaces.")
	}else{
		xvar <- view[1]

		if(xvar %in% names(cond)){
			warning(sprintf('Predictor %s specified in view and cond. Values in cond being used, rather than the whole range of %s.', xvar))
		}else{
			cond[[xvar]] <- seq(min(na.exclude(dat[,xvar])), max(na.exclude(dat[,xvar])), length=n.grid)
		}
	}

	newd <- c()
	newd <- get_difference(model, comp=comp, cond=cond, 
		print.summary=print.summary, rm.ranef=rm.ranef, f=f)


	if (is.null(main)) {
		levels1 <- paste(sapply(comp, function(x) x[1]), collapse='.')
		levels2 <- paste(sapply(comp, function(x) x[2]), collapse='.')
		main = sprintf('Difference between %s and %s', levels1, levels2)
	} 
	if(is.null(ylab)) {
		ylab = sprintf("Est. difference in %s", as.character(model$formula[[2]]))
	}
	if(is.null(xlab)) {
		xlab = xvar
	}
	if(is.null(ylim)){
		ylim <- range(newd$difference)
		if (plotCI) { 
			ylim <- with(newd, range(c(difference+CI, difference-CI)))
		}
	}

	out <- list(est=newd$difference)
	if(plotCI){
		out[['se']]<- newd$CI
	}
	out[['xVals']] <- newd[,xvar]
	out[['f']] <- f
	out[['comp']] = comp

	if(plot==TRUE){
		if(add==FALSE){
			emptyPlot(range(newd[,xvar]), ylim, 
				main=main, xlab=xlab, ylab=ylab, h0=0,
				eegAxis=eegAxis, ...)
			if(hide.label==FALSE){
	            addlabel = "difference"
	            if(!is.null(rm.ranef)){
	                if(rm.ranef !=FALSE){
	                    addlabel = paste(addlabel, "excl. random", sep=", ")
	                }
	            }
	            mtext(addlabel, side=4, line=0, adj=0, 
	                cex=.75, col='gray35', xpd=TRUE)
	        }        			
		}
		if(plotCI==TRUE){
			plot_error(newd[,xvar], newd$difference, newd$CI, shade=shade, col=col, ...)
		}else{
			lines(newd[,xvar], newd$difference, col=col, ...)
		}
		invisible(out)

	}else{
		return(out)
	}
}






#' Plot difference surface based on model predictions.
#' 
#' @export
#' @aliases plotDiff2D
#' @param model A GAMM model, resulting from the functions
#' \code{\link[mgcv]{gam}} or \code{\link[mgcv]{bam}}.
#' @param  view Name of continuous predictors that should be plotted on the x-
#'  and y-axes. Vector of two values.
#' @param comp Named list with the grouping predictor (categorical variable)
#' and the 2 levels to calculate the difference for.
#' @param cond Named list of the values to use for the other predictor terms 
#' (not in view). 
#' @param color Colorpalette
#' @param nCol Range of colors of background of contour plot.
#' @param col Line color.
#' @param add.color.legend Logical: whether or not to add a color legend. 
#' Default is TRUE. If FALSE (omitted), one could use the function
#' \code{\link{gradientLegend}} to add a legend manually at any position.
#' @param plotCI Logical: whether or not to plot confidence intervals.
#' @param f A number to scale the standard error. Defaults to 1.96, resulting 
#' in 95\% confidence intervals. For 99\% confidence intervals use a value of 
#' 2.58.
#' @param n.grid Resolution.
#' @param nlevels Levels of contour lines.
#' @param zlim A two item array giving the lower and upper limits for the z-
#' axis scale. NULL to choose automatically.
#' @param xlim A two item array giving the lower and upper limits for the x-
#' axis scale. NULL to choose automatically.
#' @param ylim A two item array giving the lower and upper limits for the y-
#' axis scale. NULL to choose automatically.
#' @param main Title of plot.
#' @param xlab Label x-axis.
#' @param ylab Label y-axis.
#' @param rm.ranef Logical: whether or not to remove random effects. 
#' Default is FALSE. Alternatively a string (or vector of strings) with the 
#' name of the random effect(s) to remove.
#' @param print.summary Logical: whether or not to print a summary.
#' Default set to the print info messages option 
#' (see \code{\link{infoMessages}}).
#' @param hide.label Logical: whether or not to hide the label 
#' (i.e., "difference"). Default is FALSE.
#' @param ... Optional arguments for \code{\link{plotsurface}}.
#' @return If the result is not being plotted, a list is 
#' returned with the estimated difference (\code{est}) and the standard error 
#' over the estimate (\code{se.est}) and the x-values (\code{x}) is returned.
#' @author Martijn Wieling, reimplemented by Jacolien van Rij
#'
#' @examples
#' data(simdat)
#' \dontrun{
#' m1 <- bam(Y ~ Group + te(Time, Trial, by=Group),
#'     data=simdat)
#' plot_diff2(m1, view=c('Time', 'Trial'), 
#'     comp=list(Group=c("Children", "Adults")))
#' }
#' # see the vignette for examples:
#' vignette("plotfunctions", package="itsadug")
#' @family functions for interpreting nonlinear effects


# plots differences in 2D plot
plot_diff2 <- function(model, view, comp, cond=NULL, 
	color='topo', nCol=100, col=NULL, add.color.legend=TRUE,
	plotCI=FALSE, f=1.96, n.grid=30, nlevels=10, 
	zlim=NULL, xlim=NULL, ylim=NULL, 
	main=NULL, xlab=NULL, ylab=NULL,
	rm.ranef=NULL,
	hide.label=FALSE,
	print.summary=getOption('itsadug_print'), ...) { 

	dat = model$model

	xvar <- NULL
	yvar <- NULL
	by_predictor <- NULL

	# check view
	if(length(view) < 2){
		stop('Provide predictors for x- and y-axes in view.')
	}else{
		xvar <- view[1]
		yvar <- view[2]

		if(xvar %in% names(cond)){
			warning(sprintf('Predictor %s specified in view and cond. Values in cond being used, rather than the whole range of %s.', xvar, xvar))
		}else{
			cond[[xvar]] <- seq(min(na.exclude(dat[,xvar])), max(na.exclude(dat[,xvar])), length=n.grid)
			if(!is.null(xlim)){
        		if(length(xlim) != 2){
            		warning("Invalid xlim values specified. Argument xlim is being ignored.")
        		}else{ 
            		cond[[xvar]] <- seq(xlim[1], xlim[2], length=n.grid)
        		}
    		}
		}
		if(yvar %in% names(cond)){
			warning(sprintf('Predictor %s specified in view and cond. Values in cond being used, rather than the whole range of %s.', yvar, yvar))
			cond[[yvar]] <- NULL
		}else{
			cond[[yvar]] <- seq(min(na.exclude(dat[,yvar])), max(na.exclude(dat[,yvar])), length=n.grid)
			if(!is.null(ylim)){
        		if(length(ylim) != 2){
            		warning("Invalid ylim values specified. Argument ylim is being ignored.")
        		}else{ 
            		cond[[yvar]] <- seq(ylim[1], ylim[2], length=n.grid)
        		}
    		}
		}
	}


	newd <- c()
	newd <- get_difference(model, comp=comp, cond=cond, 
		print.summary=print.summary, rm.ranef=rm.ranef, f=f)

	if (is.null(main)) {
		levels1 <- paste(sapply(comp, function(x) x[1]), collapse='.')
		levels2 <- paste(sapply(comp, function(x) x[2]), collapse='.')
		main = sprintf('Difference between %s and %s', levels1, levels2)
	} 
	if(is.null(ylab)) {
		ylab = view[2]
	}
	if(is.null(xlab)) {
		xlab = view[1]
	}

	if(plotCI){
		p <- plotsurface(newd, view=view, predictor="difference", valCI='CI',
			main=main, xlab=xlab, ylab=ylab, 
			zlim=zlim, 
			col=col, color=color, nCol=nCol, add.color.legend=add.color.legend,
			nlevels=nlevels, ...)
        if(hide.label==FALSE){
            addlabel = "difference"
            if(!is.null(rm.ranef)){
                if(rm.ranef !=FALSE){
                    addlabel = paste(addlabel, "excl. random", sep=", ")
                }
            }
            mtext(addlabel, side=4, line=0, adj=0, 
                cex=.75, col='gray35', xpd=TRUE)
        }     		
	}else{
		p <- plotsurface(newd, view=view, predictor="difference", 
			main=main, xlab=xlab, ylab=ylab, 
			zlim=zlim, 
			col=col, color=color, nCol=nCol, add.color.legend=add.color.legend,
			nlevels=nlevels, ...)	
        if(hide.label==FALSE){
            addlabel = "difference"
            if(!is.null(rm.ranef)){
                if(rm.ranef !=FALSE){
                    addlabel = paste(addlabel, "excl. random", sep=", ")
                }
            }
            mtext(addlabel, side=4, line=0, adj=0, 
                cex=.75, col='gray35', xpd=TRUE)
        }  	
	}
	
	p[['zlim']] <- zlim
	invisible(p)
}
