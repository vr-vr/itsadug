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
#' @param eegAxis Logical: whether or not to reverse the y-axis, plotting 
#' negative values upwards. Default is FALSE.
#' @param ylim Range of y-axis. If not specified, the function automatically 
#' generates an appropriate y-axis.
#' @param ylab Text string, alternative label for y-axis. 
#' @param main Text string, alternative title for plot.
#' @param plot Logical: whether or not to plot the difference. If FALSE, then 
#' the output is returned as a list, with the estimated difference 
#' (\code{est}) and the standard error over the estimate (\code{se.est}) and 
#' the x-values (\code{x}). Default is TRUE.
#' @param ... Optional arguments for plot.
#' @return If the result is not being plotted, a list is 
#' returned with the estimated difference (\code{est}) and the standard error 
#' over the estimate (\code{se.est}) and the x-values (\code{x}) is returned.
#' @author Martijn Wieling
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
#' # see the vignette for examples:
#' vignette("plotfunctions", package="itsadug")
#'
#' @family functions for interpreting nonlinear effects

# plots difference curve
# set binary=T if you have a by-variable occurring multiple times 
plot_diff <- function(model, view, comp, eegAxis=F, 
	ylim=NULL, ylab=NULL, main=NULL, plot=TRUE, ...) { 
	dat = model$model

	by_predictor <- NULL
	comp_levels <- NULL
	if(length(view) > 1){
		warning('Only first element of view will be used.')
		view <- view[1]
	}
	if(names(comp)[1] %in% colnames(dat)){
		if(length(comp[[1]]) < 2){
			stop(sprintf('Provide two levels for %s to calculate difference.', names(comp)[1]))
		}else{
			by_predictor <- names(comp)[1]
			comp_levels <- comp[[1]][1:2]
		}
	}else{
		stop(sprintf('Grouping predictor %s not found in model.', names(comp)[1]))
	}

	if (is.null(ylab)) { 
		ylab = as.character(model$formula[2])
	}

	convert <- 1
	if(eegAxis==T) convert <- -1
			
	yvals <- sort(convert*ylim, decreasing=F)
	
	rngX = max(na.exclude(dat[,view])) - min(na.exclude(dat[,view]))
	
	np = 100

	if (is.null(comp_levels)) { 
		vals = sort(unique(dat[,by_predictor]))
	} else { 
		vals = comp_levels
	}

	grp1 <- dat[which(is.na(dat[,by_predictor])), ] 
	grp1[1:np,] = dat[1,]
	grp1[,view] = seq(min(na.exclude(dat[,view])),max(na.exclude(dat[,view])), by=rngX/(np-1))
	
	grp1[,by_predictor] = vals[1]

	#pred1 = predict.onlyInclude(model,grp1,onlyInclude=c(view,paste(by_predictor,vals[1],sep='')))
	pred1 = predict(model,grp1,type='lpmatrix')

	grp2 = grp1
	grp2[,by_predictor] = vals[2]

	#pred2 = predict.onlyInclude(model,grp2,onlyInclude=c(view,paste(by_predictor,vals[2],sep='')))
	pred2 = mgcv::predict.gam(model,grp2,type='lpmatrix')

	res1 = grp1
	res1$Xp <- pred2 - pred1

	res1$diff <- res1$Xp%*%coef(model)
	

	res1$XXX = res1[,view]
		
	if (is.null(main)) {
		mn = paste('Difference between',vals[1],'and',vals[2])
	} else { 
		mn = main
	}

	res1 = res1[order(res1$XXX),]
	se.diff <- sqrt(rowSums((res1$Xp%*%vcov(model))*res1$Xp))

	if(plot){	
		res1$ul <- res1$diff + 1.96 * se.diff
		res1$ll <- res1$diff - 1.96 * se.diff
		res1$ul99 <- res1$diff + 2.58 * se.diff
		res1$ll99 <- res1$diff - 2.58 * se.diff

		if (is.null(ylim)) yvals <- sort(c(convert*min(res1$ll),convert*max(res1$ul)),decreasing=F)
		plot(res1$XXX,convert*res1$diff,type='l',xlab=view, main=mn, ylab=ylab, ylim=yvals, axes=F,...)
		box()
		axis(1)
		axis(2, at=axTicks(2), labels=convert*axTicks(2))

		polygon(c(res1$XXX,rev(res1$XXX)),c(convert*res1$ul,rev(convert*res1$ll)),
				col=rgb(0.25,0.25,0.25,alpha=.25),border=NA) # shaded confidence interval

		abline(h=0)
		invisible(list(est=res1$diff, se.est=se.diff, xVals=res1$XXX))
	}else{
		return(list(est=res1$diff, se.est=se.diff, xVals=res1$XXX))
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
#' @param plotCI Logical: whether or not to plot confidence intervals.
#' @param color Colorpalette
#' @param nCol Range of colors of background of contour plot.
#' @param col Line color.
#' @param add.color.legend Logical: whether or not to add a color legend. 
#' Default is TRUE. If FALSE (omitted), one could use the function
#' \code{\link{gradientLegend}} to add a legend manually at any position.
#' @param n.grid Resolution.
#' @param nlevels Levels of contour lines.
#' @param zlim A two item array giving the lower and upper limits for the z-
#' axis scale. NULL to choose automatically.
#' @param main Title of plot.
#' @return If the result is not being plotted, a list is 
#' returned with the estimated difference (\code{est}) and the standard error 
#' over the estimate (\code{se.est}) and the x-values (\code{x}) is returned.
#' @author Martijn Wieling
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
plot_diff2 <- function(model,view, comp, plotCI=F, 
	color='topo', nCol=100, col=NULL, add.color.legend=TRUE,
	n.grid=30, nlevels=10, zlim=NULL, main=NULL) { 

	dat = model$model

	by_predictor <- NULL
	comp_levels <- NULL
	if(length(view) < 2){
		stop('Provide predictors for x- and y-axes in view.')
	}
	if(names(comp)[1] %in% colnames(dat)){
		if(length(comp[[1]]) < 2){
			stop(sprintf('Provide two levels for %s to calculate difference.', names(comp)[1]))
		}else{
			by_predictor <- names(comp)[1]
			comp_levels <- comp[[1]][1:2]
		}
	}else{
		stop(sprintf('Grouping predictor %s not found in model.', names(comp)[1]))
	}

	xvar <- view[1]
	yvar <- view[2]
	
	nX = n.grid
	nY = n.grid
	rngX = max(na.exclude(dat[,xvar])) - min(na.exclude(dat[,xvar]))
	rngY = max(na.exclude(dat[,yvar])) - min(na.exclude(dat[,yvar]))
	
	np = nX * nY

	if (is.null(comp_levels)) { 
		vals = sort(unique(dat[,by_predictor]))
	} else { 
		vals = comp_levels
	}

	# datasets maken, waarbij je alleen de variabele van interest varieert, om verschilcurves te kunnen maken
	grp1 <- dat[which(is.na(dat[,by_predictor])), ] # get same columns
	grp1[1:np,] = dat[1,]
	grp1[,xvar] = rep(seq(min(na.exclude(dat[,xvar])),max(na.exclude(dat[,xvar])), length=nX ),nY)
	grp1$nr = as.numeric(rownames(grp1))-1
	grp1[,yvar] = min(na.exclude(dat[,yvar])) + floor(grp1$nr / nX)*(rngY / (nY-1))
	
	grp1[,by_predictor] = vals[1]

	#pred1 = predict.onlyInclude(model,grp1,onlyInclude=c(xvar,yvar,paste(by_predictor,vals[1],sep='')))
	pred1 = mgcv::predict.gam(model,grp1,type='lpmatrix')

	grp2 = grp1
	grp2[,by_predictor] = vals[2]

	#pred2 = predict.onlyInclude(model,grp2,onlyInclude=c(xvar,yvar,paste(by_predictor,vals[2],sep='')))
	pred2 = mgcv::predict.gam(model,grp2,type='lpmatrix')

	res1 = grp1
	res1$Xp <- pred2 - pred1
	res1$diff <- res1$Xp%*%coef(model)


	res1$XXX = res1[,xvar]
	res1$YYY = res1[,yvar]
	z = matrix(res1$diff, nX, nY)
	m1 = seq(min(na.exclude(dat[,xvar])), max(na.exclude(dat[,xvar])), length=nX)
	m2 = seq(min(na.exclude(dat[,yvar])), max(na.exclude(dat[,yvar])), length=nY)
	
	if (is.null(main)) {
		mn = paste('Difference between',vals[1],'and',vals[2])
	} else { 
		mn = main
	}

   if (color == "heat") {
        pal <- heat.colors(nCol)
        con.col <- 3
    } else if (color == "topo") {
        pal <- topo.colors(nCol)
        con.col <- 2
    } else if (color == "cm") {
        pal <- cm.colors(nCol)
        con.col <- 1
    } else if (color == "bpy") {
        if (requireNamespace("sp", quietly = TRUE)) {
            pal <- sp::bpy.colors(nCol)
            con.col <- 1
        } else {
            warning("Package 'sp' needed for bpy color palette. Using topo.colors instead (default).")
            color <- 'topo'
            pal <- topo.colors(nCol)
            con.col <- 2
        }
    } else if (color == "terrain") {
        pal <- terrain.colors(nCol)
        con.col <- 2
    } else if (color == "gray" || color == "bw") {
        pal <- gray(seq(0.1, 0.9, length = nCol))
        con.col <- 1
    } else {
    	stop("color scheme not recognised")
    }

    if (!is.null(zlim)) {
        if (length(zlim) != 2 || zlim[1] >= zlim[2]) 
            stop("Something wrong with zlim")
    } else {
        zlim <- c(min(z, na.rm = TRUE), max(z, na.rm = TRUE))
    }

	image(m1,m2,z,col=pal,main=mn,xlab=xvar,ylab=yvar, zlim=zlim)
	contour(m1,m2,z,col=ifelse(is.null(col), con.col, col),nlevels=nlevels,add=T)

    if(add.color.legend){
        gradientLegend(round(c(min(zlim), max(zlim)), 3), n.seg=3, pos=.875, color=pal)
    }

	if (plotCI) { 
		se.diff <- sqrt(rowSums((res1$Xp%*%vcov(model))*res1$Xp))
		res1$ul <- res1$diff + 1.96 * se.diff
		res1$ll <- res1$diff - 1.96 * se.diff
		res1$ul99 <- res1$diff + 2.58 * se.diff
		res1$ll99 <- res1$diff - 2.58 * se.diff
		zu = matrix(res1$ul, nX, nY)
		contour(m1,m2,zu,col='red',nlevels=nlevels,add=T,lty=2)
		zl = matrix(res1$ll, nX, nY)
		contour(m1,m2,zl,col='green',nlevels=nlevels,add=T,lty=2)
	}
}


