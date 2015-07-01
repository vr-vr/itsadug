#' Visualization of smooths.
#'
#' @export
#' @description Plots a smooth from a \code{\link[mgcv]{gam}} or 
#' \code{\link[mgcv]{bam}} model based on predictions.
#' In contrast with the default \code{\link[mgcv]{plot.gam}}, this function 
#' plots the summed effects and optionally removes the random effects.
#'
#' @param x A gam object, produced by \code{\link[mgcv]{gam}} or 
#' \code{\link[mgcv]{bam}}.
#' @param view Text string containing the name of the smooth
#' to be displayed. Note that 
#' variables coerced to factors in the model formula won't work as view 
#' variables.
#' @param cond A named list of the values to use for the other predictor terms 
#' (not in view). Used for choosing between smooths that share the same view 
#' predictors.
#' @param plot_all A vector with a name / names of model predictors, 
#' for which all levels should be plotted.
#' @param rm.ranef Logical: whether or not to remove random effects. 
#' Default is TRUE.
#' @param n.grid  The number of grid nodes in each direction used for 
#' calculating the plotted surface. 
#' @param rug Logical: when TRUE (default) then the covariate to which the 
#' plot applies is displayed as a rug plot at the foot of each plot of a 1-d 
#' smooth. Setting to FALSE will speed up plotting for large datasets. 
#' @param col The colors for the lines and the error bars of the plot.
#' @param add Logical: whether or not to add the lines to an existing plot, or 
#' start a new plot (default).
#' @param se If less than or equal to zero then only the predicted surface is 
#' plotted, but if greater than zero, then the predicted values plus 
#' confidence intervals are plotted. The value of se will be multiplied with 
#' the standard error (i.e., 1.96 results in 95\%CI and 2.58).
#' @param shade Logical: Set to TRUE to produce shaded regions as confidence 
#' bands for smooths (not avaliable for parametric terms, which are plotted 
#' using termplot).
#' @param eegAxis Logical: whether or not to reverse the y-axis, plotting the 
#' negative amplitudes upwards as traditionally is done in EEG research.
#' If eeg.axes is TRUE, labels for x- and y-axis are provided, when not 
#' provided by the user. Default value is FALSE.
#' @param print.summary Logical: whether or not to print summary.
#' Default set to the print info messages option 
#' (see \code{\link{infoMessages}}).
#' @param main Changing the main title for the plot, see also title.
#' @param xlab Changing the label for the x axis, 
#' defaults to a description of x.
#' @param ylab Changing the label for the y axis, 
#' defaults to a description of y.
#' @param ylim the y limits of the plot.
#' @param h0 A vector indicating where to add solid horizontal lines for 
#' reference. By default no values provided.
#' @param v0 A vector indicating where to add dotted vertical lines for 
#' reference. By default no values provided.
#' @param transform Function for transforming the fitted values. 
#' Default is NULL.
#' @param legend_plot_all Legend location. This could be a keyword from 
#' the list "bottomright", "bottom", "bottomleft", "left", "topleft", "top", 
#' "topright", "right" and "center", or a list with \code{x} and \code{y} 
#' coordinate (e.g., \code{list(x=0,y=0)}). 
#' @param hide.label Logical: whether or not to hide the label 
#' (i.e., "fitted values"). Default is FALSE.
#' @param ... other options to pass on to lines and plot, 
#' see \code{\link[graphics]{par}}
#' @section Notes:
#' This function plots the summed effects, including intercept and other 
#' predictors. For plotting partial effects, see the function 
#' \code{\link[mgcv]{plot.gam}}, or see the examples with 
#' \code{\link{get_modelterm}} for more flexibility (e.g., plotting using the 
#' \code{lattice} package or \code{ggplots}).
#'
#' @examples
#' data(simdat)
#' 
#' \dontrun{
#' # Model with random effect and interactions:
#' m1 <- bam(Y ~ te(Time, Trial)+s(Time, Subject, bs='fs', m=1),
#'     data=simdat)
#'
#' # Default plot produces only surface of Time x Trial:
#' plot(m1, select=1)
#' # Only the Time component:
#' plot_smooth(m1, view="Time")
#' # Note the summary that is printed.
#'
#' # without random effects:
#' plot_smooth(m1, view="Time", rm.ranef=TRUE)
#' 
#' # Plot summed effects:
#' dev.new(width=8, height=4) # use x11(,8,4) on Linux
#' par(mfrow=c(1,2))
#' fvisgam(m1, view=c("Time", "Trial"), 
#'     plot.type='contour', color='topo', main='interaction',
#'     rm.ranef=TRUE)
#' arrows(x0=0, x1=2200, y0=-5, y1=-5, col='red', 
#'     code=2, length=.1, lwd=2, xpd=TRUE)
#' plot_smooth(m1, view='Time', cond=list(Trial=-5),
#'     main='Trial=-5', rm.ranef=TRUE)
#'
#'
#' # Model with random effect and interactions:
#' m2 <- bam(Y ~ Group + s(Time, by=Group)
#'     +s(Time, Subject, bs='fs', m=1),
#'     data=simdat)
#' 
#' # Plot all levels of a predictor:
#' plot_smooth(m2, view='Time', plot_all="Group",
#'     rm.ranef=TRUE)
#' # It also possible to combine predictors in plot_all.
#' # Note: this is not a meaningfull plot, 
#' # just for illustration purposes!
#' plot_smooth(m2, view='Time', plot_all=c("Group", "Subject"))
#' }
#'
#' # and for a quick overview of plotfunctions:
#' vignette("overview", package="itsadug")
#'
#' @author Jacolien van Rij and Martijn Wieling. 
#' @seealso \code{\link[mgcv]{plot.gam}}, \code{\link{plot_diff}} 
#'
#' @family functions for interpreting nonlinear effects

plot_smooth <- function(x, view = NULL, cond = list(), 
    plot_all=NULL, rm.ranef=NULL,
    n.grid = 30, rug = TRUE, col = 'black', add=FALSE, 
    se = 1.96, shade = TRUE, eegAxis=FALSE, 
    print.summary=getOption('itsadug_print'),
    main=NULL, xlab=NULL, ylab=NULL, ylim=NULL, h0=0, v0=NULL, 
    transform=NULL, legend_plot_all=NULL, 
    hide.label=FALSE, ...) {
       
    dnm <- names(list(...))

    v.names <- names(x$var.summary)

    if (is.null(view)) {
        stop("Specify one view predictors for the x-axis.")
    } else {
        if(length(view)>1){
            warning('Only first element of view will be used.')
            view <- view[1]
        }
        if (sum(view %in% v.names) != 1) {
            stop(paste(c("View variable must be one of", v.names), collapse = ", "))
        }
        if (!inherits(x$var.summary[[view[1]]], c("numeric"))){
            stop("Don't know what to do with parametric terms that are not simple numeric variables.")
        }
    }

    if(!is.null(cond)){
        cn <- names(cond)
        test <- sapply(cn, function(x){
            if(length(unique(cond[[x]]))>1){
                stop("Do not specify more than 1 value for conditions listed in the argument cond.")
            }else{
                TRUE
            }
        })
    }
    if(!is.null(plot_all)){
        if(!is.vector(plot_all)){
            stop("Argument plot_all should be a vector with predictor names.")
        }else{
            # check if plot_all in cond
            if(any(plot_all %in% names(cond))){
                warning(sprintf("%s in cond and in plot_all. plot_all is being ignored.",
                    paste(plot_all[plot_all %in% names(cond)], collapse=', ')))
                plot_all <- plot_all[!plot_all %in% names(cond)]
            }
            # check if plot_all are column names
            if(any(! plot_all %in% v.names)){
                warning(sprintf("%s (specified in plot_all) is not a model predictor. Will be ignored.",
                    paste(plot_all[!plot_all %in% v.names], collapse=', ')))
                plot_all <- plot_all[plot_all %in% v.names]
            }
            # check length:
            if(length(plot_all)>0){
                for(i in plot_all){
                    cond[[i]] <- unique(as.character(x$model[,i]))
                }
            }else{
                plot_all <- NULL
            }           
        }

    }
    

    m1 <- seq(min(x$var.summary[[view[1]]], na.rm=TRUE), 
        max(x$var.summary[[view[1]]], na.rm=TRUE), length=n.grid)

    cond[[view[1]]] <- m1

    newd <- get_predictions(x, cond=cond, se=ifelse(se>0, TRUE, FALSE), 
        f=ifelse(se>0, se, 1.96), rm.ranef=rm.ranef,
        print.summary=print.summary)

    if(se > 0){
        newd$ul <- with(newd, fit+CI)
        newd$ll <- with(newd, fit-CI)

        if(!is.null(transform)){
            newd$ul <- sapply(newd$ul, transform)
            newd$ll <- sapply(newd$ll, transform)
        }
    }
    if(!is.null(transform)){
        newd$fit <- sapply(newd$fit, transform)
    }


    if(is.null(main)){ main <- view[1] }
    if(is.null(xlab)){ xlab <- view[1] }
    if(is.null(ylab)){ ylab <- names(x$model)[!names(x$model) %in% v.names]}
    if(is.null(ylim)){ 
        if(se>0){
            ylim <- range(c(newd$ul, newd$ll))
        } else {
            ylim <- range(newd$fit)
        }
    }
        

    if(add==FALSE){
        emptyPlot(range(newd[,view[1]]), ylim,
            main=main, xlab=xlab, ylab=ylab,
            h0=h0, v0=v0, eegAxis=eegAxis, ...)
        if(hide.label==FALSE){
            addlabel = "fitted values"
            if(!is.null(rm.ranef)){
                if(rm.ranef !=FALSE){
                    addlabel = paste(addlabel, "excl. random", sep=", ")
                }
            }
            mtext(addlabel, side=4, line=0, adj=0, 
                cex=.75, col='gray35', xpd=TRUE)

            if(!is.null(transform)){
                mtext("transformed", side=4, line=.75, adj=0, 
                cex=.75, col='gray35', xpd=TRUE)
            }
        }        
    }

    if(rug==TRUE){
        rug(x$model[,view[1]])
    }

    if(!is.null(plot_all)){
        alllevels <- c()
        plotlevels <- c()
        if(length(plot_all)>1){
            tmpname <- sub("/", "", tempfile(pattern = "event", 
                tmpdir = "plotsmooth", fileext = ""), fixed=TRUE)
            newd[,tmpname] <- interaction(newd[, plot_all])

            alllevels <- length(levels(newd[,tmpname]))
            plotlevels <- levels(newd[,tmpname])
            cnt <- 1
            for(i in levels(newd[,tmpname])){
                if(se > 0){
                    plot_error(newd[newd[,tmpname]==i,view[1]], 
                        newd[newd[,tmpname]==i,]$fit, 
                        newd[newd[,tmpname]==i,]$ul, 
                        se.fit2=newd[newd[,tmpname]==i,]$ll, 
                        shade=shade, f=1, col=rainbow(alllevels)[cnt], ...)
                }else{
                    lines(newd[newd[,tmpname]==i,view[1]], 
                        newd[newd[,tmpname]==i,]$fit, 
                        col=rainbow(alllevels)[cnt], ...)
                }
                cnt <- cnt+1
            }
            newd[, tmpname] <- NULL
        }else{
            alllevels <- length(levels(newd[,plot_all]))
            plotlevels <- levels(newd[,plot_all])
            cnt <- 1
            for(i in levels(newd[,plot_all])){
                if(se > 0){
                    plot_error(newd[newd[,plot_all]==i,view[1]], 
                        newd[newd[,plot_all]==i,]$fit, 
                        newd[newd[,plot_all]==i,]$ul, 
                        se.fit2=newd[newd[,plot_all]==i,]$ll, 
                        shade=shade, f=1, col=rainbow(alllevels)[cnt], ...)
                }else{
                    lines(newd[newd[,plot_all]==i,view[1]], 
                        newd[newd[,plot_all]==i,]$fit, col=rainbow(alllevels)[cnt], ...)
                }
                cnt <- cnt+1
            }       
        }
        # add legend:
        if(is.null(legend_plot_all)){
            gfc <- getFigCoords()
            legend(gfc[2], gfc[4],
                legend=plotlevels,
                text.col=rainbow(alllevels),
                text.font=2,
                xjust=1, yjust=1,
                bty='n', xpd=TRUE)
        }else{
            legend(legend_plot_all,
                legend=plotlevels,
                text.col=rainbow(alllevels),
                text.font=2,
                bty='n', xpd=TRUE)
        }
    }else{
        if(se > 0){
            plot_error(newd[,view[1]], newd$fit, newd$ul, se.fit2=newd$ll, shade=shade, f=1, col=col, ...)
        }else{
            lines(newd[,view[1]], newd$fit, col=col, ...)
        }
    }
    
    invisible(list(fv = newd, rm.ranef=rm.ranef, transform=transform))
}

 
