#' Visualization of the model fit for time series data.
#'
#' @export
#' @import mgcv
#' @import stats
#' @import grDevices
#' @import graphics
#' @description Plots the fitted values and the data for \code{n} 
#' trials of time series data. For example, plots \code{n} trials 
#' of the same participant.
#'
#' @param x A lm or gam object, produced by \code{\link[mgcv]{gam}} or 
#' \code{\link[mgcv]{bam}}, \code{\link[stats]{lm}}, \code{\link[stats]{glm}}.
#' @param data Data frame on which the model was fitted.
#' @param view Text string containing the predictor or column in the data 
#' to be displayed on the x-axis. 
#' Note that variables coerced to factors in the model formula 
#' won't work as view variables.  
#' @param event column name from the data 
#' that specifies the time series from which \code{n} are being plotted.
#' @param n Number of time series to plot. Default is 3. Set to -1 for plotting 
#' all time series (which may take a considerable time).
#' @param random Numeric: if set to TRUE (default), \code{n} random events are 
#' selected to plot. If set to FALSE, the first \code{n} events are selected 
#' to plot. The events could be precisely controlled with the argument 
#' \code{cond}.
#' @param cond A named list of the values to use for the other predictor terms 
#' (not in view) or to select specific trials or time series to plot. 
#' @param col Two value vector specifiying the colors for the data and 
#' the modelfit respectively.
#' @param add Logical: whether or not to add the lines to an existing plot, or 
#' start a new plot (default).
#' @param se If less than or equal to zero then only the predicted surface is 
#' plotted, but if greater than zero, then the predicted values plus 
#' confidence intervals are plotted. The value of se will be multiplied with 
#' the standard error (i.e., 1.96 results in 95\%CI and 2.58).
#' @param eegAxis Logical: whether or not to reverse the y-axis, plotting the 
#' negative amplitudes upwards as traditionally is done in EEG research.
#' If eeg.axes is TRUE, labels for x- and y-axis are provided, when not 
#' provided by the user. Default value is FALSE.
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
#' @param hide.label Logical: whether or not to hide the label 
#' (i.e., "fitted values"). Default is FALSE.
#' @param hide.legend Logical: whether or not to hide the legend. 
#' Default is FALSE.
#' @param ... other options to pass on to lines and plot, 
#' see \code{\link[graphics]{par}}
#' @section Notes:
#' This function plots the fitted effects, including intercept and other 
#' predictors. 
#'
#' @examples
#' data(simdat)
#' 
#' # Create grouping predictor for time series:
#' simdat$Event <- interaction(simdat$Subject, simdat$Trial)
#' 
#' # model without random effects:
#' m1 <- bam(Y ~ te(Time, Trial),
#'     data=simdat)
#' plot_modelfit(m1, view="Time", event=simdat$Event)
#'
#' # colorizing residuals:
#' plot_modelfit(m1, view="Time", event=simdat$Event, fill=TRUE)
#' 
#' # All trials of one subject:
#' \dontrun{
#' # this produces error:
#' plot_modelfit(m1, view="Time", event=simdat$Event, 
#'     cond=list(Subject="a01"), n=-1)
#' }
#' # instead try this:
#' simdat$Subj <- ifelse(simdat$Subject=="a01", TRUE, FALSE)
#' plot_modelfit(m1, view="Time", event=simdat$Event, 
#'     cond=list(Subject=simdat$Subj), n=-1)
#' 
#' \dontrun{
#' # Model with random intercepts for subjects:
#' m2 <- bam(Y ~ te(Time, Trial)+s(Subject, bs='re'),
#'     data=simdat)
#' # now selecting a subject works, because it is in the model:
#' plot_modelfit(m2, view="Time", event=simdat$Event, 
#'     cond=list(Subject="a01"), n=-1, ylim=c(-13,13))
#'
#' # Model with random effect and interactions:
#' m3 <- bam(Y ~ te(Time, Trial)+s(Time, Subject, bs='fs', m=1),
#'     data=simdat)
#' plot_modelfit(m3, view="Time", event=simdat$Event, 
#'     cond=list(Subject="a01"), n=-1, ylim=c(-13,13))
#' }
#' @author Jacolien van Rij
#'

plot_modelfit <- function(x, view, event=NULL, 
	n=3, random=TRUE, cond = NULL, 
   	col = c(alpha(1), 'red'), add=FALSE, eegAxis=FALSE, 
   	fill=FALSE, rug=TRUE, 
    main=NULL, xlab=NULL, ylab=NULL, ylim=NULL, h0=0, v0=NULL, 
    transform=NULL, 
    hide.label=FALSE, hide.legend=FALSE, ...) {
       
    dnm <- names(list(...))

    v.names <- names(x$var.summary)
    
    if(is.null(main)){ main <- "Model fit" }
    if(is.null(xlab)){ xlab <- view }
    if(is.null(ylab)){ ylab <- as.character(x$formula[2])}

    dat <- x$model
    y <- as.character(x$formula[2])
    viewcol <- sprintf("%s%d", "tmp", sample.int(1e+06, size = 1L))
    eventcol <- sprintf("%s%d", "ev", sample.int(1e+06, size = 1L))
    plot.events <- NULL

    missing <- missing_est(x)

    if (is.null(view)) {
        stop("Specify one view predictor for the x-axis, either the name of a model predictor or a vector.")
    } else {
        if(length(view)>1){
        	if(view[1] %in% v.names){
        		view=view[1]
        		warning("Only first element of view is being used.")
        	}else{
        		stop("View is not column name of data.")
        	}
        }else{
        	if (sum(view %in% v.names) != 1) {
            	stop(paste(c("View variable must be one of", v.names), collapse = ", "))
        	}
	        if (!inherits(x$var.summary[[view]], c("numeric"))){
	            stop("Don't know what to do with parametric terms that are not simple numeric variables.")
	        }
        }
        dat[,viewcol] <- dat[,view]
    }
    if (!is.null(event)) {
        if(length(event)>1){
        	if(length(event)==nrow(dat)){
        		dat[, eventcol] <- event 
        	}else if(length(event)==(nrow(dat)+length(missing))){
        		dat[, eventcol] <- event[-missing] 
        	}else if(event[1] %in% v.names){
        		event=event[1]
        		dat[, eventcol] <- dat[,event]
        		warning("Only first element of event is being used.")
        	}else{
        		stop("Event is not column name of data.")
        	}
        }else{
        	if (sum(event %in% v.names) != 1) {
            	stop(paste(c("Event variable must be one of", v.names), collapse = ", "))
        	}else{
        		dat[, eventcol] <- dat[,event]
        	}
        }
    }
    dat$fit <- fitted(x)

    if(!is.null(cond)){
        cn <- names(cond)
        for(icn in cn){
        	if(icn %in% v.names){
        		dat <- dat[dat[,icn] %in% cond[[icn]],]
        	}else if(length(cond[[icn]])==nrow(dat)){
        		if(is.logical(cond[[icn]])){
        			dat <- dat[cond[[icn]]==TRUE,]
        		}else{
        			stop(sprintf("%s not found in the model. Provide a column of TRUE (include) and FALSE (do not include) of the size of the data.", icn))
        		}        		
        	}else if(length(cond[[icn]])==(nrow(dat)+length(missing))){
        		if(is.logical(cond[[icn]])){
        			dat <- dat[cond[[icn]][-missing]==TRUE,]
        		}else{
        			stop(sprintf("%s not found in the model. Provide a column of TRUE (include) and FALSE (do not include) of the size of the data.", icn))
        		}    
        	}else{
        		stop(sprintf("%s not found in the model. Provide a column of TRUE (include) and FALSE (do not include) of the size of the data.", icn))
        	}
        }
    }


    # sample events:
    dat <- droplevels(dat)

    if(!is.null(event)){
	    events <- unique(dat[,eventcol], na.rm=TRUE)
	    
	    if(n < 0){
	    	n <- length(events)
	    }

	    plot.events <- events[1:n]
	    if(random){    
	    	plot.events <- sample(events, min(length(events), n))
	    }
		if(is.factor(plot.events)){
	    	plot.events <- as.character(plot.events)
	    }

	    # select events:
	    dat <- droplevels( dat[dat[,eventcol] %in% plot.events,] )
	}

    if(!is.null(transform)){
        dat$fit <- sapply(dat$fit, transform)
    }
    if(is.null(ylim)){ 
        ylim <- range(c(dat$fit, dat[,y]), na.rm=TRUE)
    }

    if(add==FALSE){
        emptyPlot(range(dat[,view]), ylim,
            main=main, xlab=xlab, ylab=ylab,
            h0=h0, v0=v0, eegAxis=eegAxis, ...)
        if(hide.label==FALSE){
            addlabel = "fitted values"
            
            mtext(addlabel, side=4, line=0, adj=0, 
                cex=.75, col='gray35', xpd=TRUE)

            if(!is.null(transform)){
                mtext("transformed", side=4, line=.75, adj=0, 
                cex=.75, col='gray35', xpd=TRUE)
            }
        }      
    }
    if(names(dev.cur())[1] %in% c("X11", "postscript", "xfig", "pictex") ){
    	if(length(col[grepl("\\#.+", col)]) > 0 ){
    		hexcol <- which( grepl("\\#.+", col) )
    		col[hexcol] <- substr(col[hexcol], 1,7)
    	}
    	if(fill){
    		warning(sprintf("%s device does not allow for transparent colors. Fill color will be set to color 1.", names(dev.cur())[1]))
    	}
    }  

    if(rug==TRUE){
        rug(x$model[,view])
    }

    if(!is.null(event)){
	    for(i in plot.events){
	    	newd <- NULL
	    	newd <- droplevels(dat[dat[, eventcol]==i,])
	    	if(nrow(newd)>1){
	    		# plot data:
	    		if(fill){
	    			fill_area(newd[,viewcol], newd$fit, from=newd[,y], 
	    				col=col[1], border=col[2], outline=FALSE)
	    			lines(newd[,viewcol], newd[,y], col=col[1], ...)
	    		}else{
	    			lines(newd[,viewcol], newd[,y], col=col[1], ...)
	    			lines(newd[,viewcol], newd$fit, col=col[2], ...)
	    		}
	    	} else if(nrow(newd)==1){
	    		# plot data:
	    		points(newd[,viewcol], newd[,y], col=col[1], ...)
	    		points(newd[,viewcol], newd$fit, col=col[2], ...)
	    	} else{
	    		if(print.summary){
	    			message(sprintf("No data for event %s. Ignored.", i))
	    		}
	    	} 
	    }
	}else{
		newd <- droplevels(dat)
    	if(nrow(newd)>1){
    		# plot data:
    		lines(newd[,viewcol], newd[,y], col=col[1], ...)
    		lines(newd[,viewcol], newd$fit, col=col[2], ...)
    	} else if(nrow(newd)==1){
    		# plot data:
    		points(newd[,viewcol], newd[,y], col=col[1], ...)
    		points(newd[,viewcol], newd$fit, col=col[2], ...)
    	} else{
    		if(print.summary){
	    		message("No data to be plotted.")
	    	}
    	} 	
	}

   	# add legend:
   	if(!hide.legend){
	   	gfc <- getFigCoords('f')
	   	legend(gfc[2], gfc[4],
	   		xjust=1, yjust=1,
	   		col=col, lwd=2, seg.len=1.5,
	   		legend=c("data", "model fit"),
	   		bty='n', cex=.85,  xpd=TRUE)
 	}
    #output

    invisible( list(events = ifelse(!is.null(event),plot.events, NA), subdat = dat[,!colnames(dat) %in% c(viewcol, eventcol)], 
    	view = dat[,viewcol], event = ifelse(!is.null(event),dat[, eventcol], NA)  ) )
}

 
