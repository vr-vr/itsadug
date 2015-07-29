#' Inspection of random factor smooths.
#' 
#' @export
#' @import mgcv
#' @import stats
#' @import grDevices
#' @import graphics
#' @param model A gam object, produced by \code{\link[mgcv]{gam}} or 
#' \code{\link[mgcv]{bam}}.
#' @param select A number, indicating the model term to be selected. 
#' @param fun A string or function description to apply to the random effects 
#' estimates. When NULL (default), the estimates for the random effects are 
#' returned. 
#' @param cond A named list of the values to restrict the estimates for the 
#' random predictor terms. When NULL (default) all levels are returned.
#' @param n.grid Number of data points estimated for each random smooth.
#' @param print.summary Logical: whether or not to print a summary of the 
#' values selected for each predictor. 
#' Default set to the print info messages option 
#' (see \code{\link{infoMessages}}).
#' @param plot Logical: whether or not to plot the random effect estimates 
#' (TRUE by default).
#' @param add Logical: whether or not to add the random effect estimates 
#' to an existing plot (FALSE by default).
#' @param main Changing the main title for the plot, see also title.
#' @param xlab Changing the label for the x axis, 
#' defaults to a description of x.
#' @param ylab Changing the label for the y axis, 
#' defaults to a description of y.
#' @param ylim Changing the y limits of the plot.
#' @param h0 A vector indicating where to add solid horizontal lines 
#' for reference. By default 0.
#' @param v0 A vector indicating where to add dotted vertical lines 
#' for reference. By default no values provided.
#' @param col Specifying the colors of the lines.
#' @param eegAxis Whether or not to reverse the y-axis 
#' (plotting negative upwards).
#' @param ... other options to pass on to \code{\link[graphics]{lines}}, 
#' see \code{\link[graphics]{par}}
#' @return A data frame with estimates for random effects
#' is optionally returned.
#' @examples
#' # load data:
#' data(simdat)
#'
#' \dontrun{
#' # Condition as factor, to have a random intercept
#' # for illustration purposes:
#' simdat$Condition <- as.factor(simdat$Condition)
#'
#' # Model with random effect and interactions:
#' m2 <- bam(Y ~ s(Time) + s(Trial)
#' + ti(Time, Trial)
#' + s(Condition, bs='re')
#' + s(Time, Subject, bs='fs', m=1),
#' data=simdat)
#'
#' # extract with wrong select value:
#' newd <- inspect_random(m2, select=4)
#' # results in warning, automatically takes select=5
#' head(newd)
#' inspect_random(m2, select=5, cond=list(Subject=c('a01','a02','a03')))
#' 
#' # Alternatively, fix random effect of Condition, and plot 
#' # random effects for subjects with lattice:
#' newd <- inspect_random(m2, select=5,
#'     cond=list(Subject=unique(simdat[simdat$Condition==0,'Subject'])),
#'     plot=FALSE)
#'
#' # Make lattice plot:
#' require(lattice)
#' lattice::xyplot(fit~Time | Subject,
#'     data=newd, type="l",
#'     xlab="Time", ylab="Partial effect")
#'
#' # Using argument 'fun':
#' inspect_random(m2, select=5, fun=mean, 
#'     cond=list(Subject=unique(simdat[simdat$Condition==0,'Subject'])))
#' inspect_random(m2, select=5, fun=mean, 
#'     cond=list(Subject=unique(simdat[simdat$Condition==2,'Subject'])),
#'     col='red', add=TRUE)
#' }
#'
#' # see the vignette for examples:
#' vignette("overview", package="itsadug")
#' @author Jacolien van Rij
#' @family functions for model predictions

inspect_random <- function(model, select=1, fun=NULL, 
	cond=NULL, n.grid=30, print.summary=getOption('itsadug_print'),
	plot=TRUE, add=FALSE,
	main=NULL, xlab=NULL, ylab=NULL, ylim=NULL, h0=0, v0=NULL, 
	col=NULL, eegAxis=FALSE, ...
	){

	if(!"lm" %in% class(model)){
		stop("This function does not work for class %s models.", class(model)[1])
	}else{

		# find random effects:
		smoothlabels <- as.data.frame( do.call('rbind', 
			lapply(model$smooth, 
				function(x){
					data.frame(Label=x[['label']], 
						Dim=x[['null.space.dim']], 
						Class = attr(x, "class")[1],
						stringsAsFactors=FALSE)
				} ) ) )
		# smoothlabels <- smoothlabels[smoothlabels$Dim==0,c("Label", "Class")]
		fslabels <- as.vector( smoothlabels[smoothlabels$Class %in% c("fs.interaction"), "Label"] )
		if(length(fslabels) == 0){
			warning("No random smooths / factor smooths found in the model.")
			return(NULL)
		}else if( !model$smooth[[select]]$label %in% fslabels){
			select.new <- which(smoothlabels$Label==fslabels[1])
			warning(sprintf("Selected modelterm %d ( '%s' ) is not a factor smooth. Modelterm %d selected instead ( '%s' ).",
				select, model$smooth[[select]]$label,
				select.new,  model$smooth[[select.new]]$label))
			select = select.new
		}

		if(length(model$smooth[[select]]$term) > 2){
			stop("This function is implemented for factor smooths with at most two terms, e.g. s(A,B, bs='fs').")
		}
	
		numterm   <- c()
		groupterm <- c()
	
		if(!is.null(fun)){
			fv <- get_modelterm(model, select=select, cond=cond, se=FALSE,
				print.summary=print.summary, as.data.frame=TRUE)
						
			fun.cond <- list()
			fun.val = list()
			
			for(j in model$smooth[[select]][['term']]){
				if(!inherits(model$model[,j],"factor")){
					fun.cond[[j]] <- fv[,j]
					numterm <- c(numterm, j)
				}else{
					groupterm <- c(groupterm, j)
				}
			}		
		    if(is.null(main)){ main <- model$smooth[[select]]$label }
		    if(is.null(xlab)){ xlab <- numterm[1] }
		    if(is.null(ylab)){ ylab <- sprintf("est. %s", names(model$model)[1])}	

			if(length(fun.cond) > 0){
				fun.val[['values']] <- aggregate(list(x=fv$fit), by=fun.cond, fun)

				if(plot==TRUE){
					if(add==FALSE){
						if(is.null(ylim)){ 
	           				ylim <- range(fv$fit)
	    				}
				    if(is.null(col)){
				    	col=1
				    }
			        emptyPlot(range(fun.val[['values']][,numterm[1]]), ylim,
			            main=main, xlab=xlab, ylab=ylab,
			            h0=h0, v0=v0, eegAxis=eegAxis)				    
					}
					lines(fun.val$values[,numterm[1]], fun.val$values[,'x'], col=col, ...)
				}

			} else {
				fun.val[['values']] <- unlist(lapply(list(fv$fit), fun))
				if(plot==TRUE){
					warning('Plotting for 1 dimensional factor smooth not implemented yet.')
				}
			}

			for(nc in names(cond)){
				fun.val[[nc]] <- cond[[nc]]
			}
			
			fun.val[['function']] <- fun
			invisible(fun.val)
		}else{
			fv <- get_modelterm(model, select=select, cond=cond,
				print.summary=print.summary, as.data.frame=TRUE)

			for(j in model$smooth[[select]][['term']]){
				if(!inherits(model$model[,j],"factor")){
					numterm <- c(numterm, j)
				}else{
					groupterm <- c(groupterm, j)
				}
			}	

			if(plot==TRUE){
				if(add==FALSE){
					if(is.null(ylim)){ 
           				ylim <- range(fv$fit)
    				}
			    
		        	emptyPlot( range(fv[,numterm[1]]), ylim,
			            main=main, xlab=xlab, ylab=ylab,
			            h0=h0, v0=v0, eegAxis=eegAxis, ...)				    
				}
				if(is.null(col)){
					count <- 1
					for(su in unique(fv[,groupterm[1]])){
						tmp <- fv[fv[,groupterm[1]]==su,]
						lines(tmp[,numterm[1]], tmp[,'fit'], col=count, lty=count, ...)
						count <- count+1
					}
				}else{
					count <- 1
					for(su in levels(fv[,groupterm[1]])){
						tmp <- fv[fv[,groupterm[1]]==su,]
						lines(tmp[,numterm[1]], tmp[,'fit'], col=col, lty=count, ...)
						count <- count+1
					}
				}
			}

			invisible(fv)
		}
	}
}

