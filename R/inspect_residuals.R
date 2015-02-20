#' Compare distribution of data with normal distribution.
#' 
#' @export
#' @param res Vector with residuals or other data for which the distribution .
#' @param col Color for filling the area. Default is black.
#' @param col.normal Color for shading and line of normal distribution.
#' @param legend.pos Position of legend, can be string (e.g., 'topleft') or an 
#' \code{\link[grDevices]{xy.coords}} object.
#' @param legend.label Text string, label for plotted data distribution.
#' @param ... Optional arguments for the lines. See \code{\link{par}}.
#' @section Note:
#' Assumes centered data as input.
#' @examples
#' set.seed(123)
#' # normal distribution:
#' test <- rnorm(1000)
#' check_normaldist(test)
#' # t-distribution:
#' test <- rt(1000, df=5)
#' check_normaldist(test)
#' # skewed data, e.g., reaction times:
#' test <- exp(rnorm(1000, mean=.500, sd=.25))
#' check_normaldist(test)
#' # center first:
#' check_normaldist(scale(test))
#' @family functions for model criticism
#' @author Jacolien van Rij

check_normaldist <- function(res, col='red', col.normal='black', 
	legend.pos='topright', legend.label='data', ...){
    x <- sort(res[!is.na(res)])
    sd.x <- sd(x)
    d <- density(x)
    parlist <- list(...)

    emptyPlot(range(d$x), range(d$y),
        main='Density', xlab=deparse(substitute(res)))
    fill_area(x, dnorm(x, sd=sd.x), col=col.normal)
    lines(x, dnorm(x, sd=sd.x), col=col.normal)
    if('lwd' %in% names(parlist)){
        lwd <- NULL
        lines(d$x, d$y, col=col, ...)
    }else{
        lines(d$x, d$y, col=col, lwd=2, ...)
    }
    
    if(!is.null(legend.pos)){
        legend(legend.pos, 
            legend=legend.label,
            col=c(col, col.normal), seg.len=1,
            lwd=c(ifelse('lwd' %in% names(parlist), parlist[['lwd']], 2), 1),
            bty='n')
    }
}



#' Inspect residuals of regression models.
#' 
#' @export
#' @param model A regression model, resulting from the functions
#' \code{\link[mgcv]{gam}} or \code{\link[mgcv]{bam}}, or \code{lm},
#' \code{glm}, \code{lmer}, or \code{glmer}.
#' @param  AR_start Defaults to NULL. 
#' Only use this when the model was run in an old versions of package 
#' \code{mgcv} and the function cannot retrieve the used AR.start values from 
#' the \code{model}. When an error is shown with newer versions of 
#' \code{mgcv}, please check the column provided as values of AR.start.
#' when using old versions of package \code{mgcv}.
#' Function will give error when it cannot find AR.start. 
#' @param split_by A names list indicating time series in the data.
#' @section Note:
#' In the ACF plots, black lines indicate ACF of standard residuals, red 
#' vertical lines indicate the ACF for the residuals corrected for the AR1 
#' model included.
#' @examples
#' data(simdat)
#'
#' \dontrun{
#' # Add start event column:
#' simdat <- start_event(simdat, event=c("Subject", "Trial"))
#' head(simdat)
#' # bam model with AR1 model (toy example, not serious model):
#' m1 <- bam(Y ~ Group + te(Time, Trial, by=Group), 
#'    data=simdat, rho=.5, AR.start=simdat$start.event)
#' # No time series specified:
#' check_resid(m1)
#' # Time series specified:
#' check_resid(m1, split_by=list(Subject=simdat$Subject, Trial=simdat$Trial))
#' # Note: residuals do not look very good.
#'
#' # This does not work (see acf_resid for similar examples):
#' check_resid(m1, split_by=c("Subject", "Trial"))
#' # However, it does work for included predictors:
#' check_resid(m1, split_by=c("Group"))
#' }
#' @author Jacolien van Rij

check_resid <- function(model, AR_start = NULL, split_by=NULL){
	par(mfrow=c(2,2), cex=1.1)

	el.narm <- NULL
	res <- resid(model)
	res.rho <- NULL

	# retrieve res.rho
	if("lm" %in% class(model)){
		el.narm <- missing_est(model)

		if(!is.null(AR_start)){
			if(length(AR_start) > length(res)){
				AR_start <- AR_start[-el.narm]
			}
		}
		if(!is.null(model$AR1.rho)){
			if(model$AR1.rho > 0){
				suppressWarnings( res.rho <- resid_gam(model, incl_na=TRUE) )
			}
		}
	}else if("lmerMod" %in% class(model)){
		res <- resid(model)
	}

	# Retrieve split predictors:

	if(!is.null(split_by)){
		if(is.list(split_by)){
			if(length(split_by[[1]]) > length(res)){
				if(is.null(el.narm)){
					stop("More values in split_by data than in the model residuals. Please exclude missing values from the split_by data.")
				}else{
					if(length(el.narm) > 0){
						suppressWarnings( split_by <- lapply(split_by, function(x){x[-el.narm]}) )
					}			
				}	
			}

			if(length(split_by[[1]]) < length(res)){
				if(is.null(el.narm)){
					res <- resid(model, na.action=na.exclude)
				}else{
					res <- res[-el.narm]
				}
			}
			if(length(split_by[[1]]) < length(res.rho)){
				if(is.null(el.narm)){
					res.rho <- NULL
				}else{
					res.rho <- res.rho[-el.narm]
				}
			}

		}else{
			dat <- c()
			if("lm" %in% class(model)){
				dat <- model$model
			}else if("lmerMod" %in% class(model)){
				dat <- model@frame
			}
			if(!all(split_by %in% colnames(dat))){
				notindata <- paste(split_by[!split_by %in% colnames(dat)], collapse=", ")
				stop(sprintf("split_by value(s) %s is / are not included as predictor in the model.", 
					notindata))
			}else{
				split_list <- list()
				for(i in split_by){
					split_list[[i]] <- as.vector(dat[,i])
				}
				split_by <- split_list

				if(length(split_by[[1]]) < length(res)){
					if(is.null(el.narm)){
						res <- resid(model, na.action=na.exclude)
					}else{
						res <- res[-el.narm]
					}
				}
				if(length(split_by[[1]]) < length(res.rho)){
					if(is.null(el.narm)){
						res.rho <- NULL
					}else{
						res.rho <- res.rho[-el.narm]
					}
				}

			}				
		}
	}

	
	# plot 1: qqnorm residuals
	qqnorm(res)
	qqline(res, col='red')

	# plot 2: density
	check_normaldist(res, legend.label=sprintf('resid(%s)', deparse(substitute(model))) )

	# plot 3: default acf plot
	if(!is.null(res.rho)){
		val <- acf(res.rho[!is.na(res.rho)], plot=F)$acf
		acf(res, col=ifelse(is.null(res.rho), 'black','darkgray'), 
			ylim=range(c(0, val, 1)),
			main=sprintf('resid(%s)', deparse(substitute(model))))
		lines(0:(length(val)-1), val, type='h', lwd=3, col=alpha('red', .5), lend=1, xpd=TRUE)
	}else{
		message("No AR1 model included.")
		acf(res, col=ifelse(is.null(res.rho), 'black','darkgray'), 
			main=sprintf('resid(%s)', deparse(substitute(model))))
	}


	if(!is.null(split_by)){
		# acf split
		if(!is.null(res.rho)){
			val <- acf_plot(res.rho, plot=F, split_by=split_by)
			acf_plot(res, col=ifelse(is.null(res.rho), 'black','darkgray'), split_by=split_by, 
				ylim=range(c(0, val, 1)),
				main=sprintf('resid(%s) - time series', deparse(substitute(model))), ylab='ACF')
			lines(0:(length(val)-1), val, type='h', lwd=3, col=alpha('red', .5), lend=1, xpd=TRUE)
		}else{
			acf_plot(res, col='black', split_by=split_by, 
				main=sprintf('resid(%s) - time series', deparse(substitute(model))), ylab='ACF')
		}
	}else{
		message("No predictors specified to split the residuals. Last plot is canceled.")
	}
}

