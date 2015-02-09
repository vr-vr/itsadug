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
#'
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
#' @param model A GAMM model, resulting from the functions
#' \code{\link[mgcv]{gam}} or \code{\link[mgcv]{bam}}.
#' @param  AR_start ADD DESCRIPTION
#' @param split_by ADD DESCRIPTION
#' @author Jacolien van Rij

check_resid <- function(model, AR_start = NULL, split_by=NULL){
	par(mfrow=c(2,2), cex=1.1)

	el.narm <- missing_est(model)
	res <- resid(model)
	res.rho <- NULL

	if("lm" %in% class(model)){
		if(!is.null(AR_start)){
			if(length(AR_start) > length(res)){
				AR_start <- AR_start[-el.narm]
			}
		}
		if(!is.null(model$AR1.rho)){
			if(model$AR1.rho > 0){
				suppressWarnings( res.rho <- resid_gam(model, AR_start=AR_start, incl_na=TRUE) )
			}
		}
		if(!is.null(split_by)){
			if(length(split_by[[1]]) > length(res)){
				suppressWarnings( split_by <- lapply(split_by, function(x){x[-el.narm]}) )
			}
		}

	}
	
	# plot 1: qqnorm residuals
	qqnorm(res)
	qqline(res, col='red')

	# plot 2: density
	check_normaldist(res, legend.label=sprintf('resid(%s)', deparse(substitute(model))) )

	# plot 3: default acf plot
	acf(res, col=ifelse(is.null(res.rho), 'black','darkgray'), 
		main=sprintf('resid(%s)', deparse(substitute(model))))

	if(!is.null(res.rho)){
		val <- acf(res.rho[!is.na(res.rho)], plot=F)$acf
		lines(0:(length(val)-1), val, type='h', lwd=3, col=alpha('red', .5), lend=1, xpd=TRUE)
	}


	# acf split
	acf_plot(res, col=ifelse(is.null(res.rho), 'black','darkgray'), split_by=split_by, 
		main=sprintf('resid(%s) - time series', deparse(substitute(model))), ylab='ACF')
	if(!is.null(res.rho)){
		val <- acf_plot(res.rho, plot=F, split_by=split_by)
		lines(0:(length(val)-1), val, type='h', lwd=3, col=alpha('red', .5), lend=1, xpd=TRUE)
	}
}

