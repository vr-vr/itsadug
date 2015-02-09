#' Generate N ACF plots of individual or aggregated time series.
#' 
#' @export
#' @aliases acf.n.plots
#' @param x A vector with time series data, typically residuals of a 
#' regression model.
#' @param n The number of plots to generate.
#' @param split_by List of vectors (each with equal length as \code{x}) that 
#' group the values of \code{x} into trials or timeseries events. 
#' Generally other columns from the same data frame.
#' @param max_lag Maximum lag at which to calculate the acf. 
#' Default is the maximum for the longest time series.
#' @param fun The function used when aggregating over time series 
#' (depending on the value of \code{split_by}).
#' @param random Logical: determine randomly which \code{n} (aggregated) time 
#' series are plotted, or use the \code{\link{quantile}} function to find a 
#' range of different time series to plot. Default is FALSE (not random).
#' @param mfrow A vector of the form c(nr, nc). The figures will be drawn in 
#' an nr-by-nc array on the device by rows.
#' @param add Logical: whether to add the plots to an exiting plot window or 
#' not. 
#' Default is FALSE.
#' @param ... Other arguments for plotting, see \code{\link[graphics]{par}}.
#' @return \code{n} ACF plots providing information about the autocorrelation 
#' in \code{x}.
#' @author Jacolien van Rij, R. Harald Baayen
#' @seealso Use \code{\link[stats]{acf}} for the original ACF function, 
#' and \code{\link{acf_plot}} for an ACF that takes into account time series 
#' in the data.
#' @examples
#' data(simdat)
#' acf_n_plots(simdat$Y, split_by=list(simdat$Subject, simdat$Trial))
#' acf_n_plots(simdat$Y, split_by=list(simdat$Subject))
#' acf_n_plots(simdat$Y, add=TRUE)
#' acf_n_plots(simdat$Y, add=TRUE, n=1)
#' acf_n_plots(simdat$Y, add=TRUE, n=1, random=TRUE)
#' acf_n_plots(simdat$Y, random=TRUE, n=4, add=TRUE, 
#'     split_by=list(simdat$Subject, simdat$Trial))
#' # see also the vignette plotfunctions for an example:
#' vignette(topic="plotfunctions")
#' 
#' #---------------------------------------------
#' # When using model residuals
#' #---------------------------------------------
#' 
#' # add missing values to simdat:
#' simdat[sample(nrow(simdat), 15),]$Y <- NA
#' # simple linear model:
#' m1 <- lm(Y ~ Time, data=simdat)
#'
#' # This will generate an error:
#' # acf_n_plots(resid(m1), split_by=list(simdat$Subject, simdat$Trial))
#'
#' # This should work:
#' el.na <- missing_est(m1)
#' acf_n_plots(resid(m1), 
#'      split_by=list(simdat[-el.na,]$Subject, simdat[-el.na,]$Trial))
#'
#' # This should also work:
#' simdat$res <- NA
#' simdat[!is.na(simdat$Y),]$res <- resid(m1)
#' acf_n_plots(simdat$res, split_by=list(simdat$Subject, simdat$Trial))
#' 
#' @family functions for model criticism

acf_n_plots <- function(x, n = 5, split_by = NULL, max_lag = NULL, fun = mean, 
    random = F, mfrow = NULL, add=FALSE, ...) {
    
    # get acf data:
    suppressWarnings( acfdat <- acf_plot(x, split_by = split_by, 
        max_lag = max_lag, fun = fun, plot = F, return_all = T) )

    if (!nrow(acfdat$acf_split) >= n) {
        warning(sprintf("Number of time series in the data (%d) is smaller than n (%d).\n", nrow(acfdat$acf_split), n))
    }

    # get lag 1A vector:
    lag1 <- acfdat$acf_split$"1"
    lag1 <- lag1[!is.na(lag1)]
    
    findNum <- rep(1, n)
    
    if (random) {
        findNum <- sample(nrow(acfdat$acf_split), size = n)
    } else {
        q <- quantile(lag1, probs = seq(0, 1, length = n))
        
        findClosestElement <- function(el, vals, num = TRUE) {
            y <- abs(vals - el)
            y <- min(y)
            if (num) {
                return(which(abs(vals - el) == y)[1])
            } else {
                return(vals[which(abs(vals - el) == y)[1]])
            }
        }
        
        findNum <- c()
        for (i in q) {
            findNum <- c(findNum, findClosestElement(i, lag1))
        }
        cat("Quantiles to be plotted:\n")
        print(q)
    }
    
    cols <- ceiling(sqrt(n))
    rows <- ceiling(n/cols)
    
    if (!is.null(mfrow)) {
        cols = mfrow[2]
        rows = mfrow[1]
    }
    
    if(add==FALSE){
        # dev.new(width = cols * 3, height = rows * 3)
        par(mfrow = c(rows, cols))
        oldmar <- par(mar = rep(3, 4))
    }

    parlist = list(...)
    if('type' %in% names(parlist)){
        parlist[['type']] <- NULL
    }

    xlab <- "Lag"
    ylab <- "ACF"
    ylim <- range(acfdat$acf_split[findNum, ], na.rm = T)
    main.new <- NULL

    if('xlab' %in% names(parlist)){
        xlab <- parlist[['xlab']]
        parlist[['xlab']] <- NULL
    }
    if('ylab' %in% names(parlist)){
        ylab <- parlist[['ylab']]
        parlist[['ylab']] <- NULL
    }
    if('ylim' %in% names(parlist)){
        ylim <- parlist[['ylim']]
        parlist[['ylim']] <- NULL
    }
    if('main' %in% names(parlist)){
        main.new <- parlist[['main']]
        parlist[['main']] <- NULL
    }

    other <- paste(sprintf('%s=%s', names(parlist), parlist), collapse=',')
    
    for (j in 1:min(n, nrow(acfdat$acf_split)) ) {
        main <- NULL
        if(is.null(main.new)){
            main <- paste("ACF of", row.names(acfdat$acf_split[findNum[j], ]))
        }else{
            if(length(main.new) >= j){
                main <- main.new[j]
            }else{
                main <- main.new[1]
            }
        }

        eval(parse(text=sprintf(
            "plot(as.numeric(colnames(acfdat$acf_split)), acfdat$acf_split[findNum[j], ], type = 'h', 
            main=main, xlab=xlab, ylab=ylab, ylim = ylim, %s)", other)))
        abline(h = 0)
    }
    
} 
