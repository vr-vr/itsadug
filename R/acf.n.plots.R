#' Generate N ACF plots of individual or aggregated time series.
#' 
#' @export
#' @param x A vector with time series data, typically residuals of a 
#' regression model.
#' @param n The number of plots to generate.
#' @param split_by List of vectors (each with equal length as \code{x}) that 
#' group the values of \code{x} into trials or timeseries events. 
#' Generally other columns from the same data frame.
#' @param cond Named list with a selection of the time series events
#' specified in \code{split_by}. Default is NULL, indicating that 
#' all time series are being processed, rather than a selection.
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
#' @param plot Logical: whether or not to produce plot. Default is TRUE.
#' @param ... Other arguments for plotting, see \code{\link[graphics]{par}}.
#' @return \code{n} ACF plots providing information about the autocorrelation 
#' in \code{x}.
#' @author Jacolien van Rij, R. Harald Baayen
#' @seealso Use \code{\link[stats]{acf}} for the original ACF function, 
#' and \code{\link{acf_plot}} for an ACF that takes into account time series 
#' in the data.
#' @examples
#' data(simdat)
#' # Separate ACF for each time series:
#' acf_n_plots(simdat$Y, split_by=list(simdat$Subject, simdat$Trial))
#' # Average ACF per participant:
#' acf_n_plots(simdat$Y, split_by=list(simdat$Subject))
#' \dontrun{
#' # Data treated as single time series. Plot is added to current window.
#' # Note: 1 time series results in 1 plot.
#' acf_n_plots(simdat$Y, add=TRUE)
#' # Plot 4 ACF plots doesn't work without splitting data:
#' acf_n_plots(simdat$Y, add=TRUE, n=4)
#' 
#' # Plot ACFs of 4 randomly selected time series:
#' acf_n_plots(simdat$Y, random=TRUE, n=4, add=TRUE, 
#'     split_by=list(simdat$Subject, simdat$Trial))
#'
#' # See also the vignette for an example:
#' vignette(topic="plotfunctions", package="itsadug")
#' }
#'
#' #---------------------------------------------
#' # When using model residuals
#' #---------------------------------------------
#' \dontrun{
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
#' }
#' @family functions for model criticism

acf_n_plots <- function(x, n = 5, split_by = NULL, 
    cond = NULL, max_lag = NULL, fun = mean, plot=TRUE,
    random = F, mfrow = NULL, add=FALSE, ...) {
    
    # get acf data:
    suppressWarnings( acfdat <- acf_plot(x, split_by = split_by, 
        cond=cond, max_lag = max_lag, fun = fun, 
        plot = FALSE, return_all = TRUE) )

    if (!nrow(acfdat$acftable) >= n) {
        warning(sprintf("Number of time series in the data (%d) is smaller than n (%d). N is reduced to %d.\n", 
            nrow(acfdat$acftable), n, nrow(acfdat$acftable)))
        n <- nrow(acfdat$acftable)
    }


    # get lag 1A vector: 
    lag1.all <- acfdat$acftable$"1"
    rn <- rownames(acfdat$acftable)
    lag1 <- lag1.all[!is.na(lag1.all)]
    nevents <- acfdat$n

    findNum <- rep(1, n)
    out <- list()
    
    if (random) {
        findNum <- sample(nrow(acfdat$acftable), size = n)
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


        for (i in 1:n) {
            findNum <- c(findNum, findClosestElement(q[i], lag1))
            if(i > 1){
                out[[i-1]] <- list(quantile=c(q[i-1], q[i]),
                    elements = 
                        data.frame(event=rn[which(lag1.all >= q[i-1] &lag1.all < q[i])],
                            lag1 = lag1.all[which(lag1.all >= q[i-1] &lag1.all < q[i])],
                            stringsAsFactors=FALSE) )
            }
         }
        out[['quantiles']] <- q
        cat("Quantiles to be plotted:\n")
        print(q)
    }
    
    if(plot){    
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
        ylim <- range(acfdat$acftable[findNum, ], na.rm = T)
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
        
        for (j in 1:min(n, nrow(acfdat$acftable)) ) {
            main <- NULL
            if(is.null(main.new)){
                main <- paste("ACF of", row.names(acfdat$acftable[findNum[j], ]))
            }else{
                if(length(main.new) >= j){
                    main <- main.new[j]
                }else{
                    main <- main.new[1]
                }
            }

            eval(parse(text=sprintf(
                "plot(as.numeric(colnames(acfdat$acftable)), acfdat$acftable[findNum[j], ], type = 'h', 
                main=main, xlab=xlab, ylab=ylab, ylim = ylim, %s)", other)))

            ci <- -(1/nevents[findNum[j],'n'])+2/sqrt(nevents[findNum[j],'n'])
            abline(h=c(-1,1)*ci, lty=2, col='blue')
            abline(h = 0)
        }
    }

    invisible(out)    
} 
