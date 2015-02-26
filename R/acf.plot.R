#' Generate an ACF plot of an aggregated time series.
#' 
#' @export
#' @param x A vector with time series data, typically residuals of a 
#' regression model. 
#' (See examples for how to avoid errors due to missing values.)
#' @param split_by List of vectors (each with equal length as \code{x}) that 
#' group the values of \code{x} into trials or timeseries events. 
#' Generally other columns from the same data frame.
#' @param max_lag Maximum lag at which to calculate the acf. 
#' Default is the maximum for the longest time series.
#' @param plot Logical: whether or not to plot the ACF. Default is TRUE.
#' @param fun The function used when aggregating over time series 
#' (depending on the value of \code{split_by}).
#' @param return_all Returning acfs for all time series.
#' @param ... Other arguments for plotting, see \code{\link[graphics]{par}}.
#' @return An aggregated ACF plot and / or optionally a list with the aggregated ACF values.
#' @author Jacolien van Rij
#' @seealso Use \code{\link[stats]{acf}} for the original ACF function, 
#' and \code{\link{acf_n_plots}} for inspection of individual time series.
#' @examples
#' data(simdat)
#'
#' # Default acf function:
#' acf(simdat$Y)
#' # Same plot with acf_plot:
#' acf_plot(simdat$Y)
#' # Average of ACFs per time series:
#' acf_plot(simdat$Y, split_by=list(simdat$Subject, simdat$Trial))
#' # Median of ACFs per time series:
#' acf_plot(simdat$Y, split_by=list(simdat$Subject, simdat$Trial), fun=median)
#'
#' # extract value of Lag1:
#' lag1 <- acf_plot(simdat$Y, 
#'    split_by=list(Subject=simdat$Subject, Trial=simdat$Trial), 
#'    plot=FALSE)['1']
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
#' \dontrun{
#' # This will generate an error:
#' acf_plot(resid(m1), split_by=list(simdat$Subject, simdat$Trial))
#' }
#' # This should work:
#' el.na <- missing_est(m1)
#' acf_plot(resid(m1), 
#'      split_by=list(simdat[-el.na,]$Subject, simdat[-el.na,]$Trial))
#'
#' # This should also work:
#' simdat$res <- NA
#' simdat[!is.na(simdat$Y),]$res <- resid(m1)
#' acf_plot(simdat$res, split_by=list(simdat$Subject, simdat$Trial))
#'
#' # see also the vignette for an example:
#' vignette(topic="plotfunctions", package="itsadug")
#' @family functions for model criticism

acf_plot <- function(x, split_by = NULL, max_lag = NULL, plot = TRUE, fun = mean, return_all = FALSE, ...) {
    
    # check x:
    if (!is.vector(x)) {
        if (dim(x)[1] == 1) {
            x <- as.vector(x[1, ])
        } else if (dim(x)[2] == 1) {
            x <- as.vector(x[, 1])
        } else {
            stop(sprintf("%s is not a vector (dim: %d x %d).\n", deparse(substitute(x)), dim(x)[1], dim(x)[2]))
        }
    }
    
    # check length x
    if (length(x) > 1) {
        # find missing values:
        n <- which(is.na(x))
        
        # check split factors:
        if (!is.null(split_by)) {
            for (i in 1:length(split_by)) {
                if (length(split_by[[i]]) != length(x)) {
                  name_split_i <- ifelse(!is.null(names(split_by)[i]), names(split_by)[i], sprintf("split_by[[%d]]", i))
                  errormessage <- sprintf("Split factor %s is not of same length as %s: %s has %d elements, %s has %d elements.\n
                    See help(acf_plot) for examples how to avoid this error.", 
                    name_split_i, 
                    deparse(substitute(x)), 
                    deparse(substitute(x)), 
                    length(x), 
                    name_split_i, 
                    length(split_by[[i]]))
                  stop(errormessage)
                }
            }
        } else {
            warning(sprintf("No argument for split provided. %s is treated as a single time series.\n", deparse(substitute(x))))
            split_by <- as.factor(rep("Factor", length(x)))
        }
        
        # split x, and calculate average acf:
        splitdat <- split(x, f = split_by, drop = T)
        allacf <- lapply(splitdat, FUN = function(x) {
            x.rmna <- x[!is.na(x)]
            return( acf(x.rmna, plot = F)$acf )
        })
        len <- max_lag
        if (is.null(len)) {
            len <- max(unlist(lapply(allacf, FUN = length)), na.rm = T)
        }
        allacf <- lapply(allacf, FUN = function(x, max = len) {
            if (length(x) < max) {
                return(c(x, rep(NA, max - length(x))))
            } else if (length(x) > max) {
                return(x[1:max])
            } else {
                return(x)
            }
        })
        allacf <- as.data.frame(do.call("rbind", allacf))
        names(allacf) <- (1:ncol(allacf)) - 1
        avgacf <- apply(allacf, 2, FUN = fun)  #apply(allacf, 2, FUN=fun, na.rm=T)
        
        if (plot) {
            # set plot arguments
            plot_default <- list(main = sprintf("ACF of %s", deparse(substitute(x))), xlab = "Lag", ylab = "ACF function (per time series)", 
                ylim = c(min(min(avgacf, na.rm=TRUE), 0), max(max(avgacf, na.rm=TRUE), 1)), col = "black", type = "h")
            
            plot_args <- list(...)
            
            # merge plot info
            plot_args_def <- c()
            for (pa in names(plot_default)[!names(plot_default) %in% names(plot_args)]) {
                value <- plot_default[[pa]]
                if (length(value) > 1) {
                  value <- sprintf("c(%s)", paste(value, collapse = ","))
                  plot_args_def <- c(plot_args_def, paste(pa, value, sep = "="))
                } else {
                  plot_args_def <- c(plot_args_def, paste(pa, ifelse(is.character(value), sprintf("'%s'", value), value), 
                    sep = "="))
                }
            }
            
            if (length(plot_args_def) > 0) {
                eval(parse(text = paste("plot(0:(len-1), avgacf, ", paste(plot_args_def, collapse = ","), ", ...)")))
            } else {
                plot(0:(len - 1), avgacf, ...)
            }
            abline(h = 0)
        }
        
        # set output:
        acf_out <- avgacf
        if (return_all) {
            acf_out <- list(acf = avgacf, acf_split = allacf, series = deparse(substitute(x)), FUN = fun)
            
        }
        
        invisible(acf_out)
        
        
    } else {
        stop(sprintf("Not sufficient data to plot ACF: %s has %d elements.\n", deparse(substitute(x)), length(x)))
    }
}
 
