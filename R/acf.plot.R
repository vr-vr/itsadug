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
#' @param cond Named list with a selection of the time series events
#' specified in \code{split_by}. Default is NULL, indicating that 
#' all time series are being processed, rather than a selection.
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
#' # see the vignette for examples:
#' vignette("acf", package="itsadug")
#'
#' @family functions for model criticism

acf_plot <- function(x, split_by = NULL, max_lag = NULL, plot = TRUE, 
    fun = mean, cond=NULL, 
    return_all = FALSE, ...) {
    
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

    xname <- deparse(substitute(x))
    plotci <- FALSE
    
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
            # warning(sprintf("No argument for split provided. %s is treated as a single time series.\n", deparse(substitute(x))))
            split_by <- as.factor(rep("Factor", length(x)))
            plotci <- TRUE
        }

        # check condition:
        if(!is.null(cond)){
            if(!is.null(split_by)){
                el <- 1:length(split_by[[1]])
                for(i in names(cond)){
                    if(i %in% names(split_by)){
                        el <- intersect(el, which(split_by[[i]] %in% cond[[i]]))
                    }else{
                        warning(sprintf('Predictor %s not specified in split_by. cond will be ignored.', i))
                    }
                }
                if(length(el) > 0){
                    x <- x[el]
                    for(i in names(split_by)){
                        if(is.factor(split_by[[i]])){
                            split_by[[i]] <- droplevels( split_by[[i]][el] )
                        }else{
                            split_by[[i]] <- split_by[[i]][el] 
                        }
                    }
                }else{
                    warning("Specified conditions not found in values of split_by. cond will be ignored.")
                }
                
            }else{
                warning('Split_by is empty, therefore cond will be ignored. Specify time series in cond for selecting specific time series.')
            }

        }
        
        # split x, and calculate average acf:
        splitdat <- split(x, f = split_by, drop = T)
        acfn <- lapply(splitdat, FUN = function(x) {
            x.rmna <- x[!is.na(x)]
            return( length(x.rmna) )
        })

        splitacf <- lapply(splitdat, FUN = function(x) {
            x.rmna <- x[!is.na(x)]
            return( acf(x.rmna, plot = F)$acf )
        })
        len <- max_lag
        if (is.null(len)) {
            len <- max(unlist(lapply(splitacf, FUN = length)), na.rm = T)
        }
        splitacf <- lapply(splitacf, FUN = function(x, max = len) {
            if (length(x) < max) {
                return(c(x, rep(NA, max - length(x))))
            } else if (length(x) > max) {
                return(x[1:max])
            } else {
                return(x)
            }
        })       

        # create wide format dfr:
        allacf <- as.data.frame(do.call("rbind", splitacf))
        names(allacf) <- (1:ncol(allacf)) - 1
        avgacf <- apply(allacf, 2, FUN = fun)  #apply(allacf, 2, FUN=fun, na.rm=T)
        
        if (plot) {
            # set plot arguments
            plot_default <- list(main = sprintf("ACF of %s", xname) , 
                xlab = "Lag", ylab = ifelse(plotci==TRUE, "ACF function (per time series)", "ACF"),
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

            if(plotci){
                ci <- -(1/acfn[[1]])+2/sqrt(acfn[[1]])
                abline(h=c(-1,1)*ci, lty=2, col='blue')   
            }
            # }else{
            #     tmpn <- length(allacf[!is.na(allacf)])
            #     ci <- -(1/tmpn)+2/sqrt(tmpn)
            #     abline(h=c(-1,1)*ci, lty=2, col='blue')
            # }
            abline(h = 0)
        }
        
        # set output:
        acf_out <- avgacf
        if (return_all) {
            # create long format dfr:
            dfracf <- do.call('rbind',
                mapply(function(x, y, z){
                    data.frame(acf=x, 
                        lag=0:(length(x)-1),
                        n = rep(y, length(x)),
                        event = rep(z, length(x))) 
                    }, splitacf, acfn, names(splitacf), SIMPLIFY=FALSE, USE.NAMES=FALSE) )
            dfracf$ci <- -(1/dfracf$n) + 2/sqrt(dfracf$n)

            # add event info:
            events <- as.data.frame(split_by)
            events$event <- apply(events, 1, function(x){gsub(" ", "", paste(x, collapse="."), fixed=TRUE)})
            events <- events[!duplicated(events),]

            dfracf <- merge(dfracf, events, by='event', all.x=TRUE, all.y=FALSE)
            acfn <- do.call('rbind', lapply(names(acfn), 
                function(x){data.frame(n=acfn[[x]], event=x)}))
           
            acf_out <- list(acf = avgacf, acftable = allacf, 
                dataframe=dfracf, n=acfn, series = deparse(substitute(x)), FUN = fun)
            
        }
        invisible(acf_out)        
    } else {
        stop(sprintf("Not sufficient data to plot ACF: %s has %d elements.\n", deparse(substitute(x)), length(x)))
    }
}
 
