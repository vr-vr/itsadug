#' Find the regions in which the smooth is significantly different from zero.
#' 
#' @export
#' @param mean A vector with smooth predictions.
#' @param se A vector with the standard error on the smooth predictions.
#' @param xVals Optional vector with x values for the smooth. 
#' When \code{xVals} is provided, the regions are returned in terms of x-
#' values, otherwise as indices.
#' @param f A number to multiply the \code{se} with, to convert the \code{se} 
#' into confidence intervals. Use 1.96 for 95\% CI and 2.58 for 99\%CI.
#' @return The function returns a list with start points of each region 
#' (\code{start}) and end points of each region (\code{end}). The logical 
#' \code{xVals} indicates whether the returned values are on the x-scale 
#' (TRUE) or indices (FALSE).
#' @examples
#' data(simdat)
#' 
#' # Use aggregate to calculate mean and standard deviation per timestamp:
#' avg <- aggregate(simdat$Y, by=list(Time=simdat$Time),
#'     function(x){c(mean=mean(x), sd=sd(x))})
#' head(avg)
#' # Note that column x has two values in each row:
#' head(avg$x)
#' head(avg$x[,1])
#' 
#' # Plot line and standard deviation:
#' emptyPlot(range(avg$Time), c(-20,20), h0=0)
#' plot_error(avg$Time, avg$x[,'mean'], avg$x[,'sd'], 
#'    shade=TRUE, lty=3, lwd=3)
#'
#' # Show difference with 0:
#' x <- find_difference(avg$x[,'mean'], avg$x[,'sd'], xVals=avg$Time)
#' # Add arrows:
#' abline(v=c(x$start, x$end), lty=3, col='red')
#' arrows(x0=x$start, x1=x$end, y0=-5, y1=-5, code=3, length=.1, col='red')
#' @author Jacolien van Rij
#'
#' @family Utility functions for plotting

find_difference <- function(mean, se, xVals = NULL,f=1) {
    if (length(mean) != length(se)) {
        stop("The vectors mean and se are not equal in length.")
    } else {
        ub <- mean + f*se
        lb <- mean - f*se
        
        n <- which(!(ub >= 0 & lb <= 0))
        if (length(n) == 0) {
            return(NULL)
        } else {
            n_prev <- c(NA, n[1:(length(n) - 1)])
            n_next <- c(n[2:length(n)], NA)
            if (!is.null(xVals) & (length(xVals) == length(mean))) {
                return(list(start = xVals[n[which(is.na(n - n_prev) | (n - n_prev) > 1)]], end = xVals[n[which(is.na(n_next - 
                  n) | (n_next - n) > 1)]], xVals = TRUE))
            } else {
                return(list(start = n[which(is.na(n - n_prev) | (n - n_prev) > 1)], end = n[which(is.na(n_next - n) | (n_next - 
                  n) > 1)], xVals = FALSE))
            }
            
        }
    }
}
 
