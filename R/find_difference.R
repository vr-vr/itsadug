#' Find the regions in which the smooth is significantly different from zero.
#' 
#' @export
#' @param mean A vector with smooth predictions.
#' @param se A vector with the standard error on the smooth predictions.
#' @param f A number to multiply the \code{se} with, to convert the \code{se} 
#' into confidence intervals. Defaults to 1.96 (95\% CI), use 2.58 for 99\%CI.
#' @param xVals Optional vector with x values for the smooth. 
#' When \code{xVals} is provided, the regions are returned in terms of x-
#' values, otherwise as indices.
#' @return The function returns a list with start points of each region 
#' (\code{start}) and end points of each region (\code{end}). The logical 
#' \code{xVals} indicates whether the returned values are on the x-scale 
#' (TRUE) or indices (FALSE).
#' @author Jacolien van Rij
#'
#' @family Utility functions for plotting

find_difference <- function(mean, se, f=1.96, xVals = NULL) {
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
 
