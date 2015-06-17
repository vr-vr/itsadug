#' Simulated time series data.
#'
#' A dataset containing the sine wave data with random noise added. 
#'
#' @format A data frame with 75600 rows and 6 variables:
#' \describe{
#'   \item{\code{Group}}{Age group of participants: Adults or Children.}
#'   \item{\code{Time}}{Time, time measure from start of each time series.}
#'   \item{\code{Trial}}{Trial in the experiment, centered around zero.}
#'   \item{\code{Condition}}{Continuous variable, ranging from -1 to 4. 
#' For example, stimulus onset asynchrony.}
#'   \item{\code{Subject}}{Code for individual participants.}
#'   \item{\code{Y}}{Time series measure. 
#' Similar to pupil size, sensor position, or voltage.}
#' }
#' @author Jacolien van Rij
"simdat"
