#' Utility unction.
#' 
#' @param model A fitted regression model (using lm, glm, gam, or bam).
#' @return The indices of the data that were not fitted by the model.
#' @author Jacolien van Rij
#' @examples 
#' data(simdat)
#'
#' # Add missing values:
#' set.seed(123)
#' simdat[sample(nrow(simdat), size=20),]$Y <- NA
#' # Fit simple linear model:
#' lm1 <- lm(Y ~ Time, data=simdat)
#' na.el <- missing_est(lm1)
#' length(na.el)
#'
#' @family utility functions

missing_est <- function(model){
	if("lm" %in% class(model)){
		el <- unique(model$na.action)
		if(length(el) > 0){
			return(sort(el))
		}else{
			return(el)
		}
	}else{
		stop("This method is currently only implemented for lm, glm, gam, and bam models.")
	}
}
