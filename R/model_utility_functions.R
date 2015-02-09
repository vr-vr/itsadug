#' Utility unction.
#' 
#' @param model A fitted regression model (using lm, glm, gam, or bam).
#' @return The indices of the data that were not fitted by the model.
#' @author Jacolien van Rij
#' @family utility functions

missing_est <- function(model){
	if("lm" %in% class(model)){
		return(sort( unique(model$na.action)))
	}else{
		stop("This method is currently only implemented for lm, glm, gam, and bam models.")
	}
}
