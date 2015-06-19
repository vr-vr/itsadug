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


#' Utility unction.
#' 
#' @param model1 A fitted regression model (using lm, glm, gam, or bam).
#' @param model2 A fitted regression model (using lm, glm, gam, or bam).
#' @return A list with model terms that are not shared by both models.
#' @author Jacolien van Rij
#' @examples 
#' data(simdat)
#'
#' # Fit simple GAM model:
#' gam1 <- bam(Y ~ s(Time), data=simdat)
#' gam2 <- bam(Y ~ Group+s(Time), data=simdat)
#' diff_terms(gam1, gam2)
#'
#' @family utility functions

diff_terms <- function(model1, model2){

	fa <- c()
	fb <- c()

	get_formula <- function(x){
		if("gam" %in% class(x)){
			return( attr(terms(x$formula), "term.labels") )
		}else{
			stop(sprintf("This function does not work for model of class %s.", class(x)))
		}
	}

	fa <- get_formula(model1)
	fb <- get_formula(model2)

	d1 <- fa[!fa %in% fb]
	d2 <- fb[!fb %in% fa]

	out <- list()
	out[[deparse(substitute(model1))]] <- d1
	out[[deparse(substitute(model2))]] <- d2
	return(out)
}

#' Utility unction.
#' 
#' @param text A text string (smooth term label) that needs to be converted 
#' to a regular expression. 
#' @return A regular expression string.
#' @author Jacolien van Rij
#' @examples
#' data(simdat)
#' # Model for illustrating coefficients:
#' m0 <- bam(Y ~ s(Time) + s(Subject, bs='re') 
#' + s(Time, Subject, bs='re'), data=simdat)
#' 
#' # get all coefficients:
#' coef(m0)
#' # to get only the Subject intercepts:
#' coef(m0)[grepl(convertNonAlphanumeric("s(Subject)"), names(coef(m0)))]
#' # to get only the Subject slopes:
#' coef(m0)[grepl(convertNonAlphanumeric("s(Time,Subject)"), names(coef(m0)))]
#'
#' @family utility functions

convertNonAlphanumeric <- function(text){
	return( gsub("([^a-zA-Z0-9])", "\\\\\\1", as.character(text)) )
}

