#' Function for comparing two GAMM models.
#' 
#' @export
#' @import mgcv
#' @import stats
#' @param model1 First model.
#' @param model2 Second model.
#' @param print.output Logical: whether or not to print the output. 
#' By default controlled globally in the package options:  
#' If the function \code{\link{infoMessages}} is set to TRUE, the output 
#' will be automatically printed.
#' Could be also set by explicitly providing TRUE or FALSE. See notes.
#' @details
#' 
#' As an Chi-Square test is performed on two times the difference in 
#' minimized smoothing parameter selection score (GCV, fREML, REML, ML), 
#' and the difference in degrees of freedom specified in the model. 
#' The degrees of freedom of the model terms are the sum of
#' 1) the number of estimated smoothing parameters for the model, 
#' 2) number of parametric (non-smooth) model terms including the intercept, 
#' and 3) the sum of the penalty null space dimensions of each smooth object.
#'
#' This method is preferred over other functions such as \code{\link{AIC}} for 
#' models that include an AR1 model or random effects (especially nonlinear 
#' random smooths using \code{bs="fs"}). CompareML also reports the AIC 
#' difference, but that value should be treated with care.
#' 
#' Note that the Chi-Square test will result in a very low p-value
#' when the difference in degrees of freedom approaches zero. Use common sense 
#' to determine if the difference between the two models is meaningful. 
#' A warning is presented when the difference in score is smaller 
#' than 5.
#'
#' The order of the two models is not important.
#' Model comparison is only implemented for the methods GCV, fREML, REML, and ML.
#'
#' @section Notes:
#' If no output is provided in the command window, set info messages to TRUE: 
#' \code{infoMessages("on")} and try again.
#' For suppressing the output and all warnings, set infoMessages to FALSE 
#' (\code{infoMessages("on")} ) and use the function 
#' \code{\link{suppressWarnings}} to suppress warning messages.
#' @return Optionally returns the Chi-Square test table.
#' @author Jacolien van Rij. With many thanks to Simon N. Wood for his feedback.
#' @seealso For models without AR1 model or random effects \code{\link{AIC}} can be used.
# help function
#' @examples
#' data(simdat)
#'
#' \dontrun{
#' infoMessages("on")
#' # some arbitrary models:
#' m1 <- bam(Y~Group + s(Time, by=Group), method="fREML", data=simdat)
#' m2 <- bam(Y~Group + s(Time), method="fREML", data=simdat)
#' 
#' compareML(m1, m2)
#'
#' m3 <- bam(Y~Group + s(Time, by=Group, k=25), method="fREML", 
#'     data=simdat)
#' compareML(m1, m3)
#' 
#' }


compareML <- function(model1, model2,
    print.output=getOption('itsadug_print')) {
    # check gam or bam model:
    if((!"gam" %in% class(model1)) | (!"gam" %in% class(model2))){
        stop("Models are not gam objects (i.e., build with bam()/gam()).")
    }
    
    # check whether models are comparable:
    if (model1$method != model2$method) {
        stop(sprintf("Models are incomparable: method model1 = %s, method model2 = %s", model1$method, model2$method))
    }
    
    
    type <- model1$method
    
    ml1 <- model1$gcv.ubre[1]
    ml2 <- model2$gcv.ubre[1]
    
    ### OLD METHOD, SIMON SAYS NOT OK! ### edf1 <- sum(model1$edf) edf2 <- sum(model2$edf) NEW METHOD: ###

    ndf1 <- length(model1$sp) + model1$nsdf + ifelse(length(model1$smooth)>0,
        sum(sapply(model1$smooth, FUN = function(x) {
            x$null.space.dim
            }, USE.NAMES = FALSE)), 0)
    ndf2 <- length(model2$sp) + model2$nsdf + ifelse(length(model2$smooth)>0,
        sum(sapply(model2$smooth, FUN = function(x) {
            x$null.space.dim
            }, USE.NAMES = FALSE)), 0)
    
    if (! model1$method %in%  c("fREML", "REML", "ML")) {
        type <- "AIC"
        ml1 <- AIC(model1)
        ml2 <- AIC(model2)
        
	    ndf1 <- length(model1$sp) + model1$nsdf + ifelse(length(model1$smooth)>0,
	        sum(sapply(model1$smooth, FUN = function(x) {
	            x$null.space.dim
	            }, USE.NAMES = FALSE)), 0)
	    ndf2 <- length(model2$sp) + model2$nsdf + ifelse(length(model2$smooth)>0,
	        sum(sapply(model2$smooth, FUN = function(x) {
	            x$null.space.dim
	            }, USE.NAMES = FALSE)), 0)
        warning(sprintf("\nCompareML is not implemented for smoothing parameter estimation method %s. AIC scores are used for model comparison. Consider running the model with REML, fREML, or ML as method.\n-----\n", model1$method))
    }
    
    
    # pchisq(4, .5, lower.tail=F) # p < .1 pchisq(-4, .5, lower.tail=F) # p = 1 pchisq(4, -.5, lower.tail=F) # NaN
    
    # Book keeping;
    info <- sprintf("%s: %s\n\n%s: %s\n", deparse(substitute(model1)), deparse(model1$formula),
        deparse(substitute(model2)), deparse(model2$formula))

    if(print.output){
        cat(info)
    }
    
    out      <- NULL
    advice   <- NULL
    warning  <- NULL
    
    # if (type != 'AIC') {
    # Situation 1: model 1 has lower score, but model 2 has lower df. Is it significantly better model than model 2?
	# Situation 0: equal df
	if(abs(round(ndf2 - ndf1)) < .5){
		if( ml1 < ml2){
            advice <- sprintf("\nModel %s preferred: lower %s score (%.3f), and equal df (%.3f).\n-----\n", 
            	deparse(substitute(model1)), 
                type, ml2 - ml1, ndf2 - ndf1)

            out <- data.frame(Model = c(deparse(substitute(model2)), deparse(substitute(model1))), 
            	Score = c(ml2, ml1), 
            	Edf = c(ndf2, ndf1), 
            	Difference = c("", sprintf("%.3f", ml2 - ml1)), 
            	Df = c("", sprintf("%.3f", ndf1 - ndf2)),
                p.value = c("",NA),
                Sign. = c("", ""))
		}else{
            advice <- sprintf("\nModel %s preferred: lower %s score (%.3f), and equal df (%.3f).\n-----\n", 
            	deparse(substitute(model2)), 
                type, ml1 - ml2, ndf1 - ndf2)
            out <- data.frame(Model = c(deparse(substitute(model1)), deparse(substitute(model2))), 
            	Score = c(ml1, ml2), 
            	Edf = c(ndf1, ndf2),
            	Difference = c("", sprintf("%.3f", ml2 - ml1)), 
            	Df = c("", sprintf("%.3f", ndf1 - ndf2)),
                p.value = c("",NA),
                Sign. = c("", ""))
		}
	# Situation 1: model 1 has lower score, but model 2 has lower df. Is it significantly better model than model 2?
    }else if ((ml1 < ml2) & (ndf2 < ndf1)) {
        
        # twice the amount of difference in likelihood
        h1 <- pchisq(2 * (ml2 - ml1), abs(ndf1 - ndf2), lower.tail = F)
        
        out <- data.frame(Model = c(deparse(substitute(model2)), deparse(substitute(model1))), 
        	Score = c(ml2, ml1), 
        	Edf = c(ndf2, ndf1), 
            Chisq = c("", sprintf("%.3f", ml2 - ml1)), 
            Df = c("", sprintf("%.3f", abs(ndf1 - ndf2))), 
            p.value = c("", ifelse(h1 < 2e-16, sprintf(" < 2e-16"), 
            	ifelse(h1 < 0.001, sprintf("%.3e", h1), 
            	ifelse(h1 < 0.01, sprintf("%.3f", h1), 
            	ifelse(h1 < 0.05, sprintf("%.3f", h1), sprintf("%.3f", h1)))))), 
            Sig. = c("", ifelse(h1 < 0.001, sprintf("***", h1), 
            	ifelse(h1 < 0.01, sprintf("** ", h1), 
            	ifelse(h1 < 0.05, sprintf("*  ", h1), sprintf("   ", h1))))))
        
    # Situation 2: model 2 has lower score, but model 1 has lower df. Is it significantly better model than model 1?
    } else if ((ml2 < ml1) & (ndf1 < ndf2)) {
        
        h1 <- pchisq(2 * (ml1 - ml2), abs(ndf1 - ndf2), lower.tail = F)
        
        out <- data.frame(Model = c(deparse(substitute(model1)), deparse(substitute(model2))), 
        	Score = c(ml1, ml2), 
        	Edf = c(ndf1, ndf2), 
        	Chisq = c("", sprintf("%.3f", ml1 - ml2)), 
        	Df = c("", sprintf("%.3f", abs(ndf1 - ndf2))), 
        	p.value = c("", ifelse(h1 < 2e-16, sprintf(" < 2e-16"), 
        		ifelse(h1 < 0.001, sprintf("%.3e", h1), 
            	ifelse(h1 < 0.01, sprintf("%.3f", h1), 
            	ifelse(h1 < 0.05, sprintf("%.3f", h1), sprintf("%.3f", h1)))))), 
            Sig. = c("", ifelse(h1 < 0.001, sprintf("***", h1), 
            	ifelse(h1 < 0.01, sprintf("** ", h1), 
            	ifelse(h1 < 0.05, sprintf("*  ", h1), sprintf("   ", h1))))))
        
    # Situation 3: model 1 has lower score, and also lower df.
    } else if ((ml1 < ml2) & (ndf1 < ndf2)) {
        advice <- sprintf("\nModel %s preferred: lower %s score (%.3f), and lower df (%.3f).\n-----\n", 
        	deparse(substitute(model1)), 
        	type, ml2 - ml1, ndf2 - ndf1)
        out <- data.frame(Model = c(deparse(substitute(model2)), deparse(substitute(model1))), 
        	Score = c(ml2, ml1), 
        	Edf = c(ndf2, ndf1), 
        	Difference = c("", sprintf("%.3f", ml2 - ml1)), 
        	Df = c("", sprintf("%.3f", ndf1 - ndf2)),
            p.value = c("",NA),
            Sign. = c("", ""))
        
    # Situation 4: model 2 has lower score, and also lower df.
    } else if ((ml2 < ml1) & (ndf2 < ndf1)) {
        advice <- sprintf("\nModel %s preferred: lower %s score (%.3f), and lower df (%.3f).\n-----\n", 
        	deparse(substitute(model2)), 
            type, ml1 - ml2, ndf1 - ndf2)
        out <- data.frame(Model = c(deparse(substitute(model1)), deparse(substitute(model2))), 
        	Score = c(ml1, ml2), 
        	Edf = c(ndf1, ndf2), 
        	Difference = c("", sprintf("%.3f", ml2 - ml1)), 
        	Df = c("", sprintf("%.3f", ndf1 - ndf2)),
            p.value = c("",NA),
            Sign. = c("", ""))
    # Other cases:
    } else {
        advice <- "No preference:\n-----\n"
        out <- data.frame(Model = c(deparse(substitute(model1)), deparse(substitute(model2))), 
        	Score = c(ml1, ml2), Edf = c(ndf1, ndf2), 
        	Difference = c("", sprintf("%.3f", ml2 - ml1)), 
        	Df = c("", sprintf("%.3f", ndf1 - ndf2)),
            p.value = c("",NA),
            Sign. = c("", ""))
    }

    if(print.output){
        if(is.null(advice)){
            cat(sprintf("\nChi-square test of %s scores\n-----\n", type))
        }else{
            cat(advice)
        }
        if(is.na(out$p.value[2])){
            print(out[,1:5])
        }else{
            print(out)
        }
        
        cat("\n")
    }
    
    if (type != "AIC") {
        if (is.null(model1$AR1.rho)) {
            rho1 = 0
        } else {
            rho1 = model1$AR1.rho
        }
        
        if (is.null(model2$AR1.rho)) {
            rho2 = 0
        } else {
            rho2 = model2$AR1.rho
        }
        
        if (rho1 == 0 & rho2 == 0) {
            # AIC is useless for models with rho
            if (AIC(model1) == AIC(model2)) {
              warning <- sprintf("AIC difference: 0.\n\n")
            } else {
              warning <- sprintf("AIC difference: %.2f, model %s has lower AIC.\n\n", 
                AIC(model1) - AIC(model2), 
                ifelse(AIC(model1) >= AIC(model2), 
                    deparse(substitute(model2)), 
                    deparse(substitute(model1))))
            }
            if(print.output){
                cat(warning)
            }
        } else {
            warning(warning <- sprintf(" AIC is not reliable, because an AR1 model is included (rho1 = %f, rho2 = %f). ", 
                rho1, rho2))
        }
    }
    
    if (abs(ml1 - ml2) <= 5) {
        warning(sprintf("Only small difference in %s...\n", type))
    }
    
    invisible( list(method=type,
    	m1=list(Model=model1$formula, Score=ml1, Df=ndf1),
    	m2=list(Model=model2$formula, Score=ml2, Df=ndf2),
    	table = out,
        advice = ifelse(is.null(advice), NA, advice),
        AIC = ifelse(is.null(warning), NA, warning) ) )
}
 
