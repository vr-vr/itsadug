#' Get model predictions for specific conditions.
#' 
#' @export
#' @param model A gam object, produced by \code{\link[mgcv]{gam}} or 
#' \code{\link[mgcv]{bam}}.
#' @param cond A named list of the values to use for the predictor terms. 
#' Variables omitted from this list will have the closest observed value to 
#' the median for continuous variables, or the reference level for factors. 
#' @param se Logical: whether or not to return the confidence interval or 
#' standard error around the estimates.
#' @param f A number to scale the standard error. Defaults to 1.96, resulting 
#' in 95\% confidence intervals. For 99\% confidence intervals use a value of 
#' 2.58.
#' @param rm.ranef Logical: whether or not to remove random effects. 
#' Default is TRUE. Alternatively a string (or vector of strings) with the 
#' name of the random effect(s) to remove.
#' @param print.summary Logical: whether or not to print a summary of the 
#' values selected for each predictor. Defaults to TRUE.
#' @return A data frame with estimates and optinally errors.
#' @examples
#' # see the vignette for examples:
#' vignette("plotfunctions", package="itsadug")
#' @author Jacolien van Rij
#' @family functions for model predictions

get_predictions <- function(model, cond=NULL, se=TRUE, f=1.96, rm.ranef=NULL, 
	print.summary=TRUE){

	if(is.null(cond)){
		stop("Please specify values for at least one predictor in the parameter 'cond'.")
	}

	if(!"gam" %in% class(model)){
		stop("This function does not work for class %s models.", class(model)[1])
	}else{

		if(!is.null(rm.ranef)){
			if(rm.ranef==FALSE){
				rm.ranef <- NULL
			}
		}

		newd <- NULL
		su <- model$var.summary
		new.cond <- list()

		for(i in names(su)){
			if(i %in% names(cond)){
				new.cond[[i]] <- cond[[i]]
			}else{
				if(class(su[[i]])=="factor"){
					new.cond[[i]] <- as.character(su[[i]][1])
				}else if(class(su[[i]])=="numeric"){
					new.cond[[i]] <- su[[i]][2]

				}
			}
		}

		newd <- expand.grid(new.cond)
		mysummary <- summary_data(newd, print=FALSE)		

		p <- mgcv::predict.gam(model, newd, type='lpmatrix')

		if(!is.null(rm.ranef)){	

			# get random effects columns:
			smoothlabels <- as.data.frame( do.call('rbind', 
				lapply(model$smooth, 
					function(x){
						data.frame(Label=x[['label']], 
							Dim=x[['null.space.dim']], 
							Class = attr(x, "class")[1],
							stringsAsFactors=FALSE)
					} ) ) )
			# smoothlabels <- smoothlabels[smoothlabels$Dim==0,c("Label", "Class")]
			smoothlabels <- as.vector( smoothlabels[smoothlabels$Class %in% c("random.effect","fs.interaction"), "Label"] )

			if(class(rm.ranef)=="logical"){
				if(rm.ranef==TRUE){
					rm.ranef <- smoothlabels			
				}else if(rm.ranef==FALSE){
					rm.ranef <- ""
				}
			}
			rm.col <- unlist(lapply(rm.ranef, 
				function(x){
					colnames(p)[grepl(x, colnames(p), fixed=TRUE)]
				}))
			rm.col <- unlist(lapply(smoothlabels,
				function(x){
					rm.col[grepl(x, rm.col, fixed=TRUE)]
				}))

			# cancel random effects
			p[,rm.col] <- 0

			# find terms that only occur in random effects:
			predictors <- do.call('rbind',
				lapply(model$smooth, 
					function(x){
						data.frame(Label=x[['label']],
							Terms=x[['term']])
					} )) 	
			test <- table(predictors$Terms) - table(predictors[predictors$Label %in% rm.ranef,]$Terms)
			for(pred in names(test[test==0])){
				mysummary[[pred]] <- paste(mysummary[[pred]], "(Canceled as random effect.)")
			}

			if(length(rm.col)>0){
				mysummary[['NOTE']] =  sprintf("The following random effects columns are canceled: %s\n", 
				paste(smoothlabels, collapse=","))
			}else{
				warning("No random effects to cancel.\n")				
			}

		}

		if(print.summary){
			print_summary(mysummary)
		}
		
		newd$fit <- p %*% coef(model)
		if(se){
			newd$CI <- f*sqrt(rowSums((p%*%vcov(model))*p))
		}
		if(!is.null(rm.ranef)){
			newd$rm.ranef <- paste(smoothlabels, collapse=",")
		}
		
		return(newd)

	}
}


#' Get model predictions for the random effects.
#' 
#' @export
#' @param model A gam object, produced by \code{\link[mgcv]{gam}} or 
#' \code{\link[mgcv]{bam}}.
#' @param fun A string or function description to apply to the random effects 
#' estimates. When NULL (default), the estimates for the random effects are 
#' returned. 
#' @param cond A named list of the values to restrict the estimates for the 
#' random predictor terms. When NULL (default) all levels are returned.
#' @param n.grid Number of data points estimated for each random smooth.
#' @param print.summary Logical: whether or not to print a summary of the 
#' values selected for each predictor. Defaults to TRUE.
#' @return A data frame with estimates for random effects.
#' @examples
#' # see the vignette for examples:
#' vignette("plotfunctions", package="itsadug")
#' @author Jacolien van Rij
#' @family functions for model predictions

get_random <- function(model, fun=NULL, cond=NULL, n.grid=30, print.summary=TRUE){

	if(!"lm" %in% class(model)){
		stop("This function does not work for class %s models.", class(model)[1])
	}else{

		# find random effects:
		smoothlabels <- as.data.frame( do.call('rbind', 
			lapply(model$smooth, 
				function(x){
					data.frame(Label=x[['label']], 
						Dim=x[['null.space.dim']], 
						Class = attr(x, "class")[1],
						stringsAsFactors=FALSE)
				} ) ) )
		# smoothlabels <- smoothlabels[smoothlabels$Dim==0,c("Label", "Class")]
		smoothlabels <- as.vector( smoothlabels[smoothlabels$Class %in% c("random.effect","fs.interaction"), "Label"] )



		if(length(smoothlabels) == 0){
			warning("No random effects specified in the model.")
			return(NULL)
		}

		su <- model$var.summary
		newd <- NULL
		new.cond <- list()

		for(i in names(su)){
			if(i %in% names(cond)){
				new.cond[[i]] <- cond[[i]]
			}else{
				if(any(grepl(i, smoothlabels))){
					if(inherits(su[[i]],"factor")){
						new.cond[[i]] <- levels(su[[i]])
					}else if(inherits(su[[i]],"numeric")){
						new.cond[[i]] <- seq(su[[i]][1], su[[i]][3], length=n.grid)
					}
				}else{
					if(inherits(su[[i]],"factor")){
						new.cond[[i]] <- as.character(su[[i]][1])
					}else if(inherits(su[[i]],"numeric")){
						new.cond[[i]] <- su[[i]][2]
					}
				}
			}
		}

		newd <- expand.grid(new.cond)

		p <- mgcv::predict.gam(model, newd, type='terms')

		select <- p[,smoothlabels]
		cnames <- smoothlabels
		if(is.null(dim(select))){
			if(length(select) != length(smoothlabels)){
				# if length smoothlabels equals 1			
				newd <- cbind(select, newd)
				colnames(newd)[1] <- paste("fit",smoothlabels, sep='.')
				cnames <- paste("fit",smoothlabels, sep='.')
			}else{
				names(select) <- paste("fit",smoothlabels, sep='.')
				cnames <- names(select)
				newd <- cbind(t(select), newd)
			}
		}else{
			colnames(select) <- paste("fit",smoothlabels, sep='.')
			cnames <- colnames(select)
			newd <- cbind(select, newd)
		}
		

		if(!is.null(fun)){


			val <- NULL
			if(length(cnames) > 1){
				val <- as.list(newd[, cnames])
			}else{
				val <- list(newd[, cnames])
			}
			
			names(val) <- smoothlabels
			el <- sapply(model$smooth, function(x){return(x[['label']])})

			fun.val = list()
			for(i in names(val)){

				fun.cond <- list()

				for(j in model$smooth[[which(el==i)]][['term']]){
					if(!inherits(su[[j]],"factor")){
						fun.cond[[j]] <- newd[,j]
					}
				}			

				if(length(fun.cond) > 0){
					fun.val[[i]] <- aggregate(list(x=val[[i]]), by=fun.cond, fun)
				} else {
					fun.val[[i]] <- unlist(lapply(list(val[[i]]), fun))

				}

			}
			fun.val[['function']] <- fun

			return(fun.val)
		}else{
			# print summary of chosen values
			if(print.summary){
				summary_data(newd)
			}
			return(newd)
		}
	}
}


#' Get model predictions for differences between conditions.
#' 
#' @export
#' @param model A gam object, produced by \code{\link[mgcv]{gam}} or 
#' \code{\link[mgcv]{bam}}.
#' @param cond1 A named list with the values to use for the first difference 
#' condition.
#' @param cond2 A named list with the values to use for the second difference 
#' condition.
#' @param cond A named list of the values to use for the other predictor 
#' terms. Variables omitted from this list will have the closest observed 
#' value to the median for continuous variables, or the reference level for 
#' factors. 
#' @param se Logical: whether or not to return the confidence interval or 
#' standard error around the estimates.
#' @param f A number to scale the standard error. Defaults to 1.96, resulting 
#' in 95\% confidence intervals. For 99\% confidence intervals use a value of 
#' 2.58.
#' @param print.summary Logical: whether or not to print a summary of the 
#' values selected for each predictor. Defaults to TRUE.
#' @return Returns a data frame with the estimates of the difference and 
#' optionally the confidence intervals around that estimate.
#' @examples
#' data(simdat)
#' 
#' # first fit a simple model:
#' m1 <- bam(Y ~ Group+te(Time, Trial, by=Group), data=simdat)
#'
#' # get difference estimates:
#' diff <- get_difference(m1, cond1=list(Group='Adults'), 
#' cond2=list(Group='Children'), cond=list(Time=seq(0,500,length=100)))
#' head(diff)
#' @author Jacolien van Rij, Martijn Wieling
#' @family functions for model predictions

get_difference <- function(model, cond1, cond2, cond=NULL,
	se=TRUE, f=1.96, print.summary=TRUE){


	if(is.null(cond)){
		stop("Please specify values for at least one predictor in the parameter 'cond'.")
	}

	if(!"lm" %in% class(model)){
		stop("This function does not work for class %s models.", class(model)[1])
	}else{
		newd <- NULL
		su <- model$var.summary

		if(names(cond1) != names(cond2)){
			stop("The arguments cond1 and cond2 do not specify different levels of the same predictors.")
		}
		if(any(names(cond) %in% names(cond1))){
			for(i in names(cond)[names(cond) %in% names(cond1)] ){
				cond[[i]] <- NULL
			}
			warning(sprintf('Predictor %s specified in cond1, cond2 and cond.', i))
		}

		new.cond1 <- list()
		new.cond2 <- list()

		for(i in names(su)){
			if(i %in% names(cond1)){
				new.cond1[[i]] <- cond1[[i]]
				new.cond2[[i]] <- cond2[[i]]
			}else if(i %in% names(cond)){
				new.cond1[[i]] <- new.cond2[[i]] <- cond[[i]]
			}else{
				if(class(su[[i]])=="factor"){
					new.cond1[[i]] <- as.character(su[[i]][1])
					new.cond2[[i]] <- as.character(su[[i]][1])
				}else if(class(su[[i]])=="numeric"){
					new.cond1[[i]] <- su[[i]][2]
					new.cond2[[i]] <- su[[i]][2]
				}
			}
		}

		newd1 <- expand.grid(new.cond1)
		newd2 <- expand.grid(new.cond2)

		p1 <- mgcv::predict.gam(model, newd1, type='lpmatrix') 
		p2 <- mgcv::predict.gam(model, newd2, type='lpmatrix')
		p <- p1 - p2

		newd <- as.data.frame(newd1[,names(cond)])
		names(newd) <- names(cond)
		newd1 <- as.data.frame(newd1[,names(cond1)])
		names(newd1) <- paste(names(cond1), '1', sep='..')
		newd2 <- as.data.frame(newd2[,names(cond2)])
		names(newd2) <- paste(names(cond2), '2', sep='..')
		newd <- cbind(newd, newd1, newd2)

		newd$difference <- as.vector(p %*% coef(model))
		if(se){
			newd$CI <- f*sqrt(rowSums((p%*%vcov(model))*p))
		}

		# print summary of chosen values
		if(print.summary){
			summary_data(newd)
		}
		
		return(newd)

	}
}



#' Get estimated for selected model terms.
#' 
#' @export
#' @param model A gam object, produced by \code{\link[mgcv]{gam}} or 
#' \code{\link[mgcv]{bam}}.
#' @param select A number, indicating the model term to be selected. 
#' @param cond A named list of the values to restrict the estimates for the 
#' random predictor terms. When NULL (default) all levels are returned.
#' Only relevant for complex interactions, which involve more than two 
#' dimensions.
#' @param n.grid Number of data points estimated for each random smooth.
#' @return A data frame with estimates for the selected smooth term.
#' @param se Logical: whether or not to return the confidence interval or 
#' standard error around the estimates.
#' @param f A number to scale the standard error. Defaults to 1.96, resulting 
#' in 95\% confidence intervals. For 99\% confidence intervals use a value of 
#' 2.58.
#' @examples
#' # see the vignette for examples:
#' vignette("plotfunctions", package="itsadug")
#' @author Jacolien van Rij
#' @family functions for model predictions

get_modelterm <- function(model, select, cond=NULL, n.grid=30, 
	se=TRUE, f=1.96){

	if(!"lm" %in% class(model)){
		stop("This function does not work for class %s models.", class(model)[1])
	}else{

		# find terms:
		smoothlabels <- model$smooth[[select[1]]][['label']]
		smoothterms <-  model$smooth[[select[1]]][['term']]
		# select right grouping predictor:
		if(model$smooth[[select[1]]][['by']] !="NA"){
			cond[[model$smooth[[select[1]]][['by']]]] <- model$smooth[[select[1]]][['by.level']]
		}

		su <- model$var.summary
		newd <- NULL
		new.cond <- list()

		for(i in names(su)){
			if((i %in% names(cond)) & any(grepl(i, smoothlabels)) ){
				new.cond[[i]] <- cond[[i]]
			}else{
				if(any(grepl(i, smoothlabels))){
					if(inherits(su[[i]],"factor")){
						new.cond[[i]] <- levels(su[[i]])
					}else if(inherits(su[[i]],"numeric")){
						new.cond[[i]] <- seq(su[[i]][1], su[[i]][3], length=n.grid)
					}
				}else{
					if(inherits(su[[i]],"factor")){
						new.cond[[i]] <- as.character(su[[i]][1])
					}else if(inherits(su[[i]],"numeric")){
						new.cond[[i]] <- su[[i]][2]
					}
				}
			}
		}

		newd <- expand.grid(new.cond)
		p <- mgcv::predict.gam(model, newd, type='terms', se.fit=se)

		fv <- list()
		
		if(se){
			fv[['fit']] <- p$fit[, smoothlabels]
			fv[['se.fit']] <- f*p$se.fit[,smoothlabels]
			fv[['f']] <- f
		}else{
			fv[['fit']] <- p[, smoothlabels]
		}
		fv[['terms']] <- newd[,smoothterms]
		return(fv)
	}
}


