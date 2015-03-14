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
#' Default is FALSE. Alternatively a string (or vector of strings) with the 
#' name of the random effect(s) to remove.
#' @param print.summary Logical: whether or not to print a summary of the 
#' values selected for each predictor. Defaults to TRUE.
#' @return A data frame with estimates and optionally errors.
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
				mysummary[[pred]] <- paste(mysummary[[pred]], "(Might be canceled as random effect, check below.)")
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
#' @param comp A named list with the two levels to compare.
#' @param cond A named list of the values to use for the other predictor 
#' terms. Variables omitted from this list will have the closest observed 
#' value to the median for continuous variables, or the reference level for 
#' factors. 
#' @param rm.ranef Logical: whether or not to remove random effects. 
#' Default is FALSE. Alternatively a string (or vector of strings) with the 
#' name of the random effect(s) to remove. 
#' (See notes.)
#' @param se Logical: whether or not to return the confidence interval or 
#' standard error around the estimates.
#' @param f A number to scale the standard error. Defaults to 1.96, resulting 
#' in 95\% confidence intervals. For 99\% confidence intervals use a value of 
#' 2.58.
#' @param print.summary Logical: whether or not to print a summary of the 
#' values selected for each predictor. Defaults to TRUE.
#' @return Returns a data frame with the estimates of the difference and 
#' optionally the confidence intervals around that estimate.
#' @section Notes:
#' Other, not specified effects and random effects are generally canceled 
#' out, when calculating the difference. When the predictors that 
#' specify the conditions to compare are involved in other interactions
#' or included as random slopes, it may be useful to specify the values
#' of other predictors with \code{cond} or remove the random effects with
#' \code{rm.ranef}. 
#' @examples
#' data(simdat)
#' 
#' # first fit a simple model:
#' m1 <- bam(Y ~ Group+te(Time, Trial, by=Group), data=simdat)
#'
#' # get difference estimates:
#' diff <- get_difference(m1, comp=list(Group=c('Adults', 'Children')), 
#'     cond=list(Time=seq(0,500,length=100)))
#' head(diff)
#' @author Jacolien van Rij, Martijn Wieling
#' @family functions for model predictions

get_difference <- function(model, comp, cond=NULL,
	rm.ranef=NULL,
	se=TRUE, f=1.96, print.summary=TRUE){


	if(!"lm" %in% class(model)){
		stop("This function does not work for class %s models.", class(model)[1])
	}else{
		newd <- NULL
		su <- model$var.summary
		dat <- model$model

		# check comp
		if(is.null(names(comp))){
			stop("Predictor specified in 'comp' unknown. Please provide a named list for 'comp', in the form of 'comp=list(Predictor=c('level1', 'level2'))'.")
		}

		if(all(names(comp) %in% colnames(dat))){
			for(i in 1:length(comp)){
				if(length(comp[[i]]) < 2){
					stop(sprintf('Provide two levels for %s to calculate difference.', names(comp)[i]))
				}else if(length(comp[[i]]) > 2){
					warning(sprintf('More than two levels provided for predictor %s. Only first two levels are being used.', names(comp)[i]))
				}
			}
		}else{
			errname <- paste(which(!names(comp)%in% colnames(dat)), collapse=", ")
			stop(sprintf('Grouping predictor(s) not found in model: %s.', errname))
		}

		if(any(names(cond) %in% names(comp))){
			for(i in names(cond)[names(cond) %in% names(comp)] ){
				cond[[i]] <- NULL
				warning(sprintf('Predictor %s specified in comp and cond. (The value in cond will be ignored.)', i))
			}
		}

		new.cond1 <- list()
		new.cond2 <- list()

		for(i in names(su)){
			if(i %in% names(comp)){
				new.cond1[[i]] <- comp[[i]][1]
				new.cond2[[i]] <- comp[[i]][2]
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

		newd <- as.data.frame(newd1[,!names(newd1) %in% names(comp)])
		mysummary <- summary_data(newd, print=FALSE)	

		# Check for random effects:
		if(class(rm.ranef)=="logical"){
			if(rm.ranef==FALSE){
				rm.ranef <- NULL			
			}
		}
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
				}
			}
			rm.col <- unlist(lapply(rm.ranef, 
				function(x){
					colnames(p1)[grepl(x, colnames(p1), fixed=TRUE)]
				}))
			rm.col <- unlist(lapply(smoothlabels,
				function(x){
					rm.col[grepl(x, rm.col, fixed=TRUE)]
				}))

			# cancel random effects
			p1[,rm.col] <- 0
			p2[,rm.col] <- 0

			# find terms that only occur in random effects:
			predictors <- do.call('rbind',
				lapply(model$smooth, 
					function(x){
						data.frame(Label=x[['label']],
							Terms=x[['term']])
					} )) 	
			test <- table(predictors$Terms) - table(predictors[predictors$Label %in% rm.ranef,]$Terms)
			for(pred in names(test[test==0])){
				if(pred %in% names(mysummary)){
					mysummary[[pred]] <- paste(mysummary[[pred]], "(Might be canceled as random effect, check below.)")
				}
			}

			if(length(rm.col)>0){
				mysummary[['NOTE']] =  sprintf("The following random effects columns are canceled: %s\n", 
				paste(smoothlabels, collapse=","))
			}else{
				warning("No random effects to cancel.\n")				
			}

		}



		# calculate the difference:
		p <- p1 - p2

		newd$difference <- as.vector(p %*% coef(model))
		if(se){
			newd$CI <- f*sqrt(rowSums((p%*%vcov(model))*p))
		}

		# print summary of chosen values
		if(print.summary==TRUE){
			print_summary(mysummary)
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
#' @param as.data.frame Logical: whether or not to return a data.frame. 
#' Default is false, and a list will be returned.
#' @return A list with two or more elements:
#' \itemize{
#' \item \code{fit}: Numeric vector with the fitted values;
#' \item \code{se.fit}: Optionally, only with \code{se=TRUE}.
#' Numeric vector with the error or confidence interval values (f*SE);
#' \item \code{f}: The multiplication factor for generating 
#' the confidence interval values;
#' \item \code{terms}: Numeric vector (for 1-dimensional smooth) 
#' or data frame (more 2- or more dimensional surfaces) 
#' with values of the modelterms.
#' \item \code{title}: String with name of the model term.
#' \item \code{xlab}, \code{ylab}, or \code{labels}:
#' Labels for x-axis and optionally y-axis. Precise structure depends 
#' on type of smooth term: for 1-dimensional smooth only x-label is provided, 
#' for 2-dimensional smooths x-label and y-label are provided, 
#' for more complex smooths a vector of of labels is provided.
#' }
#' @examples
#' data(simdat)
#'
#'\dontrun{ 
#' # Model with random effect and interactions:
#' m1 <- bam(Y ~ s(Time) + s(Trial) 
#' +ti(Time, Trial)
#' +s(Time, Subject, bs='fs', m=1),
#' data=simdat)
#' 
#' # Get list with predictions:
#' p <- get_modelterm(m1, select=1)
#' emptyPlot(range(p$terms), range(p$fit), h=0)
#' plot_error(p$terms, p$fit, p$se.fit, shade=TRUE, xpd=TRUE)
#' 
#' # Plot random effects in separate panels:
#' pp <- get_modelterm(m1, select=4, as.data.frame=TRUE)
#' require(lattice)
#' lattice::xyplot(fit~Time|Subject, 
#'     data=pp, type="l",
#'     xlab="Time", ylab="Partial effect")
#' 
#' # Plot selection of subjects:
#' pp <- get_modelterm(m1, select=4, 
#'     cond=list(Subject=c('a01', 'a03', 'c16')),
#'     as.data.frame=TRUE)
#' lattice::xyplot(fit~Time|Subject, 
#'     data=pp, type="l",
#'     xlab="Time", ylab="Partial effect")
#'
#' # Or using the package ggplot2:
#' require(ggplot2)
#' pp <- get_modelterm(m1, select=4, as.data.frame=TRUE)
#' pg <- ggplot2::qplot(Time, fit, data = pp, 
#'     geom = c("point", "line"), colour = Subject)
#' pg + ggplot2::guides(col = guide_legend(nrow = 18))
#' }
#' 
#' # see the vignette for plot examples:
#' vignette("plotfunctions", package="itsadug")
#' @author Jacolien van Rij
#' @family functions for model predictions

get_modelterm <- function(model, select, cond=NULL, n.grid=30, 
	se=TRUE, f=1.96, as.data.frame=FALSE){

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

		if(as.data.frame){
			fv <- c()
			if(length(smoothterms) > 1){
				fv <- newd[,smoothterms]
			}else{
				fv <- data.frame(st=newd[,smoothterms])
				names(fv) <- smoothterms
			}
			fv$fit <-  p$fit[, smoothlabels]
			if(se){
				fv$se.fit <- f*p$se.fit[,smoothlabels]
			}
			return(fv)
		}else{
			fv <- list()
			
			if(se){
				fv[['fit']] <- p$fit[, smoothlabels]
				fv[['se.fit']] <- f*p$se.fit[,smoothlabels]
				fv[['f']] <- f
			}else{
				fv[['fit']] <- p[, smoothlabels]
			}

			fv[['terms']] <- newd[,smoothterms]
			fv[['title']] <- model$smooth[[select]]['label']
			if(length(smoothterms)==1){
				fv[['xlab']] <- smoothterms[1]
			}else if(length(smoothterms)==2){
				fv[['xlab']] <- smoothterms[1]
				fv[['ylab']] <- smoothterms[2]
			}else{
				fv[['labels']] <- smoothterms
			}	

			return(fv)
		}
	}
}


