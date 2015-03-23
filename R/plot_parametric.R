#' Visualization of group estimates.
#'
#' @export
#' @description Plots a smooth from a \code{\link[mgcv]{gam}} or 
#' \code{\link[mgcv]{bam}} model based on predictions.
#' In contrast with the default \code{\link[mgcv]{plot.gam}}, this function 
#' plots the summed effects and optionally removes the random effects.
#'
#' @param x A gam object, produced by \code{\link[mgcv]{gam}} or 
#' \code{\link[mgcv]{bam}}.
#' @param pred A named list of the values to use for the predictor terms 
#' to plot. 
#' @param cond A named list of the values to use for the other predictor terms 
#' (not in view). Used for choosing between smooths that share the same view 
#' predictors.
#' @param parametricOnly Logical: whether or not to cancel out all smooth 
#' terms and only use the predictors in the parametric summary. 
#' @param rm.ranef Logical: whether or not to remove random effects. 
#' Default is TRUE.
#' @param col The colors for the lines and the error bars of the plot.
#' @param se If less than or equal to zero then only the predicted surface is 
#' plotted, but if greater than zero, then the predicted values plus 
#' confidence intervals are plotted. The value of se will be multiplied with 
#' the standard error (i.e., 1.96 results in 95\%CI and 2.58).
#' @param print.summary Logical: whether or not to print summary.
#' Default set to the print info messages option 
#' (see \code{\link{infoMessages}}).
#' @param main Changing the main title for the plot, see also title.
#' @param xlab Changing the label for the x axis, 
#' defaults to a description of x.
#' @param ... other options to pass on to \code{\link{dotplot_error}}, 
#' see \code{\link[graphics]{par}}
#' @section Warning:
#' Use \code{parametricOnly} with care! When set to TRUE, all smooth 
#' predictors are set to 0. Note that this might result in strange 
#' predictions, because a value of 0 does not always represents a realistic 
#' situation (e.g., body temperature of 0 is highly unlikely).  
#' Note that linear slopes are not set to zero, because they are 
#' considered as parametric terms. If \code{cond} does not specify a value for 
#' these continuous predictors, the closes value to the mean is automatically  
#' selected.
#'
#' @examples
#' data(simdat)
#' \dontrun{
#' m1 <- bam(Y ~ Group + te(Time, Trial, by=Group)
#'     + s(Time, Subject, bs='fs', m=1), data=simdat)
#' plot_parametric(m1, pred=list(Group=c('Adults', 'Children')))
#' # Note the summary that is printed.
#' 
#' # use rm.ranef to cancel random effects:
#' plot_parametric(m1, pred=list(Group=c('Adults', 'Children')),
#'     rm.ranef = TRUE)
#' 
#' # It is possible to get estimates that do not make sense:
#' out <- plot_parametric(m1, 
#'     pred=list(Group=c('Adults', 'Children'), Subject=c('a01', 'a02', 'c01')))
#' print(out)
#' }
#' 
#' # see the vignette for examples:
#' vignette("plotfunctions", package="itsadug")
#' @author Jacolien van Rij, based on a function of Fabian Tomaschek 
#' @seealso \code{\link[mgcv]{plot.gam}}, \code{\link{plot_diff}} 
#'
#' @family functions for interpreting nonlinear effects

plot_parametric <- function(x, pred, cond = list(), 
    parametricOnly = FALSE, rm.ranef=NULL, 
    col = 'black', se = 1.96, print.summary=getOption('itsadug_print'),
    main=NULL, xlab=NULL, ...) {
       
    dnm <- names(list(...))

    parTerms <- NULL
    if(parametricOnly){
        parTerms <- summary(x)$p.t
    }
    v.names <- names(x$var.summary)

    if (sum(names(pred) %in% v.names) != length(pred)) {
        stop(paste(c("Pred variable must be one of", v.names), collapse = ", "))
    }
    for(i in 1:length(names(pred))){
        if (!inherits(x$var.summary[[names(pred)[i]]], c("factor"))){
            stop("Don't know what to do with parametric terms that are not simple grouping variables.")
        }
    }

    if(!is.null(cond)){
        cn <- names(cond)
        test <- sapply(cn, function(x){
            if(length(unique(cond[[x]]))>1){
                stop("Do not specify more than 1 value for conditions listed in the argument cond.")
            }else{
                TRUE
            }
        })
    }

    for(i in names(pred)){
        cond[[i]] <- pred[[i]]
    }

    newd <- NULL
    if(parametricOnly){
        su <- x$var.summary
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


        p <- mgcv::predict.gam(x, newd, type='lpmatrix')
        rm.col <- colnames(p)[!colnames(p) %in% names(parTerms)]
        p[,rm.col] <- 0

        if(length(rm.col)==0){
            warning("No smooth terms in the model.\n")               
        }
        newd$fit <- p %*% coef(x)
        if(se>0){
            newd$CI <- se*sqrt(rowSums((p%*%vcov(x))*p))
        }

    }else{
        newd <- get_predictions(x, cond=cond, se=ifelse(se>0, TRUE, FALSE), 
            f=ifelse(se>0, se, 1.96), rm.ranef=rm.ranef,
            print.summary=print.summary)

    }

    newd$VnewCol <- NA
    newd <- droplevels(newd)
    if(length(pred)>1){
        newd$VnewCol <- interaction(newd[, names(pred)])
    }else{
        newd$VnewCol <- newd[,names(pred)[1]]
    }
    newd <- newd[order(newd$VnewCol),]


    if(is.null(main)){ main <- paste(names(pred), collapse=' x ') }
    if(is.null(xlab)){ xlab <- names(x$model)[!names(x$model) %in% v.names]}
    

    dotplot_error(x=as.vector(newd$fit), se.val=as.vector(newd$CI),
        labels=as.character(newd$VnewCol), 
        main=main, xlab=xlab, ...)
    abline(v=0, lty=3)

    newd$VnewCol <- NULL
    invisible(list(fv = newd))


}

 
