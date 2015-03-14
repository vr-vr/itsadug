#' Visualization of nonlinear interactions.
#'
#' @export
#' @aliases pvis.gam
#' @description Produces perspective or contour plot views of gam model 
#' predictions of the partial effects interactions. Combines the function 
#' \code{\link[mgcv]{plot.gam}} for interaction surfaces with the function 
#' \code{\link[mgcv]{vis.gam}}. Similar to \code{\link[mgcv]{plot.gam}}, 
#' \code{pvisgam} plots the partial interaction surface, without including 
#' values for other predictors that are not being shown. Similar to 
#' \code{\link[mgcv]{vis.gam}} the user can set the two predictors to be 
#' viewed, and colors are added behind the contours to facilitate 
#' interpretation. In contrast to \code{\link[mgcv]{plot.gam}}, this function 
#' allows to plotting of interactions with three of more continuous predictors 
#' by breaking it down in two-dimensional surfaces.
#' The code is derivated from the script for \code{\link[mgcv]{vis.gam}}.
#' @param x A gam object, produced by \code{\link[mgcv]{gam}} or 
#' \code{\link[mgcv]{bam}}.
#' @param view A two-value vector containing the names of the two main effect 
#' terms to be displayed on the x and y dimensions of the plot. Note that 
#' variables coerced to factors in the model formula won't work as view 
#' variables.
#' @param select  A number, selecting a single model term for printing. e.g. 
#' if you want the plot for the second smooth term set select=2.
#' @param cond A named list of the values to use for the other predictor terms 
#' (not in view). Used for choosing between smooths that share the same view 
#' predictors.
#' @param n.grid  The number of grid nodes in each direction used for 
#' calculating the plotted surface.
#' @param too.far Plot grid nodes that are too far from the points defined by 
#' the variables given in view can be excluded from the plot. too.far 
#' determines what is too far. The grid is scaled into the unit square along 
#' with the view variables and then grid nodes more than too.far from the 
#' predictor variables are excluded.
#' @param col The colors for the facets of the plot.
#' @param color The color scheme to use for plots. One of "topo", "heat", 
#' "cm", "terrain", "gray" or "bw". 
#' @param contour.col sets the color of contours when using plot.
#' @param add.color.legend Logical: whether or not to add a color legend. 
#' Default is TRUE. If FALSE (omitted), one could use the function
#' \code{\link{gradientLegend}} to add a legend manually at any position.
#' @param se If less than or equal to zero then only the predicted surface is 
#' plotted, but if greater than zero, then 3 surfaces are plotted, one at the 
#' predicted values minus se standard errors, one at the predicted values and 
#' one at the predicted values plus se standard errors.
#' @param type "link" to plot on linear predictor scale and "response" to plot 
#' on the response scale.
#' @param plot.type one of "contour" or "persp" (default is "contour").
#' @param zlim A two item array giving the lower and upper limits for the z-
#' axis scale. NULL to choose automatically.
#' @param nCol The number of colors to use in color schemes.
#' @param labcex Size of the contour labels.
#' @param print.summary Logical: whether or not to print summary.
#' @param ... other options to pass on to persp, image or contour. In 
#' particular ticktype="detailed" will add proper axes labeling to the plots.
#' @section Warnings:
#' In contrast to vis.gam, do not specify other predictors in \code{cond} that 
#' are not to be plotted.
#' @examples
#' data(simdat)
#' 
#' \dontrun{
#' # Model with random effect and interactions:
#' m1 <- bam(Y ~ te(Time, Trial)+s(Time, Subject, bs='fs', m=1),
#'     data=simdat)
#'
#' # Plot summed effects:
#' vis.gam(m1, view=c("Time", "Trial"), plot.type='contour', color='topo')
#' # Partial effect of interaction:
#' pvisgam(m1, view=c("Time", "Trial"), select=1)
#' # Same:
#' plot(m1, select=1, scheme=2)
#' plot(m1, select=1)
#' # Alternatives:
#' pvisgam(m1, view=c("Trial", "Time"), select=1)
#' pvisgam(m1, view=c("Trial", "Time"), select=1, zlim=c(-20,20))
#' }
#' # see the vignette for examples:
#' vignette("plotfunctions", package="itsadug")
#' @author Jacolien van Rij. Modification of \code{\link[mgcv]{vis.gam}} from 
#' package \code{\link[mgcv]{mgcv}} of Simon N. Wood.
#' @seealso \code{\link[mgcv]{vis.gam}}, \code{\link[mgcv]{plot.gam}} 
#'
#' @family functions for interpreting nonlinear effects

pvisgam <- function(x, view = NULL, select = NULL, cond = list(), n.grid = 30, 
    too.far = 0, col = NA, color = "topo", contour.col = NULL, 
    add.color.legend=TRUE,
    se = -1, type = "link", plot.type = "contour", zlim = NULL, 
    nCol = 50, labcex=.6, print.summary=TRUE,...) {
    
    # This modfication of vis.gam allows the user to specify one condition to plot as partial effect surface.  Use: 1)
    # view=c('Time','Trial') to specify which surface to plot, and 2) select=2 to select a specific smooth term (necessary
    # for distinguishing between several levels of a predictor, such as te(Time, Trial):Condition2 of the tensor te(Time,
    # Trial,by=Cond).  3) cond=list(X=5) can be used to select a specific value of a continuous predictor in a complex
    # interaction, e.g. to specify the value of X in te(Time,Trial,X, by=Cond).  Important: do not specify other predictors
    # in cond that are not to be plotted.
    
    fac.seq <- function(fac, n.grid) {
        fn <- length(levels(fac))
        gn <- n.grid
        if (fn > gn) 
            mf <- factor(levels(fac))[1:gn] else {
            ln <- floor(gn/fn)
            mf <- rep(levels(fac)[fn], gn)
            mf[1:(ln * fn)] <- rep(levels(fac), rep(ln, fn))
            mf <- factor(mf, levels = levels(fac))
        }
        mf
    }
    dnm <- names(list(...))
    v.names <- names(x$var.summary)
    if (is.null(view)) {
        k <- 0
        view <- rep("", 2)
        for (i in 1:length(v.names)) {
            ok <- TRUE
            if (is.matrix(x$var.summary[[i]])) 
                ok <- FALSE else if (is.factor(x$var.summary[[i]])) {
                if (length(levels(x$var.summary[[i]])) <= 1) 
                  ok <- FALSE
            } else {
                if (length(unique(x$var.summary[[i]])) == 1) 
                  ok <- FALSE
            }
            if (ok) {
                k <- k + 1
                view[k] <- v.names[i]
            }
            if (k == 2) 
                break
        }
        if (k < 2) 
            stop("Model does not seem to have enough terms to do anything useful")
    } else {
        if (sum(view %in% v.names) != 2) {
            stop(paste(c("view variables must be one of", v.names), collapse = ", "))
        }
        for (i in 1:2) if (!inherits(x$var.summary[[view[i]]], c("numeric", "factor"))) 
            stop("Don't know what to do with parametric terms that are not simple numeric or factor variables")
    }
    ok <- TRUE
    for (i in 1:2) if (is.factor(x$var.summary[[view[i]]])) {
        if (length(levels(x$var.summary[[view[i]]])) <= 1) 
            ok <- FALSE
    } else {
        if (length(unique(x$var.summary[[view[i]]])) <= 1) 
            ok <- FALSE
    }
    if (!ok) 
        stop(paste("View variables must contain more than one value. view = c(", view[1], ",", view[2], ").", sep = ""))
    if (is.factor(x$var.summary[[view[1]]])) {
        m1 <- fac.seq(x$var.summary[[view[1]]], n.grid)
    } else {
        r1 <- range(x$var.summary[[view[1]]])
        m1 <- seq(r1[1], r1[2], length = n.grid)
    }
    if (is.factor(x$var.summary[[view[2]]])) {
        m2 <- fac.seq(x$var.summary[[view[2]]], n.grid)
    } else {
        r2 <- range(x$var.summary[[view[2]]])
        m2 <- seq(r2[1], r2[2], length = n.grid)
    }
    v1 <- rep(m1, n.grid)
    v2 <- rep(m2, rep(n.grid, n.grid))
    newd <- data.frame(matrix(0, n.grid * n.grid, 0))
    
    # add factor to condition list
    if (is.numeric(select)) {
        if (x$smooth[[select]]$by != "NA") {
            level <- x$smooth[[select]]$by.level
            if (is.null(level)) {
                level <- 1
            }
            cond[[x$smooth[[select]]$by]] = level
        }
    }
    
    for (i in 1:length(x$var.summary)) {
        ma <- cond[[v.names[i]]]
        
        # if no value for this variable is specified in cond, then take mean
        if (is.null(ma)) {
            ma <- x$var.summary[[i]]
            if (is.numeric(ma)) 
                ma <- ma[2]
        }
        
        if (is.matrix(x$var.summary[[i]])) {
            newd[[i]] <- matrix(ma, n.grid * n.grid, ncol(x$var.summary[[i]]), byrow = TRUE)
        } else {
            newd[[i]] <- rep(ma, n.grid * n.grid)
        }
    }
    names(newd) <- v.names
    newd[[view[1]]] <- v1
    newd[[view[2]]] <- v2
    if (type == "link") {
        zlab <- paste("linear predictor")
    } else if (type == "response") {
        zlab <- type
    } else stop("type must be \"link\" or \"response\"")
    
    # -------------------------------------------- NEW:
    
    X1 <- mgcv::predict.gam(x, newdata = newd, type = "terms", se.fit = TRUE)
    
    fv <- NULL
    
    # determine select value
    n.linpred <- 0
    if (length(attr(x$pterms, "term.labels")) > 0) {
        n.linpred <- length(attr(x$pterms, "term.labels"))
    }
    
    if (is.numeric(select)) {
        fv <- data.frame(fit = X1$fit[, select + n.linpred])
        fv$se.fit <- X1$se.fit[, select + n.linpred]
        if(print.summary){
            print(paste("Tensor(s) to be plotted:", colnames(X1$fit)[select + n.linpred]))
        }
        
        
    } else {
        
        if (!is.na(view[1])) {
            
            colnamesX1 <- NA
            
            for (i in 1:length(view)) {
                if (is.na(colnamesX1[1])) 
                  colnamesX1 <- colnames(X1$fit)[grepl(view[i], colnames(X1$fit))] else colnamesX1 <- colnamesX1[colnamesX1 %in% colnames(X1$fit)[grepl(view[i], colnames(X1$fit))]]
            }
            
            if (length(colnamesX1) > 1) {
                if (length(cond) > 0) {
                  select = c()
                  for (i in 1:length(cond)) {
                    # check if cond is factor:
                    if (!is.numeric(x$var.summary[[names(cond[i])]])) {
                      test <- strsplit(colnamesX1, names(cond[i]))
                      for (j in 1:length(test)) {
                        if ((length(test[[j]]) > 1) & (grepl(test[[j]][2], as.character(cond[[i]])))) {
                          select = c(select, j)
                        }
                      }
                      colnamesX1 <- colnamesX1[select]
                    }
                  }
                }
            }
        }
        


        
        if (length(colnamesX1) == 1) {
            fv <- data.frame(fit = X1$fit[, colnamesX1])
            fv$se.fit <- X1$se.fit[, colnamesX1]
            if(print.summary){
                cat(paste("Tensor(s) to be plotted:", colnamesX1))
            }
        } else {
            stop(sprintf("More than one level is selected for plotting: %s. Use 'select=(number of smooth, found in summary). Note that the surface does not reflect a true average of their partial effects. Please consider specifying one of these conditions in the cond parameter."
                ,paste(colnamesX1, collapse = " + ")), immediate. = TRUE)
        }
    }
    
    # END NEW --------------------------------------------
    
    z <- fv$fit
    if (too.far > 0) {
        ex.tf <- mgcv::exclude.too.far(v1, v2, x$model[, view[1]], x$model[, view[2]], dist = too.far)
        fv$se.fit[ex.tf] <- fv$fit[ex.tf] <- NA
    }
    if (is.factor(m1)) {
        m1 <- as.numeric(m1)
        m1 <- seq(min(m1) - 0.5, max(m1) + 0.5, length = n.grid)
    }
    if (is.factor(m2)) {
        m2 <- as.numeric(m2)
        m2 <- seq(min(m1) - 0.5, max(m2) + 0.5, length = n.grid)
    }
    if (se <= 0) {
        old.warn <- options(warn = -1)
        av <- matrix(c(0.5, 0.5, rep(0, n.grid - 1)), n.grid, n.grid - 1)
        options(old.warn)
        max.z <- max(z, na.rm = TRUE)
        z[is.na(z)] <- max.z * 10000
        z <- matrix(z, n.grid, n.grid)
        surf.col <- t(av) %*% z %*% av
        surf.col[surf.col > max.z * 2] <- NA
        if (!is.null(zlim)) {
            if (length(zlim) != 2 || zlim[1] >= zlim[2]) 
                stop("Something wrong with zlim")
            min.z <- zlim[1]
            max.z <- zlim[2]
        } else {
            min.z <- min(fv$fit, na.rm = TRUE)
            max.z <- max(fv$fit, na.rm = TRUE)
        }
        surf.col <- surf.col - min.z
        surf.col <- surf.col/(max.z - min.z)
        surf.col <- round(surf.col * nCol)
        con.col <- 1
        if (color == "heat") {
            pal <- heat.colors(nCol)
            con.col <- 3
        } else if (color == "topo") {
            pal <- topo.colors(nCol)
            con.col <- 2
        } else if (color == "cm") {
            pal <- cm.colors(nCol)
            con.col <- 1
        } else if (color == "terrain") {
            pal <- terrain.colors(nCol)
            con.col <- 2
        } else if (color == "bpy") {
            if (requireNamespace("sp", quietly = TRUE)) {
                pal <- sp::bpy.colors(nCol)
                con.col <- 1
            } else {
                warning("Package 'sp' needed for bpy color palette. Using topo.colors instead (default).")
                color <- 'topo'
                pal <- topo.colors(nCol)
                con.col <- 2
            }
        } else if (color == "gray" || color == "bw") {
            pal <- gray(seq(0.1, 0.9, length = nCol))
            con.col <- 1
        } else stop("color scheme not recognised")
        if (is.null(contour.col)) 
            contour.col <- con.col
        surf.col[surf.col < 1] <- 1
        surf.col[surf.col > nCol] <- nCol
        if (is.na(col)) 
            col <- pal[as.array(surf.col)]
        z <- matrix(fv$fit, n.grid, n.grid)
        if (plot.type == "contour") {
            stub <- paste(ifelse("xlab" %in% dnm, "", ",xlab=view[1]"), ifelse("ylab" %in% dnm, "", ",ylab=view[2]"), ifelse("main" %in% 
                dnm, "", ",main=zlab"), ",...)", sep = "")
            if (color != "bw") {
                txt <- paste("image(m1,m2,z,col=pal,zlim=c(min.z,max.z)", stub, sep = "")
                eval(parse(text = txt))
                txt <- paste("contour(m1,m2,z,col=contour.col,zlim=c(min.z,max.z), labcex=labcex", ifelse("add" %in% dnm, "", ",add=TRUE"), 
                  ",...)", sep = "")
                eval(parse(text = txt))
            } else {
                txt <- paste("contour(m1,m2,z,col=1,zlim=c(min.z,max.z), labcex=labcex", stub, sep = "")
                eval(parse(text = txt))
            }
            if(add.color.legend){
                gradientLegend(round(c(min.z, max.z), 3), n.seg=3, pos=.875, color=pal)
            }
        }else{
             stub <- paste(ifelse("xlab" %in% dnm, "", ",xlab=view[1]"), ifelse("ylab" %in% dnm, "", ",ylab=view[2]"), ifelse("main" %in% 
                dnm, "", ",main=zlab"), ",...)", sep = "")
            if (color == "bw") {
                op <- par(bg = "white")
                txt <- paste("persp(m1,m2,z,col=\"white\",zlim=c(min.z,max.z) ", stub, sep = "")
                eval(parse(text = txt))
                par(op)
            }
            else {
                txt <- paste("persp(m1,m2,z,col=col,zlim=c(min.z,max.z)", 
                  stub, sep = "")
                eval(parse(text = txt))
            }           
        }
    } else {
        if (color == "bw" || color == "gray") {
            subs <- paste("grey are +/-", se, "s.e.")
            lo.col <- "gray"
            hi.col <- "gray"
        } else {
            subs <- paste("red/green are +/-", se, "s.e.")
            lo.col <- "green"
            hi.col <- "red"
        }
        if (!is.null(zlim)) {
            if (length(zlim) != 2 || zlim[1] >= zlim[2]) 
                stop("Something wrong with zlim")
            min.z <- zlim[1]
            max.z <- zlim[2]
        } else {
            z.max <- max(fv$fit + fv$se.fit * se, na.rm = TRUE)
            z.min <- min(fv$fit - fv$se.fit * se, na.rm = TRUE)
        }
        zlim <- c(z.min, z.max)
        z <- fv$fit - fv$se.fit * se
        z <- matrix(z, n.grid, n.grid)
        if (plot.type == "contour") 
            warning("sorry no option for contouring with errors: try plot.gam")
        stub <- paste(ifelse("xlab" %in% dnm, "", ",xlab=view[1]"), ifelse("ylab" %in% dnm, "", ",ylab=view[2]"), ifelse("zlab" %in% 
            dnm, "", ",zlab=zlab"), ifelse("sub" %in% dnm, "", ",sub=subs"), ",...)", sep = "")
        txt <- paste("persp(m1,m2,z,col=col,zlim=zlim", ifelse("border" %in% dnm, "", ",border=lo.col"), stub, sep = "")
        eval(parse(text = txt))
        par(new = TRUE)
        z <- fv$fit
        z <- matrix(z, n.grid, n.grid)
        txt <- paste("persp(m1,m2,z,col=col,zlim=zlim", ifelse("border" %in% dnm, "", ",border=\"black\""), stub, sep = "")
        eval(parse(text = txt))
        par(new = TRUE)
        z <- fv$fit + se * fv$se.fit
        z <- matrix(z, n.grid, n.grid)
        txt <- paste("persp(m1,m2,z,col=col,zlim=zlim", ifelse("border" %in% dnm, "", ",border=hi.col"), stub, sep = "")
        eval(parse(text = txt))
    }

    invisible(list(fv = fv, m1 = m1, m2 = m2))
}

 
