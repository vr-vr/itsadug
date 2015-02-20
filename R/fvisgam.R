#' Visualization of nonlinear interactions.
#'
#' @export
#' @aliases vis.gam2
#' @description Produces perspective or contour plot views of gam model 
#' predictions of the additive effects interactions.
#' The code is based on the script for \code{\link[mgcv]{vis.gam}}, 
#' but allows to cancel random effects.
#'
#' @param x A gam object, produced by \code{\link[mgcv]{gam}} or 
#' \code{\link[mgcv]{bam}}.
#' @param view A two-value vector containing the names of the two main effect 
#' terms to be displayed on the x and y dimensions of the plot. Note that 
#' variables coerced to factors in the model formula won't work as view 
#' variables.
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
#' @param plot.type one of "contour" or "persp" (default is "contour").
#' @param zlim A two item array giving the lower and upper limits for the z-
#' axis scale. NULL to choose automatically.
#' @param nCol The number of colors to use in color schemes.
#' @param rm.ranef Logical: whether or not to remove random effects. 
#' Default is TRUE.
#' @param suppressMessages Logical: whether or not to print messages.
#' @param ... other options to pass on to persp, image or contour. In 
#' particular ticktype="detailed" will add proper axes labeling to the plots.
#'
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
#' # Same plot:
#' fvisgam(m1, view=c("Time", "Trial"), rm.ranef=FALSE)
#' # Without random effects included:
#' fvisgam(m1, view=c("Time", "Trial"), rm.ranef=TRUE)
#' }
#' # see the vignette for examples:
#' vignette("plotfunctions", package="itsadug")
#' @author Jacolien van Rij and Martijn Wieling. 
#' Modification of \code{\link[mgcv]{vis.gam}} from 
#' package \code{\link[mgcv]{mgcv}} of Simon N. Wood.
#' @seealso \code{\link[mgcv]{vis.gam}}, \code{\link[mgcv]{plot.gam}} 
#'
#' @family functions for interpreting nonlinear effects

fvisgam <- function(x, view = NULL, cond = list(), 
    n.grid = 30, too.far = 0, col = NA, color = "topo", contour.col = NULL, 
    add.color.legend=TRUE,
    se = -1, plot.type = "contour", zlim = NULL, nCol = 50, 
    rm.ranef=NULL, suppressMessages=FALSE, ...) {
       
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
        stop("Specify two view predictors for the x- and y-axis.")
    } else {
        if (sum(view %in% v.names) != 2) {
            stop(paste(c("view variables must be one of", v.names), collapse = ", "))
        }
        for (i in 1:2) if (!inherits(x$var.summary[[view[i]]], c("numeric"))) 
            stop("Don't know what to do with parametric terms that are not simple numeric variables.")
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
    

    m1 <- seq(min(x$var.summary[[view[1]]], na.rm=TRUE), 
        max(x$var.summary[[view[1]]], na.rm=TRUE), length=n.grid)
    m2 <- seq(min(x$var.summary[[view[2]]], na.rm=TRUE), 
        max(x$var.summary[[view[2]]], na.rm=TRUE), length=n.grid)

    cond[[view[1]]] <- m1
    cond[[view[2]]] <- m2

    newd <- get_predictions(x, cond=cond, se=ifelse(se>0, TRUE, FALSE), 
        f=ifelse(se>0, se, 1.96), rm.ranef=rm.ranef,
        print.summary=!suppressMessages)
    newd <- newd[order(newd[,view[1]], newd[, view[2]]),]

    z <- matrix(newd$fit, byrow=TRUE, n.grid, n.grid)

    zlab <- colnames(x$model)[!colnames(x$model) %in% names(cond)][1]
   
    if (too.far > 0) {
        ex.tf <- mgcv::exclude.too.far(m1, m2, x$model[, view[1]], x$model[, view[2]], dist = too.far)
        newd$se.fit[ex.tf] <- newd$fit[ex.tf] <- NA
    }
    if (se <= 0) {
        old.warn <- options(warn = -1)
        av <- matrix(c(0.5, 0.5, rep(0, n.grid - 1)), byrow=TRUE, n.grid, n.grid - 1)
        options(old.warn)
        max.z <- max(z, na.rm = TRUE)
        z[is.na(z)] <- max.z * 10000
        z <- matrix(z, byrow=TRUE, n.grid, n.grid)
        surf.col <- t(av) %*% z %*% av
        surf.col[surf.col > max.z * 2] <- NA
        if (!is.null(zlim)) {
            if (length(zlim) != 2 || zlim[1] >= zlim[2]) 
                stop("Something wrong with zlim")
            min.z <- zlim[1]
            max.z <- zlim[2]
        } else {
            min.z <- min(newd$fit, na.rm = TRUE)
            max.z <- max(newd$fit, na.rm = TRUE)
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
        } else stop("color scheme not recognized")
        if (is.null(contour.col)) 
            contour.col <- con.col
        surf.col[surf.col < 1] <- 1
        surf.col[surf.col > nCol] <- nCol
        if (is.na(col)) 
            col <- pal[as.array(surf.col)]
        z <- matrix(newd$fit, byrow=TRUE, n.grid, n.grid)
        if (plot.type == "contour") {
            stub <- paste(ifelse("xlab" %in% dnm, "", ",xlab=view[1]"), ifelse("ylab" %in% dnm, "", ",ylab=view[2]"), ifelse("main" %in% 
                dnm, "", ",main=zlab"), ",...)", sep = "")
            if (color != "bw") {
                txt <- paste("image(m1,m2,z,col=pal,zlim=c(min.z,max.z)", stub, sep = "")
                eval(parse(text = txt))
                txt <- paste("contour(m1,m2,z,col=contour.col,zlim=c(min.z,max.z)", ifelse("add" %in% dnm, "", ",add=TRUE"), 
                  ",...)", sep = "")
                eval(parse(text = txt))
            } else {
                txt <- paste("contour(m1,m2,z,col=1,zlim=c(min.z,max.z)", stub, sep = "")
                eval(parse(text = txt))
            }
            if(add.color.legend){
                gradientLegend(round(c(min.z, max.z), 3), n.seg=3, pos=.875, 
                    color=pal)
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
            z.max <- max(newd$fit + newd$CI, na.rm = TRUE)
            z.min <- min(newd$fit - newd$CI, na.rm = TRUE)
        }
        zlim <- c(z.min, z.max)
        z <- newd$fit - newd$CI
        z <- matrix(z, byrow=TRUE, n.grid, n.grid)
        if (plot.type == "contour") 
            warning("sorry no option for contouring with errors: try plot.gam")
        stub <- paste(ifelse("xlab" %in% dnm, "", ",xlab=view[1]"), ifelse("ylab" %in% dnm, "", ",ylab=view[2]"), ifelse("zlab" %in% 
            dnm, "", ",zlab=zlab"), ifelse("sub" %in% dnm, "", ",sub=subs"), ",...)", sep = "")
        txt <- paste("persp(m1,m2,z,col=col,zlim=zlim", ifelse("border" %in% dnm, "", ",border=lo.col"), stub, sep = "")
        eval(parse(text = txt))
        par(new = TRUE)
        z <- newd$fit
        z <- matrix(z, byrow=TRUE, n.grid, n.grid)
        txt <- paste("persp(m1,m2,z,col=col,zlim=zlim", ifelse("border" %in% dnm, "", ",border=\"black\""), stub, sep = "")
        eval(parse(text = txt))
        par(new = TRUE)
        z <- newd$fit + newd$CI
        z <- matrix(z, byrow=TRUE, n.grid, n.grid)
        txt <- paste("persp(m1,m2,z,col=col,zlim=zlim", ifelse("border" %in% dnm, "", ",border=hi.col"), stub, sep = "")
        eval(parse(text = txt))
    }
    invisible(list(fv = newd, m1 = m1, m2 = m2))
}

 
