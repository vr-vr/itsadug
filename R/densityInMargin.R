#' Plot function
#' 
#' @description Plot density in margins of the plot.
#' 
#' @export
#' @import grDevices
#' @import graphics
#' @param x Density object.
#' @param side Number: 1 = bottom, 2 = left, 3 = top, 4 = left
#' @param from A number indicating the starting position (bottom) of the 
#' density plot. Defaults to 0, which is the border of the plot. 
#' Measured in proportions of the margin area available. Note that  
#' value could be negative (starting in the plot region).
#' @param maxDensityValue Number for scaling the density axis. 
#' Default is NULL (automatic scaling fitting the d)
#' @param otherDensities List with other density objects to determine 
#' the plotting scale such that they all fit. Defaults to NULL.
#' @param plot Logical: whether to plot the density (default) or not.
#' @param ... Optional arguments for the lines and fill_area. See \code{\link{par}}.
#' @author Jacolien van Rij
#' @examples
#' # density of a random sample from normal distribution:
#' val1 <- qnorm(ppoints(500))
#' val2 <- qt(ppoints(500), df = 2)
#' dens1 <- density(val1)
#' dens2 <- density(val2)
#' 
#' # setup plot window:
#' par(mfrow=c(1,1), cex=1.1)
#' 
#' # increase margin
#' oldmar <- par()$mar 
#' par(mar=oldmar + c(0,0,0,4))
#' 
#' # plot qqnorm
#' qqnorm(val2, main='t distribution',
#'        pch="*", col='steelblue',
#'        xlim=c(-3,3),
#'        bty='n')
#' qqline(val1)
#' abline(h=0, col=alpha('gray'))
#' abline(v=0, col=alpha('gray'))
#' 
#' # filled distribution in right margin:
#' marginDensityPlot(dens2, side=4, allDensities=list(dens1, dens2),
#' 	col='steelblue',lwd=2)
#' # add lines:
#' marginDensityPlot(dens2, side=4, allDensities=list(dens1, dens2),
#' 	col='steelblue',density=25, lwd=2)
#' # compare to normal:
#' marginDensityPlot(dens1, side=4, allDensities=list(dens1, dens2), 
#' 	col=NA, border=1)
#' # Other sides are also possible:
#' marginDensityPlot(dens1, side=3, allDensities=list(dens1, dens2), 
#' 	col=NA, border=alpha(1), lwd=2)
#' marginDensityPlot(dens2, side=3, allDensities=list(dens1, dens2), 
#' 	col=NA, border=alpha('steelblue'), lwd=3)
#' # adjust the starting point with argument from:
#' marginDensityPlot(dens1, side=1, 
#' 	from=.5, lwd=2)
#' marginDensityPlot(dens2, side=1, 
#' 	col='steelblue', from=-.90, lwd=2,
#'  maxDensityValue=2*max(dens2$y))
#' 
#' legend(getFigCoords('p')[2], getFigCoords('p')[3],
#' 	yjust=0,
#' 	legend=c("t distribution", "Gaussian"),
#' 	fill=c("steelblue", 'black'),
#' 	cex=.75,
#' 	xpd=TRUE, bty='n')
#' 
#' 
#' @seealso \code{\link{check_normaldist}}
#' @family Utility functions for plotting
marginDensityPlot <- function(x, y=NULL, side, from=0, 
	maxDensityValue=NULL, 
	allDensities=NULL, plot=TRUE, ...){

    if(!inherits(x, "density")){
        if(is.null(y)){
            d <- density(x, na.rm=TRUE)
            x <- d$x
            y <- d$y
            message("x converted to density object.")

        }else{
            if(length(x) != length(y)){
                stop("x and y do not have the same length.")
            }
        }
    }else{
        y <- x$y
        x <- x$x
    }
    if(is.null(maxDensityValue) & is.null(allDensities)){
        maxDensityValue = max(y, na.rm=TRUE)
    }else if (is.null(maxDensityValue) & !is.null(allDensities)){
    	maxDensityValue <- max( unlist( lapply(allDensities, function(a){ max(a$y)}) ) )
    }

    # set region:
    x0 <- y0 <- 0
    x1 <- y1 <- 1
    y.range <- 1
    y.dist <- 1

    gfc.f <- getFigCoords("f")
    gfc.p <- getFigCoords("p")

    horiz = TRUE

    if( side==1){       # bottom, going down
        x0 <- gfc.p[1]
        x1 <- gfc.p[2]
        y.range <- gfc.f[3] - gfc.p[3]
        y0 <- gfc.p[3] + from*y.range
        y1 <- gfc.f[3] - .05*y.range
        y.dist <- y1-y0
    }else if (side==2){   # left
        x0 <- gfc.p[3]
        x1 <- gfc.p[4]
        y.range <- gfc.f[1] - gfc.p[1]
        y0 <- gfc.p[1] + from*y.range
        y1 <- gfc.f[1] - .05*y.range
        y.dist <- y1-y0
        horiz = FALSE
    }else if (side==3){   # top
        x0 <- gfc.p[1]
        x1 <- gfc.p[2]
        y.range <- gfc.f[4] - gfc.p[4]
        y0 <- gfc.p[4] + from*y.range
        y1 <- gfc.f[4] - 0.05*y.range
        y.dist <- y1-y0
    }else if (side==4){   # right
        x0 <- gfc.p[3]
        x1 <- gfc.p[4]
        y.range <- gfc.f[2] - gfc.p[2]
        y0 <- gfc.p[2] + from*y.range
        y1 <- gfc.f[2] - 0.05*y.range
        y.dist <- y1-y0
        horiz = FALSE
    }

    scale <- y.dist / maxDensityValue

    if(plot){
        fill_area(x, y*scale+y0, from=y0, horiz=horiz, xpd=TRUE, ...)
    }
    
    invisible( list(plot.x=x, plot.y=y*scale+y0, x=x, y=y, scale=scale, y0=y0 ))
 
}
