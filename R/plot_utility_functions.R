#' Utility function.
#'
#' @description Wrapper around \code{\link[grDevices]{adjustcolor}}.
#' 
#' @param x A color or a vector with color values.
#' @param f A number for adjusting the transparency ranging from 0 (completely 
#' transparent) to 1 (not transparent).
#'
#' @family Utility functions for plotting

alpha <- function(x, f = 0.5) {
    if(f > 1 | f < 0){
        stop("Transparency value should be in range 0 to 1.")
    }else{
        return( adjustcolor(x, alpha.f = f) )
    }
}

#' Utility function.
#'
#' @export
#' @description Generate an color palette with changing transparency.
#' 
#' @param x A vector with color values. Could be a single value specifying a 
#' single color palette, ranging in transparency values, or a vector with 
#' different colors. 
#' @param f.seq A vector with transparency values, ranging from 0 to 1.
#' @param n Optional argument. A number specifying the number of colors in the 
#' palette. If \code{n} > 1, then N transparency values are generated ranging  
#' from the minimum of \code{f.seq} to the maximum of \code{f.seq}. \code{n} 
#' will only be used when the vector \code{f.seq} has two elements or more.
#' @return A vector with color values.
#' @author Jacolien van Rij
#' @section Warning:
#' On Linux \code{\link{x11}} devices may not support transparency. 
#' In that case, a solution might be to write the plots immediately to a file 
#' using functions such as \code{\link{pdf}}, or \code{\link{png}}.
#' @seealso 
#' \code{\link[grDevices]{palette}}, \code{\link[grDevices]{colorRampPalette}},
#' \code{\link[grDevices]{adjustcolor}}, \code{\link[grDevices]{convertColor}}
#' @examples 
#' # a palette of 5 white transparent colors:
#' alphaPalette('white', f.seq=1:5/5)
#' # the same palette:
#' alphaPalette('white', f.seq=c(.2,1), n=5)
#' # a palette with 10 colors blue, yellow and red, that differ in transparency
#' alphaPalette(c('blue', "yellow", "red"), f.seq=c(0.1,.8), n=10)
#'
#' @family Utility functions for plotting

alphaPalette <- function(x, f.seq, n=NULL) {
    out <- c()
    if(!is.null(n)){
        if(n[1]>1 & length(f.seq) > 1){
            f.seq <- seq(min(f.seq), max(f.seq), length=n)
        }else{
            n <- length(f.seq)
            warning("Argument n will be ignored.")
        }
    }else{
    	n <- length(f.seq)
    }
    if (length(x) == length(f.seq)) {
        out <- x
    } else if(length(x) == 1){
        out <- rep(x[1], length(f.seq))
    }else{
        x <- colorRampPalette(x)(n)
        out <- x
    }

    return(mapply(function(a, b) {
        alpha(a, b)
    }, out, f.seq, USE.NAMES = FALSE))
}



#' Utility function.
#'
#' @export
#' @description Add a transparency Rug to a contour plot or image.
#' 
#' @param x Observations on x-axis.
#' @param y Observations on y-axis.
#' @param n.grid Resolution of Rug. Defaults to 30, 
#' which means that the x- and y-axis are divided in 30 bins.
#' @param gradual Logical: whether or not to use the number of 
#' observations in an area, i.e., more transparent equals more 
#' observations. Default is FALSE, which means that the function only 
#' distinguishes between observations in a certain region or not, 
#' regardless how many observations.
#' @param max.alpha Maximum of transparency, number between 0 (completely 
#' transparent) and 1 (non-transparent). Defaults to .75.
#' @param col Color value. Defaults to "white".
#' @return Plots a shaded image over the contour plot or image.
#' @author Jacolien van Rij
#' @section Warning:
#' On Linux \code{\link{x11}} devices may not support transparency. 
#' In that case, a solution might be to write the plots immediately to a file 
#' using functions such as \code{\link{pdf}}, or \code{\link{png}}.
#' @seealso 
#' \code{\link[graphics]{rug}}, \code{\link[graphics]{contour}}, 
#' \code{\link[graphics]{image}}
#'
#' @family Utility functions for plotting

fadeRug <- function(x, y, n.grid = 30, gradual=FALSE, max.alpha = 0.75, col='white') {
    xlim <- c(par()$usr[1], par()$usr[2])
    ylim <- c(par()$usr[3], par()$usr[4])
    x.step <- length(xlim[1]:xlim[2])/n.grid
    y.step <- length(ylim[1]:ylim[2])/n.grid

    im <- matrix(table(factor(round((x - xlim[1])/x.step), levels = 1:n.grid), factor(round((y - ylim[1])/y.step), levels = 1:n.grid)))
    if(gradual==FALSE){
        im[im > 0] <- 1 
    }
    fadecols <- alphaPalette(col, f.seq = seq(max.alpha, 0, length = max(im) + 1))
    im <- matrix(fadecols[im + 1], byrow = TRUE, ncol = n.grid)
    rasterImage(as.raster(im), xleft = xlim[1], xright = xlim[2], ybottom = ylim[2], ytop = ylim[1], interpolate = FALSE)
}

#' Utility function.
#'
#' @description Function for positioning a legend or label in or outside the 
#' plot region based on proportion of the plot region rather than Cartesian 
#' coordinates.
#' 
#' @param pos A number indicating the proportion on the x-axis. Default is 1.1.
#' @param side Which axis to choose: 1=bottom, 2=left, 3=top, 4=right. Default is 1.
#' @author Jacolien van Rij
#'
#' @family Utility functions for plotting

getCoords <- function(pos = 1.1, side = 1) {
    p <- par()
    x.width = p$usr[2] - p$usr[1]
    y.width = p$usr[4] - p$usr[3]
    out <- rep(NA, length(pos))
    if(length(side)==1){
        side <- rep(side, length(pos))
    }
    out[which(side %in% c(1,3))] <- pos[which(side %in% c(1,3))] * x.width + p$usr[1]
    out[which(side %in% c(2,4))] <- pos[which(side %in% c(2,4))] * y.width + p$usr[3]
    return(out)
} 


#' Utility function.
#'
#' @export
#' @description Add a gradient legend to a contour plot (or other plot) to 
#' indicate the range of values represented by the color palette.
#' 
#' @param valRange Range of the values that is represented by the color 
#' palette. Normally two value-vector. If a larger vector is provided, only 
#' the min and max values are being used.
#' @param pos A number indicating the position on the axis in proportion. 
#' Using the arguments \code{length} and \code{depth} and \code{side} the 
#' position of the legend is calculated automatically. Alternatively, one 
#' could provide  a vector with 4 numbers, providing the xleft, ybottom, 
#' xright, ytop of a rectangle. These 4 points are indicated in proportions of 
#' the x- and y-axis. However, if the argument \code{coords} is set to TRUE, 
#' these positions are taken as values in the Cartesian coordinate system of 
#' the plot. Note: \code{coords} is only considered for 4-number vectors of 
#' \code{pos}.
#' @param side Which axis to choose: 1=bottom, 2=left, 3=top, 4=right.
#' Default is 4.
#' @param length Number, indicating the width of the legend as proportion with 
#' respect to the axis indicated by \code{side}. 
#' Note: when \code{pos} is defined by 4 numbers, \code{length} is ignored.
#' @param depth Number, indicating the height of the legend as proportion 
#' with respect to the axis perpendicular to \code{side}.
#' Note: when \code{pos} is defined by 4 numbers, \code{depth} is ignored.
#' @param inside Logical: whether or not to plot the legend inside or outside 
#' the plot area.
#' Note: when \code{pos} is defined by 4 numbers, \code{inside} is ignored.
#' @param coords Logical: whether or not \code{pos} is defined as coordinates. 
#' When FALSE, the default, \code{pos} is defined in proportions. 
#' Note: when \code{pos} is defined by 1 number, \code{inside} is ignored.
#' @param color Name of color palette to use ('topo', 'terrain', 'heat', 
#' 'rainbow'). Custom color palettes can also be provided, but then the 
#' argument \code{nCol} is ignored.
#' @param nCol Number of colors in the color palette.
#' @param n.seg Number of ticks and markers on the scale.
#' @param border.col Color of the border and the ticks.
#' @author Jacolien van Rij
#' @examples
#' data(simdat)
#' # simple GAM model:
#' m1 <- bam(Y~te(Time, Trial), data=simdat)
#'
#' # The functions pvisgam and fvisgam automatically plot legend,
#' # but vis.gam does not:
#' vis.gam(m1, view=c("Time", "Trial"), plot.type='contour', color='topo',
#' zlim=c(-14,14) )
#' gradientLegend(valRange=c(-14,14),pos=.5, side=3)
#' gradientLegend(valRange=c(-14,14),pos=.125, side=4, inside=FALSE)
#' gradientLegend(valRange=c(-14,14),pos=.75, length=.5,
#' color=alphaPalette('white', f.seq=seq(0,1, by=.1)), border.col='white')
#' # when defining custom points, it is still important to specify side:
#' gradientLegend(valRange=c(-14,14), pos=c(500,-5,1250,-4), coords=TRUE, 
#' border.col='red', side=1)
#' 
#' @family Utility functions for plotting

gradientLegend <- function(valRange, color='topo', nCol=30, 
    pos=.5, side=4, length=.25, depth=.05, inside=TRUE, coords=FALSE,
    n.seg=3, border.col = 'black'){

    loc   <- c(0,0,0,0)
    sides <- c(0,0,0,0)

    # case when length pos = 1: relative to axis and parallel, width, height
    if(length(pos)==1){
        pos.other <- ifelse(side>2, 1, 0)   

        if(side %in% c(1,3)){
            switch <- ifelse(inside, 0,1)
            switch <- ifelse(side>2, 1-switch, switch)
            loc <- getCoords(c(pos-.5*length,
                pos.other-switch*depth,
                pos+.5*length,
                pos.other+(1-switch)*depth), 
                side=c(side,2, side,2))
        }else if(side%in% c(2,4)){
            switch <- ifelse(inside, 0,1)
            switch <- ifelse(side>2, 1-switch, switch)
            loc <- getCoords(
                c(pos.other-switch*depth,
                  pos-.5*length,
                  pos.other+(1-switch)*depth,
                  pos+.5*length),
                side=c(1, side,1,side))
        }
    }else if(length(pos)==4){
        if(coords){
            loc <- pos
        }else{
            loc <- getCoords(pos, side=c(1,2,1,2))
        }
    }

    mycolors <- c()

    if(length(color)>1){
        mycolors <- color
    }else if (!is.null(nCol)){
        if(color=='topo'){
            mycolors <- topo.colors(nCol)
        }else if(color=='heat'){
            mycolors <- heat.colors(nCol)
        }else if(color=='terrain'){
            mycolors <- terrain.colors(nCol)
        }else if(color=='rainbow'){
            mycolors <- rainbow(nCol)
        }else{
            warning('Color %s not recognized. A palette of topo.colors is used instead.')
            mycolors <- topo.colors(nCol)
        }
    }else{
        stop('No color palette provided.')
    }


    vals <- seq(min(valRange), max(valRange), length=length(mycolors))

    im <- as.raster(mycolors[matrix(1:length(mycolors), ncol=1)])

    if(side%%2==1){
        rasterImage(t(im), loc[1],loc[2], loc[3], loc[4], col=mycolors, xpd=T)
        rect(loc[1],loc[2], loc[3], loc[4], border=border.col, xpd=T)
        ticks <- seq(loc[1],loc[3], length=n.seg)
        text(y=loc[4], x=ticks, 
            labels=seq(min(valRange),max(valRange), length=n.seg),
            pos=3, cex=.8, xpd=T)
        segments(x0=ticks, x1=ticks, y0=rep(loc[2], n.seg), y1=rep(loc[4], n.seg), col=border.col, xpd=TRUE)
    }else{
        rasterImage(rev(im), loc[1],loc[2], loc[3], loc[4], col=mycolors, xpd=T)
        rect(loc[1],loc[2], loc[3], loc[4], border=border.col, xpd=T)
        ticks <- seq(loc[2],loc[4], length=n.seg)
        text(y=ticks, x=loc[3], 
            labels=seq(min(valRange),max(valRange), length=n.seg),
            pos=4, cex=.8, xpd=T)
        segments(x0=rep(loc[1], n.seg), x1=rep(loc[3], n.seg), y0=ticks, y1=ticks, col=border.col, xpd=TRUE)
    }

}


#' Utility function.
#'
#' @export
#' @description Add vertical error bars.
#' 
#' @param x Vector with x-values.
#' @param mean Vector with means.
#' @param se Vector with errors or confidence bands.
#' @param minmax Optional argument, vector with two values indicating the 
#' minimum and maximum value for the error bars. If NULL (default) the error 
#' bars are not corrected.
#' @param ... Optional graphical parameters (see \code{\link[graphics]{par}}).
#' @author Jacolien van Rij
#'
#' @family Utility functions for plotting

errorBars <- function(x, mean, se, minmax=NULL, ...){
    cu <- mean+se
    cl <- mean-se
    if(!is.null(minmax)){
        cu[!is.na(cu) & cu>minmax[2]] <- minmax[2]
        cl[!is.na(cl) & cl<minmax[1]] <- minmax[1]
    }
    dnm <- names(list(...))
    pars <- list()

    if(!"length" %in% dnm){
        pars[['length']] <- .1
    }
    if(!"angle" %in% dnm){
        pars[["angle"]] <- 90
    }
    if(!"code" %in% dnm){
        pars[["code"]] <- 3
    }

    if(length(pars) > 0){
        pars <- paste(paste(names(pars),pars, sep='='), collapse=',')
        eval(parse(text=paste('arrows(x0=x, x1=x, y0=cl, y1=cu,', pars ,',...)', sep='')))
    }else{
        arrows(x0=x, x1=x, y0=cl, y1=cu,...)
    }        
}

#' Utility function.
#'
#' @description Add horizontal error bars.
#' 
#' @param y Vector with y-values.
#' @param mean Vector with means.
#' @param se Vector with errors or confidence bands.
#' @param minmax Optional argument, vector with two values indicating the 
#' minimum and maximum value for the error bars. If NULL (default) the error 
#' bars are not corrected.
#' @param ... Optional graphical parameters (see \code{\link[graphics]{par}}).
#' @author Jacolien van Rij
#'
#' @family Utility functions for plotting

horiz_error <- function(y, mean, se, minmax=NULL, ...){
    cr <- mean+se
    cl <- mean-se
    if(!is.null(minmax)){
        cr[!is.na(cr) & cr>minmax[2]] <- minmax[2]
        cl[!is.na(cl) & cl<minmax[1]] <- minmax[1]
    }
    dnm <- names(list(...))
    pars <- list()

    if(!"length" %in% dnm){
        pars[['length']] <- .1
    }
    if(!"angle" %in% dnm){
        pars[["angle"]] <- 90
    }
    if(!"code" %in% dnm){
        pars[["code"]] <- 3
    }

    if(length(pars) > 0){
        pars <- paste(paste(names(pars),pars, sep='='), collapse=',')
        eval(parse(text=paste('arrows(x0=cl, x1=cr, y0=y, y1=y,', pars ,',...)', sep='')))
    }else{
        arrows(x0=cl, x1=cr, y0=y, y1=y,...)
    }        
}


#' Utility function.
#'
#' @description Add horizontal or vertical interval indications. 
#' This function can also be used to plot asymmetric (non-parametric)
#' error bars or confidence intervals. Basically a wrapper around arrows.
#' 
#' @param pos Vector with x- or y-values (depending on \code{horizontal}).
#' @param lowVals Vector with low values, .
#' @param highVals Vector with errors or confidence bands.
#' @param horiz Logical: whether or not to plot the intervals horizontally.
#' Defaults to TRUE (horizontal intervals).
#' @param minmax Optional argument, vector with two values indicating the 
#' minimum and maximum value for the error bars. If NULL (default) the error 
#' bars are not corrected.
#' @param length Number, size of the edges in inches.
#' @param ... Optional graphical parameters (see \code{\link[graphics]{par}})
#' to be forwarded to the function \code{\link[graphics]{arrows}}.
#' @author Jacolien van Rij
#' @examples
#' emptyPlot(1000,5, xlab='Time', ylab='Y')
#' # add interval indication for Time=200 to Time=750:
#' addInterval(1, 200, 750, lwd=2, col='red')
#'
#' # zero-length intervals also should work:
#' addInterval(pos=521, lowVals=c(1.35, 1.5, 4.33), highVals=c(1.15,1.5, 4.05),
#'     horiz=FALSE, length=.1, lwd=4)
#'
#' # combine with getCoords for consistent positions with different axes:
#' par(mfrow=c(2,2))
#' # 1st plot:
#' emptyPlot(1000,c(-1,5), h0=0)
#' addInterval(getCoords(.1,side=2), 200,800, 
#'     col='red', lwd=2)
#' addInterval(getCoords(.5,side=1), 1,4, horiz=FALSE,
#'     col='blue', length=.15, angle=100, lwd=4)
#' abline(h=getCoords(.1, side=2), lty=3, col='red', xpd=TRUE)
#' abline(v=getCoords(.5, side=1), lty=3, col='blue', xpd=TRUE)
#' # 2nd plot:
#' emptyPlot(1000,c(-250, 120), h0=0)
#' addInterval(getCoords(.1,side=2), 750,1200, 
#'     col='red', lwd=2, minmax=c(0,1000))
#' abline(h=getCoords(.1, side=2), lty=3, col='red', xpd=TRUE)
#' # 3rd plot:
#' emptyPlot(c(-50,50),c(20,120), h0=0)
#' addInterval(getCoords(.5,side=1), 80,120, horiz=FALSE,
#'     col='blue', code=2, length=.15, lwd=4, lend=1)
#' abline(v=getCoords(.5, side=1), lty=3, col='blue', xpd=TRUE)
#'
#' # Plot boxplot hinges with medians:
#' data(simdat)
#' b <- boxplot(simdat$Y ~ simdat$Condition, plot=FALSE)$stats
#' emptyPlot(c(1,6), range(b[c(2,4),]), h0=0)
#' addInterval(1:6,b[2,], b[4,], horiz=FALSE)
#' # reset
#' par(mfrow=c(1,1))
#' @family Utility functions for plotting


addInterval <- function(pos, lowVals, highVals, horiz=TRUE, minmax=NULL, length=.05,...){
    convert2num <- function(x){  
        if(!any(c("numeric", "integer") %in% class(x))){
            if("factor" %in% class(x)){
                return( as.numeric( as.character(x)) )
            }else{
                return( as.vector(unlist(x)) )
            }
        }else{
            return(x)
        }
    }

    pos <- convert2num(pos)
    lowVals <- convert2num(lowVals)
    highVals <- convert2num(highVals)

    if(!is.null(minmax)){
        lowVals[!is.na(lowVals) & lowVals < minmax[1]] <- minmax[1]
        highVals[!is.na(highVals) & highVals > minmax[2]] <- minmax[2]
    }
    if(length(lowVals) != length(highVals)){
        if(length(lowVals)==1){
            lowVals <- rep(lowVals, length(highVals))
        }else if(length(highVals)==1){
            highVals <- rep(highVals, length(lowVals))          
        }else{
            stop('Vectors lowVals and highVals do not have same length.')
        }
    }
    if(length(pos)==1){
        pos <- rep(pos, length(lowVals))
    }else if(length(pos) != length(lowVals)){
        stop('Vector pos should be of the same length as lowVals and highVals.')
    }

    dnm <- names(list(...))
    pars <- list()

    if(!"angle" %in% dnm){
        pars[["angle"]] <- 90
    }
    if(!"code" %in% dnm){
        pars[["code"]] <- 3
    }

    if(length(pars) > 0){
        pars <- paste(paste(names(pars),pars, sep='='), collapse=',')
    }else{
        pars <- NULL
    }

    len.check <- highVals - lowVals
    len.check <- which(len.check == 0)

    if(length(len.check)>0){
        usr <- par()$usr
        pin <- par()$pin
        
        if(horiz){
            len <- ((usr[4]-usr[3])*length) / pin[2]
            segments(x0=lowVals[len.check], x1=lowVals[len.check], y0=pos-len, y1=pos+len, ...)
        }else{
            len <- ((usr[2]-usr[1])*length) / pin[1]
            segments(y0=lowVals[len.check], y1=lowVals[len.check], x0=pos-len, x1=pos+len, ...)
        }

        pos <- pos[-len.check]
        lowVals <- lowVals[-len.check]
        highVals <- highVals[-len.check]
    }

    if(horiz){
        if(is.null(pars)){
            arrows(x0=lowVals, x1=highVals, y0=pos, y1=pos, length=length, ...)
        }else{
            eval(parse(text=paste('arrows(x0=lowVals, x1=highVals, y0=pos, y1=pos,length=length,', pars ,',...)', sep='')))
        }
    }else{
        if(is.null(pars)){
            arrows(y0=lowVals, y1=highVals, x0=pos, x1=pos, length=length, ...)
        }else{
            eval(parse(text=paste('arrows(y0=lowVals, y1=highVals, x0=pos, x1=pos, length=length,', pars ,',...)', sep='')))
        }        
    }
}





