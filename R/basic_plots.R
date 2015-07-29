#' Utility function
#' 
#' @description Generate an empty plot window.
#' 
#' @export
#' @import grDevices
#' @import graphics
#' @param xlim A one- or two-value vector indicating the range of the x-axis. 
#' If \code{xlim} is a number, then it is assumed that the other value is 0. 
#' Thus, \code{xlim=3000} wil result in a x-axis ranging from 0 to 3000, 
#' and \code{xlim=-3} will result in a x-axis ranging from -3 to 0.
#' @param ylim A one- or two-value vector indicating the range of the y-axis. 
#' (See \code{xlim}) for more information.
#' @param main Title for the plot. Empty by default. 
#' Note: the title can be added later using \code{\link[graphics]{title}}.
#' @param xlab Label for x-axis. Empty by default. If no label is provided, 
#' use \code{\link[graphics]{mtext}} for adding axis labels.
#' @param ylab Label for y-axis. Empty by default. (See \code{xlab}.)
#' @param h0 A vector indicating where to add solid horizontal lines for 
#' reference. By default no values provided.
#' @param v0 A vector indicating where to add dotted vertical lines for 
#' reference. By default no values provided.
#' @param bty A character string which determined the type of box which is 
#' drawn about plots. If bty is one of "o", "l", "7", "c", "u", or "]" the 
#' resulting box resembles the corresponding upper case letter. A value of 
#' "n"  (the default) suppresses the box.
#' @param eegAxis Logical: whether or not to reverse the y-axis, plotting the 
#' negative amplitudes upwards as traditionally is done in EEG research.
#' If eeg.axes is TRUE, labels for x- and y-axis are provided, when not 
#' provided by the user. Default value is FALSE.
#' @param ... Other arguments for plotting, see \code{\link[graphics]{par}}.
#' @return An empty plot window.
#' @author Jacolien van Rij
#' @seealso Use \code{\link[graphics]{title}} and 
#' \code{\link[graphics]{mtext}}  for drawing labels and titles; 
#' use  \code{\link[graphics]{lines}} and \code{\link[graphics]{points}} 
#' for plotting the data; 
#' use \code{\link[graphics]{legend}} for adding a legend.
#' and \code{\link{acf_n_plots}} for inspection of individual time series.
#' @examples
#' data(simdat)
#' test <- simdat[simdat$Subject=='a10' & simdat$Trial==10,]
#' 
#' emptyPlot(range(test$Time), range(test$Y),
#' main='Data', ylab='Y', xlab='Time', 
#' h0=0, v0=c(0,1000,2000))
#' # Note that this is the same as:
#' emptyPlot(range(test$Time), range(test$Y))
#' title(main='Data', ylab='Y', xlab='Time')
#' abline(h=0)
#' abline(v=c(0,1000,2000), lty=3)
#' 
#' # To add data, use lines() and points()
#' lines(test$Time, test$Y)
#'
#' @family Utility functions for plotting

emptyPlot <- function(xlim, ylim, 
    main=NULL, xlab=NULL, ylab=NULL, h0=NULL, v0=NULL, 
    bty='n', eegAxis=FALSE, ...){

    if(length(xlim)==1){
        xlim <- sort(c(0,xlim))
    }
    if(length(ylim)==1){
        ylim <- sort(c(0,ylim))
    }
    if(is.null(main)){
        main=''
    }
    if(eegAxis){
        ylim <- sort(ylim, decreasing=T)
        if(is.null(ylab)){
            ylab=expression(paste('Amplitude (', mu, V,')', sep=''))
        }
        if(is.null(xlab)){
            xlab="Time (ms)"
        }
    }else{
        if(is.null(ylab)){
            ylab=''
        }
        if(is.null(xlab)){
            xlab=""
        }
    }

    plot(range(xlim), range(ylim), type='n',
        xlim=xlim, ylim=ylim, 
        main=main, xlab=xlab, ylab=ylab,
        bty=bty, ...)
    if(!is.null(h0)){
        abline(h=h0)
    }
    if(!is.null(v0)){
        abline(v=v0, lty=3)
    }
}

#' Utility function
#' 
#' @description Plot line with confidence intervals.
#' 
#' @export
#' @import grDevices
#' @import graphics
#' @param x Vector with values on x-axis.
#' @param fit Vector with values on y-axis.
#' @param se.fit Vector with standard error; or when \code{se.fit2}
#' is provided, \code{se.fit} specifies upper values confidence
#' interval.
#' @param se.fit2 Optional: lower values confidence interval.
#' @param shade Logical: whether or not to produce shaded regions as 
#' confidence bands.
#' @param f Factor for converting standard error in confidence intervals. 
#' Defaults to 1. Use 1.96 for 95\% CI, and 2.58 for 99\% CI.
#' @param col Color for lines and confindence bands.
#' @param border A number, indicating the linewidth of the border around 
#' shaded area. No border with value NA (default). Only applies when 
#' \code{shade=TRUE}. 
#' @param alpha Transparency of shaded area. Number between 0 
#' (completely transparent) and 1 (not transparent). 
#' @param ... Optional arguments for the lines. See \code{\link{par}}.
#' @author Jacolien van Rij
#'
#' @examples
#' data(simdat)
#' 
#' # Use aggregate to calculate mean and standard deviation per timestamp:
#' avg <- aggregate(simdat$Y, by=list(Time=simdat$Time),
#'     function(x){c(mean=mean(x), sd=sd(x))})
#' head(avg)
#' # Note that column x has two values in each row:
#' head(avg$x)
#' head(avg$x[,1])
#' 
#' # Plot line and standard deviation:
#' emptyPlot(range(avg$Time), c(-20,20), h0=0)
#' plot_error(avg$Time, avg$x[,'mean'], avg$x[,'sd'])
#' # Change layout:
#' emptyPlot(range(avg$Time), c(-20,20), h0=0)
#' plot_error(avg$Time, avg$x[,'mean'], avg$x[,'sd'], 
#'    shade=TRUE, lty=3, lwd=3)
#'
#' # Show difference with 0:
#' x <- find_difference(avg$x[,'mean'], avg$x[,'sd'], xVals=avg$Time)
#' # Add arrows:
#' abline(v=c(x$start, x$end), lty=3, col='red')
#' arrows(x0=x$start, x1=x$end, y0=-5, y1=-5, code=3, length=.1, col='red')
#'
#' # Use of se.fit2:
#' avg$cu <- avg$x[,'mean'] + avg$x[,'sd']
#' avg$cl <- avg$x[,'mean'] - avg$x[,'sd']
#' emptyPlot(range(avg$Time), c(-20,20), h0=0)
#' plot_error(avg$Time, avg$x[,'mean'], se.fit=avg$cu, se.fit2=avg$cl, col='red')
#' # Change layout:
#' emptyPlot(range(avg$Time), c(-20,20), h0=0)
#' plot_error(avg$Time, avg$x[,'mean'], avg$x[,'sd'], 
#'    shade=TRUE, lty=3, lwd=3)
#'
#' @family Utility functions for plotting

plot_error <- function(x, fit, se.fit, se.fit2=NULL, 
    shade=FALSE, f=1, col='black', border=NA, alpha=.25,  ...){
    if(shade){
        xval <- c(x, rev(x))
        yval <- NA
        if(is.null(se.fit2)){
            yval <- c(fit+f*se.fit, rev(fit-f*se.fit))
        }else{
            yval <- c(se.fit, rev(se.fit2))
        }
        polygon(x=xval, y=yval, border=border, col=alpha(col, f=alpha))
    }else{
        if(is.null(se.fit2)){
            lines(x, fit+f*se.fit, lty=2, col=col, ...)
            lines(x, fit-f*se.fit, lty=2, col=col, ...)
        }else{
            lines(x, se.fit, lty=2, col=col, ...)
            lines(x, se.fit2, lty=2, col=col, ...)
        } 
    }
    lines(x, fit, col=col, ...)

}

#' Utility function
#' 
#' @description Fill area under line or plot.
#' 
#' @export
#' @import grDevices
#' @import graphics
#' @param x Vector with values on x-axis.
#' @param y Vector with values on y-axis.
#' @param from A number indicating until which value on the y-axis the graph 
#' is colored. Defaults to 0.
#' @param col Color for filling the area. Default is black.
#' @param alpha Transparency of shaded area. Number between 0 
#' (completely transparent) and 1 (not transparent). Default is .25.
#' @param border A color, indicating the color of the border around 
#' shaded area. No border with value NA (default). 
#' @param na.rm Logical: whether or not to remove the missing values in 
#' \code{x} and \code{y}. Defaults to TRUE. If set to FALSE, missing values 
#' may cause that the filled area is split in various smaller areas.
#' @param horiz Logical: whether or not to plot with respect to the 
#' x-axis (TRUE) or y-xis (FALSE). Defaults to TRUE.
#' @param ... Optional arguments for the lines. See \code{\link{par}}.
#' @author Jacolien van Rij
#' @examples
#' # density of a random sample from normal distribution:
#' test <- density(rnorm(1000))
#' emptyPlot(range(test$x), range(test$y))
#' fill_area(test$x, test$y)
#' fill_area(test$x, test$y, from=.1, col='red')
#' fill_area(test$x, test$y, from=.2, col='blue', density=10, lwd=3)
#' lines(test$x, test$y, lwd=2)
#' 
#' @seealso \code{\link{check_normaldist}}
#' @family Utility functions for plotting

fill_area <- function(x, y, from=0, col='black', alpha=.25,  border=NA, na.rm=TRUE, 
    horiz=TRUE, outline=FALSE, ...){
    el.narm <- c()
    if(na.rm){
        el.narm <- which(is.na(x) | is.na(y))
        if(length(from)==length(x)){
            from = from[!is.na(x) | !is.na(y)]
        }
        x <- x[!is.na(x) | !is.na(y)]
        y <- y[!is.na(x) | !is.na(y)]
    }

    xval <- c(x, rev(x))
    yval <- c()
    if(length(from)==1){
        yval <- c(y, rep(from, length(y)))
    }else if(length(from)==length(x)){
        yval <- c(y, rev(from))
    }else{
        warning("Argument from has more than 1 element. Only first element being used.")
        yval <- c(y, rep(from, length(y)))
    }

    if(names(dev.cur())[1] %in% c("X11", "postscript", "xfig", "pictex") ){
        alpha = 1
    }

    line.args <- list2str(x=c("type", "pch", "lty", "bg", "cex", "lwd", "lend", "ljoin", "lmitre"), inputlist=list(...))
    fill.args <- list2str(x= c("density", "angle", "lty", "fillOddEven", "lwd", "lend", "ljoin", "lmitre"), inputlist=list(...))
    

    if(horiz){  
        if(!is.na(border) ){
            if( outline==TRUE){
                eval(parse(text=sprintf("polygon(x=xval, y=yval, border=border, col=alpha(col, f=alpha), %s, xpd=TRUE)", fill.args )  ))
            }else{
                eval(parse(text=sprintf("polygon(x=xval, y=yval, border=NA, col=alpha(col, f=alpha), %s, xpd=TRUE)", fill.args )  ))
                eval(parse(text=sprintf("lines(x=x, y=y, col=border, %s, xpd=TRUE)", line.args )  ))
            }
        }else{
            eval(parse(text=sprintf("polygon(x=xval, y=yval, border=NA, col=alpha(col, f=alpha), %s, xpd=TRUE)", fill.args )  ))
        }
    }else{  
        if(!is.na(border) ){
            if( outline==TRUE){
                eval(parse(text=sprintf("polygon(x=yval, y=xval, border=border, col=alpha(col, f=alpha), %s, xpd=TRUE)", fill.args )  ))
            }else{
                eval(parse(text=sprintf("polygon(x=yval, y=xval, border=NA, col=alpha(col, f=alpha), %s, xpd=TRUE)", fill.args )  ))
                eval(parse(text=sprintf("lines(x=y, y=x, col=border, %s, xpd=TRUE)", line.args )  ))
            }
        }else{
            eval(parse(text=sprintf("polygon(x=yval, y=xval, border=NA, col=alpha(col, f=alpha), %s, xpd=TRUE)", fill.args )  ))
        }     
    }
}


#' Utility function
#' 
#' @export
#' @export
#' @import grDevices
#' @import graphics
#' @import stats
#' @description Adjusted version of the a Cleveland dot plot implemented in 
#' \code{\link[graphics]{dotchart}} with the option to add confidence 
#' intervals.
#' 
#' @param x  either a vector or matrix of numeric values (NAs are allowed). 
#' If x is a matrix the overall plot consists of juxtaposed dotplots for each 
#' row. Inputs which satisfy is.numeric(x) but not is.vector(x) || is.matrix(
#' x) are coerced by as.numeric, with a warning.
#' @param se.val a vector or matrix of numeric values representing the 
#' standard error or confidence bands.
#' @param labels a vector of labels for each point. For vectors the default is 
#' to use names(x) and for matrices the row labels dimnames(x)[[1]].
#' @param groups  an optional factor indicating how the elements of x are 
#' grouped. If x is a matrix, groups will default to the columns of x.
#' @param gdata data values for the groups. This is typically a summary such 
#' as the median or mean of each group.
#' @param cex the character size to be used. Setting cex to a value smaller
#' than one can be a useful way of avoiding label overlap. Unlike many other 
#' graphics functions, this sets the actual size, not a multiple of par("cex").
#' @param pch the plotting character or symbol to be used.
#' @param gpch the plotting character or symbol to be used for group values.
#' @param bg  the background color of plotting characters or symbols to be 
#' used; use par(bg= *) to set the background color of the whole plot.
#' @param color the color(s) to be used for points and labels.
#' @param gcolor the single color to be used for group labels and values.
#' @param lcolor the color(s) to be used for the horizontal lines.
#' @param xlim horizontal range for the plot, see plot.window, e.g.
#' @param main overall title for the plot, see title.
#' @param xlab x-axis annotation as in title.
#' @param ylab y-axis annotation as in title.
#' @param lwd with of error bars.
#' @param ... graphical parameters can also be specified as arguments
#' see \code{\link[graphics]{par}}
#' @author This function is a slightly adjusted version of the function 
#' \code{\link[graphics]{dotchart}} of the package \code{\link{graphics}} 
#' version 3.1.1
#' @examples
#' data(simdat)
#' avg <- with(simdat, 
#'      aggregate(list(Y=Y), list(Group=Group, Condition=Condition), 
#'          function(x){c(mean=mean(x), sd=sd(x))}))
#' dotplot_error(avg$Y[,'mean'], se.val=avg$Y[,'sd'], 
#'     groups=avg$Group, labels=avg$Condition)
#' @seealso \code{\link[graphics]{dotchart}}
#' @family Utility functions for plotting

dotplot_error <- function (x, se.val=NULL, labels = NULL, groups = NULL, 
    gdata = NULL, cex = par("cex"), 
    pch = 21, gpch = 21, bg = "black", color = par("fg"), gcolor = par("fg"), 
    lcolor = "gray", xlim = NULL, main = NULL, 
    xlab = NULL, ylab = NULL, lwd=1, ...) 
{
    opar <- par("mai", "mar", "cex", "yaxs")
    on.exit(par(opar))
    par(cex = cex, yaxs = "i")
    if (!is.numeric(x)) 
        stop("'x' must be a numeric vector or matrix")
    n <- length(x)
    if(!is.null(se.val)){
        if(length(x) != length(se.val)){
            warning("se.val not equal in length as x. se.val will be ignored.")
            se.val <- NULL
        }
    }
    if (is.matrix(x)) {
        if (is.null(labels)) 
            labels <- rownames(x)
        if (is.null(labels)) 
            labels <- as.character(1L:nrow(x))
        labels <- rep_len(labels, n)
        if (is.null(groups)) 
            groups <- col(x, as.factor = TRUE)
        glabels <- levels(groups)
    }
    else {
        if (is.null(labels)) 
            labels <- names(x)
        glabels <- if (!is.null(groups)) 
            levels(groups)
        if (!is.vector(x)) {
            warning("'x' is neither a vector nor a matrix: using as.numeric(x)")
            x <- as.numeric(x)
        }
        if(! is.null(se.val)){
            if (!is.vector(se.val)) {
                warning("'se.val' is neither a vector nor a matrix: using as.numeric(se.val)")
                se.val <- as.numeric(se.val)
            }
        }

    }
    if(is.null(xlim)){
        xlim <- range(x[is.finite(x)])
        if(!is.null(se.val)){
            xlim <- range(c(x[is.finite(x)]-se.val[is.finite(se.val)], x[is.finite(x)]+se.val[is.finite(se.val)]))
        }
    }
    plot.new()
    linch <- if (!is.null(labels)) 
        max(strwidth(labels, "inch"), na.rm = TRUE)
    else 0
    if (is.null(glabels)) {
        ginch <- 0
        goffset <- 0
    }
    else {
        ginch <- max(strwidth(glabels, "inch"), na.rm = TRUE)
        goffset <- 0.4
    }
    if (!(is.null(labels) && is.null(glabels))) {
        nmai <- par("mai")
        nmai[2L] <- nmai[4L] + max(linch + goffset, ginch) + 
            0.1
        par(mai = nmai)
    }
    if (is.null(groups)) {
        o <- sort.list(as.numeric(x), decreasing = TRUE)
        x <- x[o]
        y <- 1L:n
        ylim <- c(0, n + 1)
    }
    else {
        o <- group_sort(x, group=groups, decreasing = TRUE)
        x <- x[o]
        if(!is.null(se.val)){
            se.val <- se.val[o]
        }
        groups <- groups[o]
        color <- rep_len(color, length(groups))[o]
        lcolor <- rep_len(lcolor, length(groups))[o]
        bg <- rep_len(bg, length(groups))[o]
        offset <- cumsum(c(0, diff(as.numeric(groups)) != 0))
        y <- 1L:n + 2 * offset
        ylim <- range(0, y + 2)
    }
    plot.window(xlim = xlim, ylim = ylim, log = "")
    lheight <- par("csi")
    if (!is.null(labels)) {
        linch <- max(strwidth(labels, "inch"), na.rm = TRUE)
        loffset <- (linch + 0.1)/lheight
        labs <- labels[o]
        mtext(labs, side = 2, line = loffset, at = y, adj = 0, 
            col = color, las = 2, cex = cex, ...)
    }
    abline(h = y, lty = "dotted", col = lcolor)
    if(!is.null(se.val)){
        segments(x0=x-se.val, x1=x+se.val, y0=y, y1=y, col=color, lwd=lwd)
    }
    points(x, y, pch = pch, col = color, bg = bg)
    if (!is.null(groups)) {
        gpos <- rev(cumsum(rev(tapply(groups, groups, length)) + 
            2) - 1)
        ginch <- max(strwidth(glabels, "inch"), na.rm = TRUE)
        goffset <- (max(linch + 0.2, ginch, na.rm = TRUE) + 0.1)/lheight
        mtext(glabels, side = 2, line = goffset, at = gpos, adj = 0, 
            col = gcolor, las = 2, cex = cex, ...)
        if (!is.null(gdata)) {
            abline(h = gpos, lty = "dotted")
            points(gdata, gpos, pch = gpch, col = gcolor, bg = bg, 
                ...)
        }
    }
    axis(1)
    box()
    title(main = main, xlab = xlab, ylab = ylab, ...)
    invisible()
}


#' Creates a contour plot with colored background.
#'
#' @description This function is a wrapper around \code{\link[graphics]{image}}
#' and \code{\link[graphics]{contour}}. See \code{vignette("plotfunctions")} 
#' for an example of how you could use \code{\link[graphics]{image}} and 
#' \code{\link[graphics]{contour}}.
#'
#' @export
#' @import grDevices
#' @import graphics
#' @param x Locations of grid lines at which the values in z are measured. 
#' These must be in ascending order. By default, equally spaced values from 0 
#' to 1 are used. If x is a list, its components x$x and x$y are used for x 
#' and y, respectively. If the list has component z this is used for z.
#' @param y Locations of grid lines at which the values in z are measured. 
#' @param z a matrix containing the values to be plotted (NAs are allowed). 
#' Note that x can be used instead of z for convenience.
#' @param main Text string, an overall title for the plot.
#' @param xlab Label for x axis. Default is name of first \code{view} variable.
#' @param ylab Label for y axis. Default is name of second \code{view} 
#' variable.
#' @param xlim x-limits for the plot.
#' @param ylim y-limits for the plot.
#' @param zlim z-limits for the plot.
#' @param col Color for the  contour lines and labels.
#' @param color a list of colors such as that generated by 
#' \code{\link[grDevices]{rainbow}}, \code{\link[grDevices]{heat.colors}}
#' \code{\link[grDevices]{colors}}, \code{\link[grDevices]{topo.colors}}, 
#' \code{\link[grDevices]{terrain.colors}} or similar functions.
#' @param nCol The number of colors to use in color schemes.
#' @param add.color.legend Logical: whether or not to add a color legend. 
#' Default is TRUE. If FALSE (omitted), one could use the function
#' \code{\link{gradientLegend}} to add a legend manually at any position.
#' @param ... Optional parameters for \code{\link[graphics]{image}}
#' and \code{\link[graphics]{contour}}.
#' @author Jacolien van Rij
#' @seealso \code{\link[graphics]{image}}, \code{\link[graphics]{contour}},
#' \code{\link[graphics]{filled.contour}}. See \code{\link{plotsurface}}
#' for plotting model predictions using \code{color_contour}.
#' @examples
#'
#' # Volcano example of R (package datasets)
#' color_contour(z=volcano)
#' # change color and lines:
#' color_contour(z=volcano, color='terrain', col=alpha(1), lwd=2, lty=5)
#' # change x-axis values and zlim:
#' color_contour(x=seq(500,700, length=nrow(volcano)),
#'     z=volcano, color='terrain', col=alpha(1), lwd=2, zlim=c(0,200))
#'
#' \dontrun{
#' # compare with similar functions:
#' filled.contour(volcano, color.palette=terrain.colors)
#' }
#' # without contour lines:
#' color_contour(z=volcano, color='terrain', lwd=0, drawlabels=FALSE)
#' # without background:
#' color_contour(z=volcano, color=NULL, add.color.legend=FALSE)
#' @family Utility functions for plotting


color_contour <- function(x = seq(0, 1, length.out = nrow(z)),
    y = seq(0, 1, length.out = ncol(z)),
    z,
    main=NULL, xlab=NULL, ylab=NULL, 
    xlim=NULL, ylim=NULL, zlim=NULL,
    col=NULL, color=topo.colors(50), nCol=50,
    add.color.legend=TRUE, ...){


    # check input:
    if(is.null(dim(z))){
        stop('z should be a matrix.')
    }
    if(length(x) != nrow(z)){
        stop(sprintf('x should have %d values, because z has %d rows.', nrow(z), nrow(z)))
    }
    if(length(y) != ncol(z)){
        stop(sprintf('y should have %d values, because z has %d columns.', ncol(z), ncol(z)))
    }
        

    ## Check plot settings
    if(is.null(main)){
        main=""
    }
    if(is.null(xlab)){
        xlab=""
    }
    if(is.null(ylab)){
        ylab=""
    }
    if(is.null(xlim)){
        xlim=range(x)
    }
    if(is.null(ylim)){
        ylim=range(y)
    }   
    if(is.null(zlim)){
        zlim=range(z)
    }   

    # colors:
    if(is.null(color)){
        color <- alphaPalette('white', f.seq=c(0,0), n=nCol)
    }
    if (color[1] == "heat") {
        color <- heat.colors(nCol)
        if(is.null(col)){
            col <- 3
        }
    } else if (color[1] == "topo") {
        color <- topo.colors(nCol)
        if(is.null(col)){
            col <- 2
        }
    } else if (color[1] == "cm") {
        color <- cm.colors(nCol)
        if(is.null(col)){
            col <- 1
        }
    } else if (color[1] == "terrain") {
        color <- terrain.colors(nCol)
        if(is.null(col)){
            col <- 2
        }
    } else if (color[1] == "bpy") {
        if (requireNamespace("sp", quietly = TRUE)) {
            color <- sp::bpy.colors(nCol)
            if(is.null(col)){
                col <- 3
            }
        } else {
            warning("Package 'sp' needed for bpy color palette. Using topo.colors instead (default).")
            color <- topo.colors(nCol)
            col <- 2
        }
    } else if (color[1] == "gray" || color[1] == "bw") {
        color <- gray(seq(0.1, 0.9, length = nCol))
        col <- 1
    } 
    if (is.null(col)){
        col <- 'black'
    } 


    dnm <- list(...)
    parlist <- names(dnm)

    type2string <- function(x){
        out <- ""
        if(length(x)>1){
            if(is.character(x)){
                out <- sprintf("c(%s)", paste(sprintf("'%s'", x), collapse=','))
            }else{
                out <- sprintf("c(%s)", paste(x, collapse=','))
            }
        }else{
            if(is.character(x)){
                out <- sprintf("'%s'", x)
            }else{
                out <- sprintf("%s", x)
            }
        }
        return(out)
    }

    # check contour input:
    cpar <- c()
    contourarg <- c('nlevels', 'levels', 'labels', 'labcex', 'drawlabels', 'method', 'lty', 'lwd')
    for(i in parlist[parlist %in% contourarg] ){
        cpar <- c(cpar, sprintf("%s=%s", i, type2string(dnm[[i]])))
    }
    cpar <- paste(",", paste(cpar, collapse=','))

    cpar2 <- c()
    for(i in parlist[parlist %in% c('nlevels', 'levels', 'method')] ){
        cpar2 <- c(cpar2, sprintf("%s=%s", i, type2string(dnm[[i]])))
    }
    cpar2 <- paste(",", paste(cpar2, collapse=','))

    # check image input:
    ipar <- c()
    contourarg <- c('nlevels', 'levels', 'labels', 'labcex', 'drawlabels', 'method', 'lty', 'lwd')

    for(i in parlist[!parlist %in% contourarg] ){
        ipar <- c(ipar, sprintf("%s=%s", i, type2string(dnm[[i]])))
    }

    ipar <- paste(",", paste(ipar, collapse=','))

    eval(parse(text=sprintf("image(x, y, z, col=color, xlim=xlim, ylim=ylim, zlim=zlim, main=main, xlab=xlab, ylab=ylab, add=FALSE%s)", ipar)))
    eval(parse(text=sprintf("contour(x, y, z, col=col, add=TRUE%s)",
        cpar)))

    if(add.color.legend){
        gradientLegend(round(zlim, 3), n.seg=3, pos=.875, 
            color=color)
    }
}