#' Utility function.
#' 
#' @export
#' @import stats
#' @param x A numeric vector.
#' @param element Logical: whether or not to return the value (FALSE, default) 
#' or the index (TRUE).
#' @return The value or index of the element closest to zero (absolute 
#' minimum).
#' @examples
#' (test <- seq(-25,25, by=3))
#' min(test[test>0])
#' max(test[test<0])
#' min(abs(test))
#' findAbsMin(test)
#' @author Jacolien van Rij
#' @family Data utility functions

findAbsMin <- function(x, element = FALSE) {
    abs_x <- abs(x)
    el <- max(which(abs_x == min(abs_x)))
    if (element) {
        return(el)
    } else {
        return(x[el])
    }
}

#' Utility function.
#' 
#' @export
#' @import stats
#' @param x A numeric vector.
#' @param dec Number of decimal points for rounding using function 
#' \code{\link[base]{round}}. Applied after argument 
#' \code{step}. If NULL (default), no rounding is applied.
#' @param step Round the 
#' @param n.seg Numeric value, number of values in the equally spaced sequence. 
#' Default is 2 (min, max).
#' @return vector, range of equally spaced sequence.
#' @examples
#' zlim <- c(-2.5, 3.01)
#' # does not change anything:
#' getRange(zlim)
#' # create a range of 5 numbers: 
#' # (basically just using seq )
#' getRange(zlim, n.seg=5)
#' # rounds the numbers:
#' getRange(zlim, dec=0)
#' getRange(zlim, n.seg=5, dec=0)
#' # extreme values are multiplications of 5
#' # that contains zlim values:
#' getRange(zlim, step=5)
#' getRange(zlim, step=5, n.seg=5)
#' # similar, but not the same:
#' getRange(zlim, n.seg=5, dec=0)
#' getRange(zlim, n.seg=5, step=1)
#' # combining:
#' getRange(zlim, n.seg=5, step=1, dec=0)
#' 
#' @author Jacolien van Rij
#' @family Data utility functions
getRange <- function(x, dec=NULL, step=NULL, n.seg=2){
	vals <- seq(min(x), max(x), length=n.seg)
    if (!is.null(step)){
        vals <- seq(floor(min(x)/step)*step, ceiling(max(x)/step)*step, length=n.seg)
    }
    if (!is.null(dec)){
        vals <- round(vals, dec)
    }
    return(vals)
}


#' Utility function.
#' 
#' @export
#' @import stats
#' @param x A numeric vector.
#' @return Number of decimals
#' 
#' @author Jacolien van Rij
#' @family Data utility functions
getDec <- function(x){
	dec <- 0

	dec <- sapply(x, function(a){
		numstr <- format(a, scientific=TRUE)
		numstr <- as.numeric( tail( unlist( strsplit(numstr, split='[e\\+]', fixed=FALSE) ), 1) )
		if(numstr < 0){
			return( abs(numstr)+1 )
		}else{
			return( 0 )
		}
	})
	
	return(dec)
}

#' Utility function.
#'
#' @export
#' @description Combine list values as string.
#' 
#' @param x A vector with the names or numbers of list elements to be combined.
#' @param inputlist A (named) list with information, e.g., graphical parameter settings.
#' @return String
#' @family Data utility functions

list2str <- function(x, inputlist) {
    out <- c()

    for(i in x){
        name.i <- NULL
        val.i  <- NULL

        if(is.numeric(i)){
            if(i > 0 & i <= length(inputlist)){
                name.i <- sprintf("el%.0f", i)
                val.i  <- inputlist[[i]]
            }
        }else if(i %in% names(inputlist)){
            name.i <- i
            val.i  <- inputlist[[i]]
        }
        if(! is.null(name.i)){
            if(inherits(val.i, c("numeric", "logical"))){
                out <- c(out, sprintf("%s=c(%s)", name.i, paste(val.i, collapse=",")))
            }else if(inherits(val.i, c("character", "factor"))){
                out <- c(out, sprintf("%s=c(%s)", name.i, paste(sprintf("'%s'",val.i), collapse=",")))
            }else{
                warning(sprintf("Class %s is not supported, element %s is ignored.",
                    class(name.i)[1], name.i))
            } 
        }
    }
    return(paste(out, collapse=", "))
}


#' Utility function.
#' 
#' @export
#' @import grDevices
#' @import graphics
#' @param el A numeric vector.
#' @param n Number indicating how many points around the elements of \code{el} 
#' need to be selected.
#' @param max The maximum value of the returned elements.
#' @return A vector with the elements of x surrounded by n points.
#' @examples
#' vectorIndices <- 1:1000
#' indOutliers <- c(2,10, 473, 359, 717, 519)
#' fn3 <- find_n_neighbors(indOutliers, n=3, max=max(vectorIndices))
#' fn20 <- find_n_neighbors(indOutliers, n=20, max=max(vectorIndices))
#'
#' # check fn3:
#' print(fn3)
#'
#' # Plot:
#' emptyPlot(c(-10,1000), c(-1,1), h0=0, v0=indOutliers)
#' points(fn3, rep(.5, length(fn3)), pch='*')
#' points(fn20, rep(-.5, length(fn20)), pch='*')
#'
#' @author Jacolien van Rij
#' @family Data utility functions

find_n_neighbors <- function(el, n, max) {
    if (length(el) > 0) {
        new_el <- sort(unique(unlist(lapply(el, FUN = function(x) {
            return(sort(unique(c(x, (x - n):x, x:(x + n)))))
        }))))
        new_el <- new_el[new_el >= 1 & new_el <= max]
        return(new_el)
    } else {
        return(NULL)
    }
}


#' Utility function.
#' 
#' @export
#' @param x A vector.
#' @param n Number indicating how many steps the vector should shift forward 
#' (N > 0) or backward (n < 0).
#' @param na_value The value to replace the empty cells with (e.g., the first 
#' or last points). Defaults to NA.
#' @return A vector with the same length of \code{x}, all moved \code{n} steps.
#' @examples
#' data(simdat)
#' (test <- simdat[1:20,])
#' test$Y.prev <- move_n_point(test$Y)
#' test$change <- test$Y - test$Y.prev
#' test$Y.post5 <- move_n_point(test$Y, n=-5)
#'
#' emptyPlot(nrow(test), range(test$Y))
#' lines(test$Y)
#' lines(test$Y.prev, col='red')
#' lines(test$Y.post5, col='blue')
#'
#' @author Jacolien van Rij
#' @family Data utility functions

move_n_point <- function(x, n = 1, na_value = NA) {
    x <- as.vector(x)
    
    if (length(x) > abs(n)) {
        if (n > 0) {
            return(c(rep(na_value, n), x[1:(length(x) - n)]))
        } else {
            return(c(x[(abs(n) + 1):(length(x))], rep(na_value, abs(n))))
        }
    } else if (length(x) == abs(n)) {
        return(NA)
    } else {
        return(NULL)
    }
} 

#' Calculate standard error of the mean.
#' 
#' @export
#' @param x A vector.
#' @return Standard Error of the mean.
#' @author Jacolien van Rij

se <- function(x){
    s <- sd(x)
    if (is.na(s)){
        warning("Problem in calculating SD.")
        return(NA)
    }
    else{ return(sd(x)/sqrt(length(x))) } 
}


#' Utility function.
#' 
#' @export
#' @description The function prints a summary of the data. 
#' Similar to the function \code{\link[utils]{str}}, but easier readable.
#' @param data A data frame.
#' @param print Logical: whether or not to print the summary.
#' @param n Number: maximum number of values being mentioned in the summary.
#' If NULL all values are being mentioned. Defaults to 10.
#' @return Optionally returns a named list with info.
#' @author Jacolien van Rij
#' @examples
#' data(simdat)
#' summary_data(simdat)
#' @family Data utility functions

summary_data <- function(data, print=TRUE, n=10){    
    
    labelColumns <- function(x, data){
        out <- NULL
        cn <- ifelse(is.numeric(x), colnames(data)[x], x)
        cl <- class(data[,cn])
        if(inherits(data[,cn],'factor')){
            vals <- sort(unique(as.character(data[,x])))
            n.cur <- length(vals)+1
            if(!is.null(n)){
                n.cur <- n
            }
            if(length(vals)>n.cur){
                out <- sprintf("factor with %d values; set to the value(s): %s, ...", 
                    length(vals),
                    paste( vals[1:n.cur], collapse=", ") )  
            }else{
                out <- sprintf("factor; set to the value(s): %s.", 
                    paste( vals, collapse=", ") )
            }
        }else if(inherits(data[,cn],'numeric')){
            if(length(unique(data[,x])) > 2){
                out <- sprintf("numeric predictor; with %d values ranging from %f to %f.", 
                    length(unique(data[,x])),
                    min(data[,x], na.rm=TRUE), max(data[,x], na.rm=TRUE)) 
            }else{
                out <- sprintf("numeric predictor; set to the value(s): %s.", paste(unique(data[,x]), collapse=", ") )
            }
        }else if(inherits(data[,cn], "matrix")){
            if(length(unique(data[,x])) > 2){
                out <- sprintf("a matrix predictor; with %d values ranging from %f to %f.", 
                    length(unique(data[,x])),
                    min(data[,x], na.rm=TRUE), max(data[,x], na.rm=TRUE))
            }else{
                out <- sprintf("matrix predictor; set to the value(s): %s.", paste(unique(data[,x]), collapse=", ") ) 
            }
        }else{
            vals <- sort(unique(data[,x]))
            n.cur <- length(vals)+1
            if(!is.null(n)){
                n.cur <- n
            }
            if(length(vals)>n.cur){
                out <- sprintf("%s vector with %d values; set to the value(s): %s, ...", 
                    class(data[,cn])[1],
                    length(vals),
                    paste( vals[1:n.cur], collapse=", ") )
            }else{
                out <- sprintf("%s vector; set to the value(s): %s.", 
                    class(data[,cn])[1],
                    paste( vals, collapse=", ") )
            }
        }
        return(out)
    }

    mysummary <- sapply(colnames(data), function(x){labelColumns(x, data)})
    if(print){
        print_summary(mysummary)
    }  
    invisible(mysummary)
}


#' Utility function.
#' 
#' @export
#' @param sumlist Named list, output of \code{\link{summary_data}}.
#' @param title Optional, text string that will be printed as title.
#' @author Jacolien van Rij
#' @family Data utility functions

print_summary <- function(sumlist, title=NULL){
    if(is.null(title)){
        cat("Summary:\n")
    }else{
        cat(title, "\n")
    }
    for(x in names(sumlist)){
        cat("\t*", x, ":", sumlist[[x]], "\n")
    }
}


#' Utility function.
#' 
#' @export
#' @description Function uses \code{\link[base]{sort.list}} to return indices
#' of of a vector, sorted per group.
#' @param x A vector to be sorted.
#' @param group A names list that specify the different groups to split the 
#' data.
#' @param decreasing Logical: whether or not the sort order should be 
#' decreasing.
#' @return Indices indicating the order of vector x per group.
#' @author Jacolien van Rij
#' @examples
#' data(simdat)
#' range(simdat$Y)
#' ind <- group_sort(simdat$Y, 
#'     group=list(Group=simdat$Group, Trial=simdat$Trial))
#' head(simdat[ind,])
#' @seealso \code{\link[base]{sort.list}}
#' @family Data utility functions

group_sort <- function(x, group=NULL, decreasing=FALSE){
    if(is.null(group)){
        return(sort.list(as.numeric(x), decreasing=decreasing))
    }else{
        el <- which(is.na(x))
        tmp <- data.frame(x=x, i=1:length(x))
        x.split <- split(tmp, f=group, drop=TRUE)
        x.split <- as.vector(unlist(lapply(x.split, 
            function(x){
                x[sort.list(as.numeric(x$x), decreasing=decreasing),'i']
            })))
        if(length(el) > 0){
            x.split <- c(x.split, el[!el %in% x.split])
        }
        return(x.split)
    }
}

#' Defining timebins.
#' 
#' @export
#' @description Function for calculating timebins.
#' @param x Numerical vector with timestamp information.
#' @param binsize Size of the timebin, measured in the same units (often ms) 
#' as \code{x}.
#' @param pos Numerical value that determines the label of the binsize 
#' as proportion of the binsize. A value of 0 will provide the minimum 
#' timestamp within the bin as label, a value of 1 will provide the maximum 
#' value within the bin as label. Defaults to 0.5, the center of the bin.
#' @return Anumerical vector of the same size as \code{x} with timebin 
#' information.
#' @author Jacolien van Rij
#' @family Data utility functions
#' @examples
#' data(simdat)
#' # grouping Time values in bins:
#' simdat$Timebin <- timeBins(simdat$Time, 200)
#' head(simdat)
#' 
#' # different labels:
#' simdat$Timebin2 <- timeBins(simdat$Time, 200, pos=0)
#' head(simdat)

timeBins <- function(x, binsize, pos=.5){
  return( ( floor( x / binsize )+ pos ) * binsize ) 
}

#' Utility function.
#' 
#' @export
#' @description Function to change cbind columns into separated columns.
#' @param x A cbind-vector to be split in columns.
#' @param data (Optional) the data frame.
#' @param values A vector with new column names.
#' @return A data.frame with 
#' @author Jacolien van Rij
#' @seealso \code{\link[base]{table}}
#' @family Functions for binomial count data.
#' @section Note: 
#' Will be moved to other package.
#' # simulate some gaze data:
#' dat <- data.frame(
#'  Subject = rep(1:3, 500),
#'  Timestamp = rep(1:500, 3),
#'  AOI = rep( rep( c('other','competitor', 'target'), 3), 
#'  c(184, 172, 144, 51, 264, 185, 127, 2, 371)) )
#' # add missing data:
#' dat[sample(nrow(dat), size=15),]$AOI <- NA
#'
#' # add timebins:
#' dat$Timebin <- timeBins(dat$Timestamp, 100)
#'
#' # calculate counts:
#' c1 <- getCountData('AOI', data=dat, split_by=c('Subject', 'Timebin'))
#' head(c1)
#'
#' # make columns:
#' c1 <- cbindToColumn("AOI", data=c1)
#' head(c1)

cbindToColumn <- function(x, data=NULL, values=NULL){

    colname <- NULL
    if(is.vector(x)){
        if(is.null(data)){
            stop("If x is a column name, provide data.")
        }else{
            colname <- x[1]
            if(length(x)>1){
                warning("Only first element of x will be used.")
            }
            colname <- x[1]
            x <- data[,colname]
            data[,colname] <- NULL
        }
    }

    if(is.matrix(x)){
        if(!is.null(values)){
            if(length(values)< ncol(x)){
                warning("Length of 'values' less than number of columns in 'x'.")
                values <- c(values, paste("NewColumn", 1:(ncol(x)-length(values)), sep='.'))
            }else if(length(values)> ncol(x)){
                warning(sprintf( "Length of 'values' more than number of columns in 'x'. Only %.0f columns will be used.",
                    ncol(x)) )
                values <- values[1:ncols(x)]
            }
        }else{
            if(!is.null(colname)){
                values <- paste(colname, colnames(x), sep=".")
            }else{
                values <- colnames(x)
            }
        }

        newdat <- as.data.frame(x)
        if(is.null(values)){
            colnames(newdat) <- paste("count", 1:ncol(x), sep='.')
        }else{
            colnames(newdat) <- values
        }

        if(!is.null(data)){
            if(nrow(data)==nrow(newdat)){
                return(cbind(data, newdat))
            }
        }

        return(newdat)
    }else {
        stop("Format of x is not matrix.")
    }
}


#' Utility function.
#' 
#' @export
#' @description Function uses \code{\link[base]{table}} with factors
#' to count the occurrances of each value in the predictor.
#' @param x A vector to be counted.
#' @param values A vector with all possible group names.
#' @param incl_na Logical: whether or not to return a count of missing values.
#' @return A table with count information.
#' @section Note:
#' Values that are not specified in \code{values} will be ignored.
#' @author Jacolien van Rij
#' @seealso \code{\link[base]{table}}
#' @family Functions for binomial count data.
#' @section Note: 
#' Will be moved to other package. 

countValues <- function(x, values, incl_na=FALSE){
    counts <- factor(x, levels=values)
    return( table(counts, 
        useNA=ifelse(incl_na==TRUE, "always", "no"), 
        dnn=list(deparse(substitute(x)))) )
    
}



#' Reducing fixations to count data.
#' 
#' @export
#' @import stats
#' @description Function uses \code{\link[base]{table}} with factors
#' to count the occurrances of each value in the predictor.
#' @param x Name of a column in \code{data} or vector with values 
#' to be counted.
#' @param data Optional: data frame.
#' @param split_by Vector with column names in \code{data} or named list with 
#' grouping predictors.
#' @param values Values of \code{x} that should be counted. 
#' If NULL (default) all values are counted.
#' @param incl_na Logical: whether or not to return a count of missing values.
#' @return A data frame with cbinded count information.
#' @section Note:
#' Values that are not specified in \code{values} will be ignored.
#' @seealso \code{\link[base]{table}}
#' @family Functions for gaze data
#' @author Jacolien van Rij
#' @section Note: 
#' Will be moved to other package. 
#' @examples
#' # simulate some gaze data:
#' dat <- data.frame(
#'  Subject = rep(1:3, 500),
#'  Timestamp = rep(1:500, 3),
#'  AOI = rep( rep( c('other','competitor', 'target'), 3), 
#'  c(184, 172, 144, 51, 264, 185, 127, 2, 371)) )
#' # add missing data:
#' dat[sample(nrow(dat), size=15),]$AOI <- NA
#'
#' # add timebins:
#' dat$Timebin <- timeBins(dat$Timestamp, 100)
#'
#' # calculate counts:
#' c1 <- getCountData('AOI', data=dat, split_by=c('Subject', 'Timebin'))
#' head(c1)
#' # calculating proportions:
#' c1$prop <- c1$AOI[,'target'] / ( c1$AOI[,'competitor']+c1$AOI[,'other'])
#'
#' # calculating counts for specific values, including missing data.
#' # Note that 'distractor' is not in the data:
#' c2 <- getCountData('AOI', data=dat, split_by=c('Subject', 'Timebin'),
#' values=c('target', 'distractor', 'competitor', 'other'), incl_na=TRUE)
#' head(c2)

getCountData <- function(x, split_by, data=NULL, values=NULL, 
    incl_na=FALSE){
    
    dat <- c()
    if(is.null(data)){
        dat <- data.frame(x=x)
        if(is.null( names(split_by) ) ){
            warning('Split_by is a not a named list. Names are given automatically.')
            names(split_by) <- paste('V', 1:length(split_by), sep='')
        }
        for(i in names(split_by)){
            dat[,i] <- split_by[[i]]
        }
    }else{
        split_by <- unlist(split_by)
        cn <- colnames(data)
        if(! x %in% cn){
            stop(sprintf('Column %s not found in data %s.', 
                deparse(substitute(x)),
                deparse(substitute(data))))
        }
        el <- which(!split_by %in% cn)
        if( length(el)> 0){
            stop(sprintf('Column(s) %s not found in data %s.', 
                paste(split_by[el], collapse=', '),
                deparse(substitute(data))))
        }
        dat <- data.frame( data[,c(x, split_by)], row.names=NULL)
        names(dat)[names(dat)==x] <- 'x'
    } 

    if(is.null(values)){
        values <- sort(unique(dat$x))
    }

    out <- values
    if(incl_na){
        out <- c(out, NA)
    }

    split <- list()
    if(length(split_by)>1){
        split <-  as.list(dat[,colnames(dat)!='x'])
    }else{
        split <-  list(dat[,colnames(dat)!='x'])
    }

    # To avoid dropping levels, treat NA in predictors as category:
    split <- lapply(split, function(x){factor(as.character(x), exclude="")})

    ## make more efficient with dplyr or data.table
    newd <- aggregate(dat$x, by=split, 
        function(x){unlist(countValues(x,values=values, incl_na=incl_na))})

    ## Convert numeric split predictors back to numeric:
    for(i in names(split)){
        if(is.numeric(dat[,i])){
            newd[,i] <- as.numeric(as.character(newd[,i]))
        }
    }
    
    colnames(newd)[colnames(newd)=='x'] <- deparse(substitute(x))
    colnames(newd) <- gsub('"','', colnames(newd))
    return(newd)
}
