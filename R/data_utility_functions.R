#' Utility function.
#' 
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
#' @description The function prints a summary of the data. 
#' Similar to the function \code{\link[utils]{str}}, but easier readable.
#' @param data A data frame.
#' @author Jacolien van Rij
#' @examples
#' data(simdat)
#' summary_data(simdat)
#' @family Data utility functions

summary_data <- function(data){
    cat("Summary:\n")
    
    labelColumns <- function(x, data){
        out <- NULL
        cn <- ifelse(is.numeric(x), colnames(data)[x], x)
        cl <- class(data[,cn])
        if(inherits(data[,cn],'factor')){
            cat( sprintf("\t* %s is a factor; set to the value(s): %s.\n", cn, 
                paste( sort(unique(as.character(data[,x]))), collapse=",") ) )
        }else if(inherits(data[,cn],'numeric')){
            if(all(data[,x] %in% c(0,1,NA))){
                 cat( sprintf("\t* %s is a binomial predictor; set to the value(s): %s.\n", cn, paste(unique(data[,x]), collapse=",") ) )
            }else{
                if(length(unique(data[,x])) > 2){
                    cat( sprintf("\t* %s is a numeric predictor; with %d values ranging from %f to %f.\n", 
                        cn, length(unique(data[,x])),
                        min(data[,x], na.rm=TRUE), max(data[,x], na.rm=TRUE)) )
                }else{
                    cat( sprintf("\t* %s is a numeric predictor; set to the value(s): %s.\n", cn, paste(unique(data[,x]), collapse=",") ) ) 
                }
            }
        }else if(inherits(data[,cn], "matrix")){
            if(length(unique(data[,x])) > 2){
                cat( sprintf("\t* %s is a matrix predictor; with %d values ranging from %f to %f.\n", 
                    cn, length(unique(data[,x])),
                    min(data[,x], na.rm=TRUE), max(data[,x], na.rm=TRUE)) )
            }else{
                cat( sprintf("\t* %s is a matrix predictor; set to the value(s): %s.\n", cn, paste(unique(data[,x]), collapse=",") ) ) 
            }
        }else{
            cat( sprintf("\t* %s is a %s vector; set to the value(s): %s.\n", 
                cn, class(data[,cn])[1],
                paste( sort(unique(data[,x])), collapse=",")) ) 
        }
        return(out)
    }

    test <- sapply(colnames(data), function(x){labelColumns(x, data)})
}



#' Utility function.
#' 
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
