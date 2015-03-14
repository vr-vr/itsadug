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


#' Utility function.
#' 
#' @description Function uses \code{\link[base]{table}} with factors
#' to count the occurrances of each value in the predictor.
#' @param x A vector to be counted.
#' @param values A vector with all possible group names.
#' @param incl_na Logical: whether or not to return a count of missing values.
#' @return A table with count information.
#' @section Note:
#' Values that are not specified in \code{values} will be ignored.
#' @seealso \code{\link[base]{table}}
#' @family Functions for binomial count data.

countValues <- function(x, values, incl_na=FALSE){
    counts <- factor(x, levels=values)
    return( table(counts, 
        useNA=ifelse(incl_na==TRUE, "always", "no"), 
        dnn=list(deparse(substitute(x)))) )
    
}



#' Utility function.
#' 
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
#' @family Functions for binomial count data.

getCountData <- function(x, split_by, data=NULL, values=NULL, incl_na=FALSE){
    
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

#' Utility function.
#' 
#' @description Create cbind of successes and failures.
#'
#' @param success Vector with column(s) of count data 
#' defining the number of successes.
#' @param failure Vector with column(s) of count data 
#' defining the number of failures.
#' @param data Data frame.
#' @param name Text string, name of column with the cbind of 
#' successes and failures. Only if \code{data} is specified.
#' @return cbind count data, which sums up to 2 for each row.
#' @author Jacolien van Rij
#' @family Functions for binomial count data.

convertToBinom <- function(success, failure, data=NULL, name="cbind"){
    # check colnames:
    sval <- c()
    fval <- c()
    if(is.list(success)){
        if(length(success)==1){
            sval <- success[[1]]
        }else{
            sval <- rowSums(as.data.frame(success))
        }
    }else if(!is.null(dim(success))){
        if(dim(success)[2]==1){
            sval <- success[,1]
        }else{
            sval <- rowSums(success)
        }
    }else{
        sval <- success
    }
    if(is.list(failure)){
        if(length(failure)==1){
            fval <- failure[[1]]
        }else{
            fval <- rowSums(as.data.frame(failure))
        }
    }else if(!is.null(dim(failure))){
        if(dim(failure)[2]==1){
            fval <- failure[,1]
        }else{
            fval <- rowSums(failure)
        }
    }else{
        fval <- failure
    }

    if(is.null(data)){
        return(t( mapply(function(x,y){cbind(x,y)},
            sval, fval) ))
    }else{
        data[,name] <- t( mapply(function(x,y){cbind(x,y)},
            sval, fval) ) 
        return(data)
    }
}


#' Reducing cbind binomial count data to avoid autocorrelation issues.
#'
#' @param x cbind data column (leave \code{y} NULL), or a vector 
#' specify \code{y}).
#' @param y Vector, in case \code{x} is vector, otherwise NULL.
#' @return cbind count data, which sums up to 2 for each row.
#' @author Jacolien van Rij
#' @family Functions for binomial count data.

reduceCounts <- function(x,y=NULL){
  check <- function(a,b){
    if(any(is.na(c(a,b)))){
      return(NA)
    }else if(a==b){
        if(a==0){
            return(cbind(0,0))
        }else{
            return(cbind(1,1))
        }
    }else if(a > b){
      return(cbind(2,0))
    }else if(a < b){
      return(cbind(0,2))
    }else{
      return(NA)
    }
  }  
  if(is.null(y)){
    return( t(mapply(function(x,y){check(x,y)}, x[,1], x[,2])) )
  }else{
    return( t(mapply(function(x,y){check(x,y)}, x, y)) )
  }
}
