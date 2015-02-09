#' Determine the starting point for each time series.
#' 
#' @export
#' @param data A data frame.
#' @param column Test string, name of the column that describes the order 
#' withing the time series. Default is "Time".
#' @param event A text string or vector indicating the columns that define the 
#' unique time series. Default is "Event".
#' @param label The name of the new column with the start point of each time 
#' series. Default is "start.event".
#' @param label.event In case \code{event} is not a single column, providing a 
#' text string will add a column with this name that defines unique time 
#' series. Default is NULL (no new column for time series is created).
#' @param order Logical: whether or not to order each time series. 
#' Default is TRUE, maybe set to FALSE with large data frames that are already ordered.
#' @return Data frame.
#' @author Jacolien van Rij
#' @examples 
#' data(simdat)
#' head(simdat)
#' test <- start_event(simdat, event=c("Subject", "Trial"), label.event="Event") 
#' head(test)

start_event <- function(data, column="Time", event="Event", label="start.event", label.event=NULL, order=TRUE){
	if(is.null(label)){
		stop("No output column specified in argument 'label'.")
	}

	if(!column %in% names(data)){
		stop(sprintf("No column '%s' in data frame %s.", column, deparse(substitute(data))))
	}
	if(!all(event %in% names(data))){
		el <- paste( which(event[!event %in% names(data)]), collapse=',' )
		stop(sprintf("Column name(s) '%s' not found in data frame %s.", el, deparse(substitute(data))))
	}
	if(length(column) > 1){
		warning(sprintf("Argument column has %d elements. Only first element is used.", length(column)))
		column <- column[1]
	}
	if(label %in% names(data)){
		warning(sprintf("Column %s already exists, will be overwritten.", ))
	}
	if(!is.null(label.event)){
		if(label.event %in% names(data)){
			warning(sprintf("Column %s already exists, will be overwritten.", ))
		}
	}
	tmp <- NULL
	if(!is.null(label.event) & length(event) > 1){
		data[, label.event] <- interaction(data[,event])
		tmp <- split(data, f=list(data[,label.event]), drop=TRUE)
	} else if (length(event) > 1){
		tmp <- split(data, f=as.list(data[,event]), drop=TRUE)
	} else {
		tmp <- split(data, f=list(data[,event]), drop=TRUE)
	}

	tmp <- lapply(tmp, function(x){
		min.x <- min(x[,column])
		x[,label] <- x[,column]==min.x
		if(order){
			x <- x[order(x[,column]),]
		}
		return(x)
	})
	tmp <- do.call('rbind', tmp)
	row.names(tmp) <- NULL
	return(tmp)
}