.onAttach <- function(...) {
  if(is.null(getOption('itsadug_print'))){
  	options(itsadug_print=TRUE)
  }
  if(getOption('itsadug_print')){
  	packageStartupMessage('Loaded package itsadug 1.0 (see \'help("itsadug")\' ).')
  }
}

#' Information on how to cite this package
#' 
#' @param input Optional parameter. Normally (NULL) the citation info is 
#' printed. If value "version" then only the version is printed.
#' @examples
#' info()
#' info("version")
#' citation(package="itsadug")
#' # To get info about R version:
#' R.version.string
#' @seealso
#' \code{\link[utils]{citation}}, \code{\link[base]{R.version}},
#' \code{\link[utils]{sessionInfo}}

info <- function(input=NULL){
	if(is.null(input)){
		citation(package="itsadug")
	}else if(input=="version"){
		cat(sprintf("Package itsadug, version %s\n", 
			packageVersion("itsadug")))
	}else{
		help(package="itsadug")
	}
}

#' Turn on or off information messages.
#' 
#' @param input Input variable indicating to print info messages 
#' ("on", or 1, or TRUE) or not ("off", 0, or FALSE).
#' @examples
#' # To turn on the info messages (all the same):
#' infoMessages("on")
#' infoMessages(1)
#' infoMessages(TRUE)
#' # To turn off the info messages (all the same):
#' infoMessages("off")
#' infoMessages(0)
#' infoMessages(FALSE)
#' # checking output:
#' (out <- infoMessages(FALSE))

infoMessages <- function(input){
	if(is.logical(input)){
		options(itsadug_print=input)
	}else if(is.numeric(input)){
		options(itsadug_print=ifelse(input<=0, FALSE, TRUE))
	}else if(is.character(input)){
		options(itsadug_print=ifelse(input=="off", FALSE, 
			ifelse(input=="on",TRUE, getOption('itsadug_print'))))
	}else{
		stop(sprintf("Cannot interpret input value %s. Try to use logical values TRUE or FALSE.", input))
	}
	invisible(list(value=getOption('itsadug_print'), 
		effect=ifelse(getOption('itsadug_print')==TRUE, "messages printed", "no messages printed")))
}

