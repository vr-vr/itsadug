.onAttach <- function(...) {
  packageStartupMessage('Loaded package itsadug 0.4 (see \'help(package="itsadug")\' ).')
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
		cat(sprintf("Package: itsadug\nVersion: %s\n", 
			packageVersion("itsadug")))
	}else{
		help(package="itsadug")
	}
	
}