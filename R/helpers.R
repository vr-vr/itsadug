#' Convert the summary  into a Latex table.
#' 
#' @export
#' @param model A GAM(M) model build in the package \code{\link[mgcv]{mgcv}} 
#' using \code{\link[mgcv]{gam}} or \code{\link[mgcv]{bam}}.
#' @param caption A string with the caption for the table.
#' @param label A string for the label to refer to the table in the markdown 
#' document.
#' @param pnames A vector with labels to relabel the rows in the parametric 
#' part of the summary.
#' @param snames A vector with labels to relabel the rows in the smooth
#' part of the summary.
#' @param ptab A vector with labels to relabel the column names of the 
#' parametric summary.
#' @param stab A vector with labels to relabel the column names of the 
#' smooth summary.
#' @param ... Optional additional arguments which are passed to 
#' \code{xtable} (see 'help(xtable)').
#' @return A vector with color values.
#' @author R. Harald Baayen
#' @section Note: 
#' This function is useful for markdown documents using the package 
#' \code{knitr} to integrate R code with Latex and Sweave. This 
#' function requires the package \code{xtable}.
#' @seealso \code{\link[mgcv]{summary.gam}}, \code{\link[mgcv]{gam}}, 
#' \code{\link[mgcv]{bam}}.
# help function

gamtabs = function(model, caption=" ", label="tab.gam", pnames = NA, snames = NA, 
                   ptab=NA, stab=NA, ...) 
{
  if (!requireNamespace("xtable", quietly = TRUE)) {
    stop("Package 'xtable' needed for this function to work. Please install it.",
      call. = FALSE)
  }
  if (is.na(ptab[1])) ptab = as.data.frame(summary(model)$p.table)
  if (is.na(stab[1])) stab = as.data.frame(summary(model)$s.table)
  if (!is.na(pnames[1])) rownames(ptab)=pnames
  if (!is.na(snames[1])) rownames(stab)=snames
  colnames(ptab)[4]="p-value"
  colnames(ptab)[3]="t-value"
  ptab.cnames = colnames(ptab)
  stab.cnames = colnames(stab)
  stab.cnames[3]="F-value"

  colnames(ptab)=c("A", "B", "C", "D")
  colnames(stab)=colnames(ptab)
  tab = rbind(ptab, stab)
  colnames(tab)=ptab.cnames
  tab = round(tab, 4)
  m = data.frame(matrix(0, nrow(tab), ncol(tab)))
  for (i in 1:nrow(tab)) {
    for (j in 1:4) {
      if ((j==4) & (tab[i,j] < 0.0001)) m[i,j]="< 0.0001"
      else m[i,j]=sprintf("%3.4f", tab[i,j])
    }
  }
  colnames(m) = colnames(tab)
  rownames(m) = rownames(tab)
  tab=m
  tab2 = rbind(c(ptab.cnames),
             tab[1:nrow(ptab),],
             c(stab.cnames),
             tab[(nrow(ptab)+1):nrow(tab),])
  rownames(tab2)[(nrow(ptab)+2)]="B. smooth terms"
  rownames(tab2)[1]="A. parametric coefficients"
  for (i in 1:nrow(tab2)) {
    if (tab2[i,4] == "0") tab2[i,4] = "< 0.0001"
    if (length(grep("\\.", tab2[i,2]))==0) tab2[i,2]=paste(tab2[i,2], ".0000", sep="")
  }
  print(xtable::xtable(tab2, caption=caption, label = label, align="lrrrr"), 
      include.colnames=FALSE, #include.rownames=FALSE,
      hline.after=c(0, (nrow(ptab)+1), nrow(tab2)), ... )
}
