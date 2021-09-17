################################################################################
## function_medianreps.r
## 01/01/2021
## ipediez93@gmail.com
##
## Used to get median genename for multiple matches genename-probe when
## annotating
################################################################################

#' Median over irregular replicate probes
#'
#'\code{medianReps} Median over irregular replicate probes
#'
#'@description Condense a microarray data objecto so that values for
#'within-array replicate probes are replaced with their median
#'@param matriz a matrix-like object
#'@return A data object of the same class as x with a row for each unique
#'value of ID.
#'@import stats
#'@export

medianReps <- function(matriz){
  ID <- as.character(rownames(matriz))
  ID <- factor(ID, levels = unique(ID))
  df <- by(matriz, ID, function(x) apply(x, 2, stats::median))
  mat <- do.call("rbind", df)
  return(mat)
}
