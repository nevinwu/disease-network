################################################################################
## function_arbol.r
## 10/03/2010
## dmontaner@cipf.es
##
## To be called after "hclust"
## $clase ELEMENT has to be introduced in the result of "hclust"
##
## EXAMPLE:
## distancia <- dist (t (datos), method = "euclidean") #dist works in rows
## he <- hclust (distancia)
## he$clase <- batch  ## COLOUR ACCORNG TO batch
## arbol (cluster = he, main = "euclidean distance")
################################################################################

#' Function to plot a clustering tree
#'
#'\code{arbol} Function to plot a clustering tree
#'
#'@description Function to plot a clustering tree
#'@param cluster result of hclust() function
#'@return Plots a clustering tree
#'@import
#'@export

arbol <- function (cluster, ...) {
  ##
  plot (cluster, hang = 0.1, axes = FALSE, ann = FALSE)
  title (...)
  ##
  ##color
  clase <- as.character (cluster$clase)
  clase.unica <- unique (clase)
  clr <- rainbow (length (clase.unica))
  names (clr) <- clase.unica
  ##
  posiciones <- 1:length (cluster$labels)
  labels.ordenados <- cluster$labels[cluster$order]
  clase.ordenada <- clase[cluster$order]
  ##
  for (cls in clase.unica) {
    touse <- clase.ordenada %in% cls
    axis (1, posiciones[touse], labels = labels.ordenados[touse],
          las = 3, col.axis = clr[cls])
  }
  ##
  try (legend (x = 3, y = max (cluster$height), legend = names (clr), fill = clr))
  #try (legend (x = "topright", legend = names (clr), fill = clr))
}
