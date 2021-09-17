################################################################################
## function_pca.r
## 23/11/2009
## dmontaner@cipf.es
##
## Function to plot pca from arrays. Modified version of Paco's idea.
## Check and put dimnames in results' matrix
################################################################################

#' Plots PCA from array
#'
#'\code{pcaGenes} Plots PCA from array
#'
#'@description A very useful function to obtain principal components analysis to
#'a matrix that has more variables than individuals.R can not apply princomp in
#'such case and when there are a lot of variables eigen(t(X)%*%X) can not be
#'computed. First we center the matrix by columns (Xoff) and then we obtain the
#'eigenvalues and the eigenvectors of the matrix Xoff%*%t(Xoff) and we use the
#'equivalences between the loadings and scores to obtain the solution.
#'@param X a matrix that has on ROWS the genes considered as variables in the
#'PCA analysis.
#'@return A data object of the same class as x with a row for each unique value
#'of ID.
#'@import
#'@export

pcaGenes <- function (X) {
  # scores1 and loadings1 are what we are looking for, and scores2 and loadings2
  # are scores and loadings of the transposed.

  # Transpose
  X <- t (X)

  # Paco's code
  X <- as.matrix (X)
  n <- ncol(X)
  p <- nrow(X)
  offset <- apply (X, 2, mean)
  Xoff <- X - (cbind (matrix(1, p, 1)) %*% rbind (offset))

  # eigen command sorts the eigenvalues in decreasing orden.
  eigen <- eigen (Xoff %*% t(Xoff) / (p-1))
  var <- cbind (eigen$values / sum(eigen$values),
                cumsum (eigen$values / sum(eigen$values)))

  loadings2 <- eigen$vectors
  scores2 <- t(Xoff) %*% loadings2

  normas2 <- sqrt (apply (scores2^2, 2, sum))

  scores1 <- loadings2 %*% diag (normas2)
  loadings1 <- scores2 %*% diag (1/normas2)

  # Transpose Xoff
  Xoff <- t (Xoff)

  output <- list (eigen, var, scores1, loadings1, Xoff)
  names(output) <- c ("eigen","var.exp","scores","loadings", "Xoff")
  return (output)
}

# To plot results
plot.pca.genes <- function (mi.pca,
                            addNames = TRUE,
                            col = "black",
                            col.class = NULL, ...) {
  par (mfrow = c(1,2))
  ##
  plot (mi.pca$scores[,1], mi.pca$score[,2], col = col, ..., #main = main,
        xlab = paste("PC1: ", round (mi.pca$var.exp[1,1] * 100),
                     "% explained variance", sep = ""),
        ylab = paste("PC2: ", round (mi.pca$var.exp[2,1] * 100),
                     "% explained variance", sep = ""))
  ##legend
  if (!is.null (col.class)) {
    clr <- unique (col)
    names (clr)<- unique (col.class)
    try (legend (x = min (mi.pca$scores[,1]),
                 y = max (mi.pca$score[,2]),
                 legend = names (clr),
                 fill = clr))
    #try (legend (x = 0, y= 200, legend = names (clr), fill = clr))
  }
  ##names
  if (addNames) {
    text (mi.pca$scores[,1],
          mi.pca$score[,2],
          labels = colnames (mi.pca$Xoff),
          col = col, ...)
  }
  ##
  plot (mi.pca$scores[,3], mi.pca$score[,2], col = col, ..., #main = main,
        xlab = paste("PC3: ", round (mi.pca$var.exp[3,1] * 100),
                     "% explained variance", sep = ""),
        ylab = paste("PC2: ", round (mi.pca$var.exp[2,1] * 100),
                     "% explained variance", sep = ""))
  ##legend
  if (!is.null (col.class)) {
    clr <- unique (col)
    names (clr)<- unique (col.class)
    try (legend (x = min (mi.pca$scores[,3]),
                 y = max (mi.pca$score[,2]),legend = names (clr), fill = clr))
    #try (legend (x = -10, y= 10, legend = names (clr), fill = clr))


  }
  ##names
  if (addNames) {
    text (mi.pca$scores[,3], mi.pca$score[,2],
          labels = colnames (mi.pca$Xoff), col = col, ...)
  }
}
