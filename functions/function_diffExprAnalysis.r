################################################################################
## function_diffExprAnalysis.r
## 15/03/2021
## FIVI
##
## Function to perform a differential expression analysis.
################################################################################

#' Performs a differential expression analysis
#'
#'\code{diffExprAnalysis} diffExprAnalysis
#'
#'@description Function to perform a differential expression analysis
#'@param dat a expression dataset
#'@return Performs a differential expression analysis
#'@import
#'@export

diffExprAnalysis <- function (dat,
                              ed,
                              condition,
                              paired = FALSE,
                              pair = NULL)
{
  group = ed[colnames(dat), condition]
  if (paired) {
    pair = ed[colnames(dat), pair]
    design <- model.matrix(~0 + group + pair)
    rownames(design) <- colnames(dat)
    fit <- lmFit(dat, design = design)
    res.limma <- eBayes(fit)
    res = list(topTable(res.limma, coef = 1, number = nrow((dat))))
    names(res) = paste(gsub("group", "", colnames(design)[1:2]),
                       collapse = "-")
  }
  else {
    design <- model.matrix(~0 + group)
    rownames(design) <- colnames(dat)
    combs = combn(x = colnames(design), m = 2, simplify = TRUE)
    cons2 = apply(combs, 2, function(y) {
      paste(y, collapse = "-")
    })
    cons = gsub("group", "", cons2)
    contrasts <- makeContrasts(contrasts = cons2, levels = design)
    fit <- lmFit(dat, design = design)
    fit.cont <- contrasts.fit(fit, contrasts)
    res.limma <- eBayes(fit.cont)
    res = list()
    for (i in 1:length(cons2)) {
      res[[cons2[i]]] = topTable(res.limma, coef = i, number = nrow(dat))
    }
    names(res) = gsub("group", "", names(res))
  }
  res2 = lapply(names(res), function(x) {
    aux = unlist(strsplit(x, split = "-", fixed = T))
    condition1 = aux[1]
    condition2 = aux[2]
    res[[x]]$FC = unlist(lapply(rownames(res[[x]]), function(y) {
      A = mean(as.numeric(dat[y, group == condition1]))
      B = mean(as.numeric(dat[y, group == condition2]))
      FC = 2^(A - B)
      if (FC < 1) {
        FC = -1/FC
      }
      FC
    }))
    res[[x]]
  })
  names(res2) = names(res)
  return(res2)
}
