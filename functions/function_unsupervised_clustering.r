################################################################################
## function_unsupervised_clustering.r
## 15/03/2021
## FIVI
##
## Function to plot a guess optimal cluster number in a clustering.
################################################################################

#' Performs an unsupervised clustering
#'
#'\code{unsupervised_clustering} unsupervised_clustering
#'
#'@description Function to perform an unsupervised clustering
#'@param dat a expression dataset
#'@return Performs an unsupervised clustering
#'@import
#'@export

unsupervised_clustering <- function (dat,
                                     method = "kmeans",
                                     nclusters = NULL,
                                     cluster_colors = NULL,
                                     groups = NULL,
                                     group_colors = NULL,
                                     components = c(1, 2),
                                     labels = T,
                                     ellipses = T,
                                     title = "",
                                     legendpos = "bottom",
                                     seed = 123)
{
  set.seed(seed)
  clusterPalette = c("#7FC97F", "#BEAED4", "#FDC086", "#FFFF99",
                     "#386CB0", "#F0027F", "#BF5B17", "#666666")
  groupPalette = c("#E69F00", "#56B4E9", "#009E73", "#F0E442",
                   "#0072B2", "#D55E00", "#CC79A7", "#9370DB", "#B3EE3A",
                   "#FFC1C1", "#999999", "#000000")
  if (is.null(nclusters)) {
    cat("You must especify the optimal number of clusters\n")
    cat("Try before optimal_nclusters() function\n")
  }

  method = match.arg(method, choices = c("kmeans", "pam"))

  if (method == "kmeans") {
    fit = kmeans(t(dat), centers = nclusters)
    clusters = fit$cluster
  }

  else {
    fit = pam(t(dat), k = nclusters)
    clusters = fit$clustering
  }

  if (is.null(cluster_colors)) {
    cluster_colors = clusterPalette[1:nclusters]
  }

  ed_aux = data.frame(sample = colnames(dat), cluster = factor(clusters),
                      stringsAsFactors = F)

  if (labels) {
    show_labels = "sample"
  }

  else {
    show_labels = NULL
  }

  if (!is.null(groups)) {
    ed_aux$group = groups
    if (is.null(group_colors)) {
      group_colors = groupPalette[1:length(unique(groups))]
    }
    p = plot_PCAscores(dat = dat, ed = ed_aux, components = components,
                       condition1 = "group", colors = group_colors,
                       title = title,
                       legendpos = legendpos, labels = show_labels,
                       ellipses = F)
    pca = prcomp(t(dat))
    toplot_ell = data.frame(pca$x[, components], cluster = factor(clusters),
                            stringsAsFactors = F)
    colnames(toplot_ell) = c("X", "Y", "cluster")
    p = p + stat_ellipse(data = toplot_ell, aes(x = X, y = Y,
                                                fill = cluster),
                         geom = "polygon", alpha = 0.2)
  }
  else {
    p = plot_PCAscores(dat = dat, ed = ed_aux, components = components,
                       condition1 = "cluster", colors = cluster_colors,
                       title = title, legendpos = legendpos,
                       labels = show_labels,
                       ellipses = ellipses)
  }
  res = list(method = method, cluster_results = fit, cluster_group = clusters,
             cluster_plot = p)
  return(res)
}
