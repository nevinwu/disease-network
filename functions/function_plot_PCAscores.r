################################################################################
## function_plot_PCAscores.r
## 20/02/2021
## FIVI
##
## Function to plot a PCA results.
################################################################################

#' Plots PCA
#'
#'\code{plot_PCAscores} Plots PCA
#'
#'@description Function to plot a PCA results. Obtained from the FIVI
#'bioinformatics group-
#'@param dat a expression dataset
#'@return plots PCA results
#'@import
#'@export

plot_PCAscores <- function (dat,
                            ed,
                            components = c(1, 2),
                            condition1,
                            colors = NULL,
                            condition2 = NULL,
                            shapes = NULL,
                            title = "",
                            legendpos = "bottom",
                            labels = NULL, ellipses = F,
                            ellipse_level = 0.95,
                            ellipse_alpha = 0.2)

{
  pca = prcomp(t(dat))
  Var = round(summary(pca)$importance[2, components] * 100,
              1)
  toplot = data.frame(pca$x[, components], stringsAsFactors = F)
  lim = max(abs(c(min(toplot), max(toplot))))
  if (ellipses) {
    lim = lim + 0.8 * lim
  }
  axis_limits = c(-lim, lim)
  toplot$color = ed[colnames(dat), condition1]
  cbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442",
                 "#0072B2", "#D55E00", "#CC79A7", "#9370DB", "#B3EE3A",
                 "#FFC1C1", "#999999", "#000000")
  if (is.null(colors)) {
    colors = cbPalette[1:length(unique(ed[, condition1]))]
  }
  if (!is.null(condition2)) {
    toplot$shape = ed[colnames(dat), condition2]
    if (is.null(shapes)) {
      shapes = 15:(length(unique(ed[, condition2])) + 15)
    }
  }
  if (!is.null(labels)) {
    toplot$labels = ed[colnames(dat), labels]
  }
  p = ggplot(toplot, aes_string(x = colnames(toplot)[1],
                                y = colnames(toplot)[2]))
  if (!is.null(condition2)) {
    p = p + geom_point(aes(color = color, shape = shape),  size = 3) +
      scale_color_manual(name = condition1, values = colors) +
      scale_shape_manual(name = condition2, values = shapes)
  }
  else {
    p = p + geom_point(aes(color = color), size = 3) +
      scale_color_manual(name = condition1, values = colors)
  }
  if (!is.null(labels)) {
    p = p + geom_text_repel(aes(label = labels, color = toplot$color),
                            size = 4, point.padding = unit(0.2, "lines"),
                            segment.size = 0)
  }
  if (ellipses) {
    p = p + stat_ellipse(aes(fill = color), geom = "polygon",
                         level = ellipse_level, alpha = ellipse_alpha,
                         show.legend = FALSE) +
      scale_fill_manual(values = colors)
  }
  p = p + xlab(paste0("PC", components[1], ": ", Var[1], "%")) +
    ylab(paste0("PC", components[2], ": ", Var[2], "%")) +
    ggtitle(title)
  p = p + xlim(axis_limits) + ylim(axis_limits)
  p = p + theme_light() + theme(legend.position = legendpos,
                                axis.title = element_text(size = 18),
                                axis.text = element_text(size = 15),
                                plot.title = element_text(size = 22,
                                                          hjust = 0.5),
                                legend.title = element_blank(),
                                legend.text = element_text(size = 13))
  return(p)
}
