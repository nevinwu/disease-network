################################################################################
## function_optimal_nclusters.r
## 15/03/2021
## FIVI
##
## Function to plot a guess optimal cluster number in a clustering.
################################################################################

#' Guesses optimal cluster number in a clustering
#'
#'\code{optimal_nclusters} optimal_nclusters
#'
#'@description Function to plot a guess optimal cluster number in a clustering
#'@param dat a expression dataset
#'@return Guesses optimal cluster number in a clustering
#'@import
#'@export

optimal_nclusters <- function (dat,
                                method = "silhouette",
                                max_nclusters = 10,
                                title = "")
{
  method = match.arg(method, choices = c("elbow", "silhouette"))
  if (max_nclusters > ncol(dat)) {
    max_nclusters = ncol(dat) - 1
  }
  if (method == "silhouette") {
    sil <- rep(0, 8)
    for (i in 2:10) {
      km.res <- kmeans(t(dat), centers = i, nstart = 25)
      ss <- silhouette(km.res$cluster, dist(t(dat)))
      sil[i] <- mean(ss[, 3])
    }
    silPlot = data.frame(sil = sil, clusters = 1:10)
    p = ggplot(silPlot,
               aes(x = factor(clusters),
                   y = sil,
                   fill = "1")) +
        geom_point(shape = 1, size = 4) +
        geom_line(data = silPlot, aes(x = clusters, y = sil),
                  size = 1, colour = "azure4") +
        geom_vline(xintercept = silPlot$clusters[which.max(silPlot$sil)],
                   linetype = 2, color = "steelblue") +
        theme_bw() +
        xlab("Number of clusters (k)") +
        ylab("Average silhouette coefficient") +
        ggtitle(title) +
        theme(legend.position = "none",
              plot.title = element_text(size = 22, hjust = 0.5),
              axis.text = element_text(size = 15, color = "dimgray"),
              axis.title = element_text(size = 18, hjust = 0.5))
  }
  else {
    wss <- sapply(1:10, function(k) {
      kmeans(t(dat), k, nstart = 10)$tot.withinss
    })
    screePlot = data.frame(wss = wss, clusters = 1:10)
    p = ggplot(screePlot,
               aes(x = factor(clusters),
                   y = wss, fill = "1")) +
        geom_point(shape = 1, size = 4) +
        geom_line(data = screePlot, aes(x = clusters, y = wss),
                  size = 1, colour = "azure4") +
        theme_bw() +
        xlab("Number of clusters (k)") +
        ylab("Total within groups sum of squares") +
        ggtitle(title) +
        theme(legend.position = "none",
              plot.title = element_text(size = 22),
              axis.text = element_text(size = 15, color = "dimgray"),
              axis.title = element_text(size = 18, hjust = 0.5))
  }
  return(p)
}
