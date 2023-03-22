#' FindClustersForUniroot
#' @param x numeric value for resolution setting
#' @param object A \code{Seurat} object.
#' @param nCluster number of clusters
#' @param graph.name Name of graph to use for the clustering algorithm
#' @param ... Arguments passed to \code{\link{FindClusters}}
#' @importFrom Seurat FindClusters
#' @return an integer number
#' @export
#'
FindClustersForUniroot <- function(x,
                                   object,
                                   nCluster=6,
                                   graph.name=NULL,
                                   ...) {
  print(x)
  object <- FindClusters(object,
                         resolution = x/100,
                         graph.name = graph.name, ...)
  currentN <- length(unique(object$seurat_clusters))
  return(currentN - nCluster)
}

#' FindNClusters
#' @param object A \code{Seurat} object.
#' @param nCluster number of clusters
#' @param graph.name Name of graph to use for the clustering algorithm
#' @param resolution numeric value for resolution setting
#' @param resolutionMin numeric value for minimum resolution setting
#' @param resolutionMax numeric value for maximum resolution setting
#' @param verbose logical value for verbose message. FALSE as default
#' @param ... arguments passed to FindClustersForUniroot or FindClusters
#' @importFrom ssanv uniroot.integer
#' @importFrom Seurat FindClusters
#' @return a seurat object
#' @export
#'
FindNClusters <- function(object,
                          nCluster = 6,
                          graph.name = NULL,
                          resolution = 0.8,
                          resolutionMin = max(0.2, resolution - 1),
                          resolutionMax = min(resolution + 1, 3),
                          verbose = FALSE,
                          ...) {
  resolutionRange <- c(resolutionMin, resolutionMax)

  #check max or min nCluster in resolution range
  resultMin <-
    FindClustersForUniroot(object = object,
                           nCluster = nCluster,
                           x = resolutionMin*100,
                           graph.name = graph.name, ...)
  resultMax <-
    FindClustersForUniroot(object = object,
                           nCluster = nCluster,
                           x = resolutionMax*100,
                           graph.name = graph.name, ...)

  if (resultMin >= 0) { #smallest N equal or still larger than nCluster
    object <- FindClusters(object,
                           resolution = resolutionMin,
                           graph.name = graph.name, ...)
    if (resultMin > 0) {
      warning(paste0("Get nCluster=", resultMin + nCluster,
                     " with smallest resolution=", resolutionMin,
                     ". Can't get nCluster=", nCluster))
    }
    return(object)
  } else if (resultMax <= 0) { #largest N still less than nCluster
    object <- FindClusters(object,
                           resolution = resolutionMax,
                           graph.name = graph.name, ...)
    if (resultMax < 0) {
      warning(paste0("Get nCluster=", resultMax + nCluster,
                     " with largest resolution=", resolutionMax,
                     ". Can't get nCluster=", nCluster))
    }
    return(object)
  }
  result <- ssanv::uniroot.integer(
    f = FindClustersForUniroot,
    interval = resolutionRange*100,
    step.power = 5,
    object=object,nCluster = nCluster,
    graph.name = graph.name,
    ...)
  selectedResolution <- result$root/100
  object <- FindClusters(object,
                         resolution = selectedResolution,
                         graph.name = graph.name, ...)
  return(object)
}

