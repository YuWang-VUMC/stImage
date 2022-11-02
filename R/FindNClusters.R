#' Title
#'
#' @param x
#' @param object
#' @param nCluster
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
FindClustersForUniroot <- function(x,
                                   object,
                                   nCluster=6,
                                   graph.name=NULL,
                                   ...) {
  print(x)
  #browser()
  object <- FindClusters(object, resolution=x/100, graph.name=graph.name, ...)
  currentN <- length(unique(object$seurat_clusters))
  return(currentN - nCluster)
}

#' Title
#'
#' @param object
#' @param nCluster
#' @param resolution
#' @param resolutionMin
#' @param resolutionMax
#' @param verbose
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
FindNClusters <- function(object,
                          nCluster = 6,
                          graph.name = NULL,
                          resolution = 0.8,
                          resolutionMin = max(0.2, resolution - 1),
                          resolutionMax = min(resolution + 1, 3),
                          verbose = FALSE,
                          ...) {
  #browser()
  resolutionRange <- c(resolutionMin, resolutionMax)

  #check max or min nCluster in resolution range
  resultMin <- FindClustersForUniroot(object = object, nCluster = nCluster,
                                      x = resolutionMin*100,graph.name=graph.name, ...)
  resultMax <- FindClustersForUniroot(object = object, nCluster = nCluster,
                                      x = resolutionMax*100,graph.name=graph.name, ...)
  #browser()

  if (resultMin >= 0) { #smallest N equal or still larger than nCluster
    object <- FindClusters(object, resolution = resolutionMin,graph.name=graph.name,...)
    if (resultMin > 0) {
      warning(paste0("Get nCluster=", resultMin + nCluster, " with smallest resolution=", resolutionMin, ". Can't get nCluster=", nCluster))
    }
    return(object)
  } else if (resultMax <= 0) { #largest N still less than nCluster
    object <- FindClusters(object, resolution = resolutionMax,graph.name=graph.name, ...)
    if (resultMax < 0) {
      warning(paste0("Get nCluster=", resultMax+nCluster, " with largest resolution=", resolutionMax, ". Can't get nCluster=", nCluster))
    }
    return(object)
  }
  result=ssanv::uniroot.integer(f=FindClustersForUniroot,interval=resolutionRange*100,
                                step.power=5,
                                object=object,nCluster=nCluster,
                                graph.name=graph.name,
                                ...)
  selectedResolution=result$root/100
  #browser()
  object <- FindClusters(object, resolution = selectedResolution, graph.name=graph.name, ...)
  return(object)
}

#temp=FindNClusters(object, nCluster=nCluster,graph.name = "ImageFeaturePCA_snn")




