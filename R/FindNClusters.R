# FindNClusters=function(object,n=6,
#                       resolution=0.8,
#                       resolutionMin = max(0.4,resolution-1),resolutionMax = min(resolution+1,2),resolutionStep=0.1,
#                       verbose=FALSE,
#                       ...) {
#   object=FindClusters(object,resolution=resolution,...)
#
#   while(resolution>=resolutionMin & resolution<=resolutionMax) {
#     currentN=length(unique(object$seurat_clusters))
#     if (currentN==n) {
#       return(object)
#     } else if (currentN>n) { #too many clusters
#       resolution=resolution-resolutionStep
#       if (verbose) {
#         print(paste0("Too few clusters. Decrease resolution to ",resolution))
#       }
#     } else if (currentN<n) { #too few clusters
#       resolution=resolution+resolutionStep
#       if (verbose) {
#         print(paste0("Too many clusters. Increase resolution to ",resolution))
#       }
#     }
#     object=FindClusters(object,resolution=resolution,...)
#   }
#
#   return(object)
#
# }



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
                                   ...) {
  print(x)
  object <- FindClusters(object, resolution=x/100, ...)
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
                          resolution = 0.8,
                          resolutionMin = max(0.2, resolution - 1),
                          resolutionMax = min(resolution + 1, 3),
                          verbose = FALSE,
                          ...) {
  #browser()
  resolutionRange <- c(resolutionMin, resolutionMax)

  #check max or min nCluster in resolution range
  resultMin <- FindClustersForUniroot(object = object, nCluster = nCluster, x = resolutionMin*100, ...)
  resultMax <- FindClustersForUniroot(object = object, nCluster = nCluster, x = resolutionMax*100, ...)
  #browser()

  if (resultMin >= 0) { #smallest N equal or still larger than nCluster
    object <- FindClusters(object, resolution = resolutionMin,...)
    if (resultMin > 0) {
      warning(paste0("Get nCluster=", resultMin + nCluster, " with smallest resolution=", resolutionMin, ". Can't get nCluster=", nCluster))
    }
    return(object)
  } else if (resultMax <= 0) { #largest N still less than nCluster
    object <- FindClusters(object, resolution = resolutionMax, ...)
    if (resultMax < 0) {
      warning(paste0("Get nCluster=", resultMax+nCluster, " with largest resolution=", resolutionMax, ". Can't get nCluster=", nCluster))
    }
    return(object)
  }
<<<<<<< HEAD
  result=ssanv::uniroot.integer(f=FindClustersForUniroot,interval=resolutionRange*100,
                                step.power=5,
                                object=object,nCluster=nCluster,
                                ...)
  selectedResolution=result$root/100
=======
  result <- ssanv::uniroot.integer(f = FindClustersForUniroot,
                                   interval = resolutionRange * 100,
                                   step.power = 2,
                                   object = object,
                                   nCluster = nCluster,
                                   ...)
  selectedResolution <- result$root/100
>>>>>>> d0fba68aee14797ae32a632982950b6a6252c3e8
  #browser()
  object <- FindClusters(object, resolution = selectedResolution, ...)
  return(object)
}

#temp=FindNClusters(object, nCluster=nCluster,graph.name = "ImageFeaturePCA_snn")




