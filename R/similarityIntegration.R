#' distToNeighborObj
#'
#' @param distMatrix distance matrix
#' @param k.param number of k
#'
#' @return a Neighbor object
#' @export
#'
distToNeighborObj <- function(distMatrix,
                              k.param = 20) {

  #distMatrix=as.matrix(bayesSpaceResultDist)
  n.cells <- nrow(distMatrix)
  knn.mat <- matrix(data = 0, ncol = k.param, nrow = n.cells)
  knd.mat <- knn.mat
  for (i in 1:n.cells) {
    knn.mat[i, ] <- order(distMatrix[i, ])[1:k.param]
    knd.mat[i, ] <- distMatrix[i, knn.mat[i, ]]
  }
  # nn.ranked <- knn.mat[, 1:k.param]
  # nn.idx=t(apply(nnMatrix,1,function(x) which(x>0)))
  # nn.dist=matrix(NA,ncol=ncol(nn.idx),nrow=nrow(nn.idx))
  # for (i in 1:nrow(nn.idx)) {
  #   nn.dist[i,]=distMatrix[i,nn.idx[i,]]
  # }

  objNeighbor <- new(
    Class = 'Neighbor',
    nn.idx = knn.mat,
    nn.dist = knd.mat,
    alg.info = list(),
    cell.names = rownames(x = distMatrix)
  )
  return(objNeighbor)
}

#' NeighborDist
#'
#' @param distMatrix distance matrix
#' @param neighborList list of neighbor
#' @param removeFirst logical value to remove the spot itself. TRUE as default
#'
#' @return a dist dataframe
#' @export
#'
NeighborDist <- function(distMatrix,
                         neighborList,
                         removeFirst = TRUE) {
  if (class(neighborList) == "Neighbor") { #make it a list
    neighborList <- lapply(data.frame(t(neighborList@nn.idx)), function(x) x)
  }
  n.cells <- length(neighborList)
  k <- length(neighborList[[1]])
  distSummary <- structure(rep(NA,n.cells), names = row.names(distMatrix))
  distToNearest <- structure(rep(NA,n.cells), names = row.names(distMatrix))
  if (removeFirst) {
    allNeighborInd = 2:k
    nearestInd = 2
  } else {
    allNeighborInd = 1:k
    nearestInd = 1
  }
  for (i in 1:n.cells) {
    #distSummary[i]=mean(distMatrix[i,objNeighbor@nn.idx[i,]])
    distSummary[i] <- mean(distMatrix[i,neighborList[[i]]][allNeighborInd])
    #need to remove self
    distToNearest[i] <- (distMatrix[i, neighborList[[i]]][nearestInd])
    #1 is self
  }
  return(data.frame(distSummary, distToNearest))
}


#' IntegrationByDistance
#'
#' @param dataObj A \code{Seurat} object.
#' @param distMatrix a list of distance Matrices if \code{dataObj} is NULL
#' @param snnMatrix a list of SNN Matrices; provided by user or from graphs slot
#' @param name naming prefix to store results related to IntegrationByDistance
#' @param graphs.name  Name of (s)nn output in graphs slot
#' @param nn.name  Name of knn output on which to run UMAP
#' @param integrated.name name of column in meta data to store the clustering
#' results
#' @param integrated.weight.name name of column in meta data to store the
#' weights of distance
#' @param k.param  number of k
#' @param minDiff minimum resdiue to add when distance calculation
#'
#' @return a seurat object
#' @export
#'
IntegrationByDistance <- function(
    dataObj = NULL,
    distMatrix = c("SCTBayesSpace_dist_Cluster10",
                 "ImageFeatureBayesSpace_dist_Cluster10"),
    snnMatrix = c("SCTBayesSpace_snn_Cluster10",
                "ImageFeatureBayesSpace_snn_Cluster10"),
    name = "DistIntegrated",
    graphs.name = paste0(name, c("_nn", "_snn")),
    nn.name = paste0(name, "_Neighbor"),
    integrated.name = paste0(name, "_Dist"),
    integrated.weight.name = paste0(name, "_DistWeight"),
    k.param = 20,
    minDiff = 10^-4) {

  if (is.null(dataObj)) {
    #not defined dataObj, distMatrix should be a list of distance Matrices
    distMatrix1 <- distMatrix[[1]]
    distMatrix2 <- distMatrix[[2]]
    snnMatrix1 <- snnMatrix[[1]]
    snnMatrix2 <- snnMatrix[[2]]
  } else { #defined dataObj, get distMatrix from graphs
    distMatrix1 <- as.matrix(dataObj[[distMatrix[1]]])
    distMatrix2 <- as.matrix(dataObj[[distMatrix[2]]])

    snnMatrix1 <- dataObj@graphs[[snnMatrix[1]]]
    snnMatrix2 <- dataObj@graphs[[snnMatrix[2]]]
  }
  snnMatrix1[snnMatrix1 == 0] <- NA
  snnMatrix2[snnMatrix2 == 0] <- NA

  distMatrix1NeighborObj <- distToNeighborObj(distMatrix1, k.param = k.param)
  distMatrix2NeighborObj <- distToNeighborObj(distMatrix2, k.param = k.param)

  distMatrix1InNeighbor1 <- NeighborDist(distMatrix1, distMatrix1NeighborObj)
  distMatrix1InNeighbor2 <- NeighborDist(distMatrix1, distMatrix2NeighborObj)
  distMatrix2InNeighbor1 <- NeighborDist(distMatrix2, distMatrix1NeighborObj)
  distMatrix2InNeighbor2 <- NeighborDist(distMatrix2, distMatrix2NeighborObj)

  snnMatrix1SmallIndexList <- list()
  for (i in 1:nrow(snnMatrix1)) {
    snnMatrix1SmallIndexList[[i]] <-
      order(snnMatrix1[i, ])[1:k.param]
    #subject Ind with smallest non zero SNN (Jarcard index)
  }
  snnMatrix2SmallIndexList <- list()
  for (i in 1:nrow(snnMatrix2)) {
    snnMatrix2SmallIndexList[[i]] <-
      order(snnMatrix2[i, ])[1:k.param]
    #subject Ind with smallest non zero SNN (Jarcard index)
  }

  distMatrix1InSmallIndex1 <- NeighborDist(distMatrix1,
                                           snnMatrix1SmallIndexList,
                                           removeFirst = FALSE)
  #distance to furthest neighbors. don't need to remove first (self)
  distMatrix2InSmallIndex2 <- NeighborDist(distMatrix2,
                                           snnMatrix2SmallIndexList,
                                           removeFirst = FALSE)

  #caculate values
  theta1InNeighbor1 <-
    exp(-pmax((distMatrix1InNeighbor1[,1] - distMatrix1InNeighbor1[,2]), 0) /
          (distMatrix1InSmallIndex1[,1] - distMatrix1InNeighbor1[,2] + minDiff))
  theta1InNeighbor2 <-
    exp(-pmax((distMatrix1InNeighbor2[,1] - distMatrix1InNeighbor1[,2]), 0) /
          (distMatrix1InSmallIndex1[,1] - distMatrix1InNeighbor1[,2] + minDiff))

  theta2InNeighbor1 <-
    exp(-pmax((distMatrix2InNeighbor1[,1] - distMatrix2InNeighbor2[,2]), 0) /
          (distMatrix2InSmallIndex2[,1] - distMatrix2InNeighbor2[,2] + minDiff))
  theta2InNeighbor2 <-
    exp(-pmax((distMatrix2InNeighbor2[,1] - distMatrix2InNeighbor2[,2]), 0) /
          (distMatrix2InSmallIndex2[,1] - distMatrix2InNeighbor2[,2] + minDiff))

  sDist1 <- theta1InNeighbor1 / (theta1InNeighbor2 + minDiff)
  sDist2 <- theta2InNeighbor2 / (theta2InNeighbor1 + minDiff)

  #weightDist1 <- exp(sDist1)/(exp(sDist1) + exp(sDist2))
  #weightDist2 <- exp(sDist2)/(exp(sDist1) + exp(sDist2))
  ##Changes to deal with very large sDist. Will get the same result
  weightDist1 <- 1 / (1 + exp(sDist2 - sDist1))
  weightDist2 <- 1 / (1 + exp(sDist1 - sDist2))

  ##Don't need this after changing caculating weightDist method
  # allNaInd=which(is.na(weightDist1) & is.na(weightDist2))
  # if (any(allNaInd)) {
  #   weightDist1[allNaInd] <- 0.5
  #   weightDist2[allNaInd] <- 0.5
  # }
  # weightDist1[is.na(weightDist1)] <- 1
  # weightDist2[is.na(weightDist2)] <- 1

  distMatrix1Weighted <- distMatrix1 * weightDist1
  distMatrix2Weighted <- distMatrix2 * weightDist2
  distMatrixWeighted <- distMatrix1Weighted + distMatrix2Weighted
  distMatrixWeighted_Weight <- data.frame(weightDist1, weightDist2)
  row.names(distMatrixWeighted_Weight) <- row.names(distMatrixWeighted)
  if (is.null(dataObj)) {
    #not defined dataObj, distMatrix should be a list of distance Matrixs
    return(distMatrixWeighted)
  } else {
    dataObj@graphs[[integrated.name]] <- distMatrixWeighted
    dataObj@graphs[[paste0(distMatrix[1],"_Weighted")]] <- distMatrix1Weighted
    dataObj@graphs[[paste0(distMatrix[2],"_Weighted")]] <- distMatrix2Weighted
    dataObj@graphs[[integrated.weight.name]] <- distMatrixWeighted_Weight
    dataObj <- FindNeighborsByDistance(dataObj,distMatrixWeighted,
                                       graphs.name = graphs.name,
                                       nn.name = nn.name,
                                       k.param = k.param)

    return(dataObj)
  }
}


#' FindNeighborsByDistance
#' @param dataObj  A \code{Seurat} object.
#' @param distMatrix  distance matrix in graphs slot
#' @param prune.SNN Sets the cutoff for acceptable Jaccard index when
#' computing the neighborhood overlap for the SNN construction. Any edges with
#' values less than or equal to this will be set to 0 and removed from the SNN
#' graph. Essentially sets the stringency of pruning (0 --- no pruning, 1 ---
#' prune everything).
#' @param assay assay name
#' @param graphs.name graph names: nn and snn
#' @param nn.name Name of knn output on which to run UMAP
#' @param reduction.name Name of projected UMAP to store in the query
#' @param reduction.key Value for the projected UMAP key
#' @param k.param number of k
#' @importFrom Seurat FindNeighbors DefaultAssay DefaultAssay<- RunUMAP
#' @return a seurat object
#' @export
#'
FindNeighborsByDistance <- function(dataObj = NULL,
                                    distMatrix,
                                    prune.SNN = 1/15,
                                    assay = NULL,
                                    graphs.name = c("nn", "snn"),
                                    nn.name = "neighbor",
                                    reduction.name = paste0(nn.name, 'UMAP'),
                                    reduction.key = paste0(nn.name, 'UMAP_'),
                                    k.param = 20) {
  if (!is.null(dataObj) & is.character(distMatrix)) {
    #defined dataObj, extract distMatrix from graphs
    distMatrix <- as.matrix(dataObj@graphs[[distMatrix]])
  }
  distMatrixNeighborsList <- FindNeighbors(distMatrix,
                                           distance.matrix = TRUE,
                                           prune.SNN = prune.SNN,
                                           force.recalc = TRUE,
                                           k.param = k.param)

  if (!is.null(assay)) {
    DefaultAssay(distMatrixNeighborsList[[1]]) <- assay
    DefaultAssay(distMatrixNeighborsList[[2]]) <- assay
  }

  if (is.null(dataObj)) {
    return(distMatrixNeighborsList)
  } else { #defined dataObj, add nn, snn, neighbor to object
    dataObj[[graphs.name[1]]] <- distMatrixNeighborsList[["nn"]]
    dataObj[[graphs.name[2]]] <- distMatrixNeighborsList[["snn"]]

    #get NeighborObj for umap
    objNeighbor <- distToNeighborObj(distMatrix, k.param = k.param)
    dataObj@neighbors[[nn.name]] <- objNeighbor
    dataObj <- RunUMAP(dataObj,
                       nn.name = nn.name,
                       reduction.name = reduction.name,
                       reduction.key = reduction.key)
    return(dataObj)
  }
}


