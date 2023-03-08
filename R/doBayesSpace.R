#' doBayesSpace
#' @inheritParams SingleCellExperiment
#' @inheritParams BayesSpace
#' @param dataObj A \code{Seurat} object.
#' @param platform ST or Visium platform
#' @param assay assay of the object used
#' @param name column name in meta data to store the clustering results
#' @param graphs.name naming parameter for stored (S)NN graph
#' @param graphNameSimilarity naming parameter for stored similarity
#' @param graphNameDist naming parameter for stored distance
#' @param graphNameChain naming parameter for stored chain matrix
#' @param nn.name Name of knn output on which to run UMAP
#' @param reduction dimension reduction matrix in \code{reductions} slot
#' @param nCluster  the number of clusters
#' @importFrom SingleCellExperiment SingleCellExperiment
#' @importFrom BayesSpace spatialPreprocess spatialCluster
#' @importFrom rhdf5 h5read
#' @return
#' @export
#'
#' @examples
#' \dontrun{
#'   object <- doBayesSpace(
#'     object,
#'     platform="ST",
#'     assay="SCT",
#'     nCluster=clusterNum)
#' }
#'
doBayesSpace <- function(dataObj,
                         platform = c("ST", "Visium"),
                         assay = c("SCT", "Spatial", "ImageFeature"),
                         name = paste0(assay, "BayesSpace_Cluster", nCluster),
                         graphs.name = paste0(name, c("_nn", "_snn")),
                         graphNameSimilarity = paste0(name, c("_Similarity")),
                         graphNameDist = paste0(name, c("_Dist")),
                         graphNameChain = paste0(name, c("_ChainMatrix")),
                         nn.name = paste0(name, c("_Neighbor")),
                         reduction = NULL,
                         nCluster = 4) {

  colData <- dataObj@images[[1]]@coordinates
  if (("imagecol" %in% colnames(colData)) &
      "imagerow" %in% colnames(colData) ) { #10X Visum data
    colData <- colData[, c("row","col")]#use row, not imagerow for BayesSpace
    platform="Visium"
  } else { #ST data
    platform="ST"
  }
  colnames(colData)=gsub("y","row",colnames(colData))
  colnames(colData)=gsub("x","col",colnames(colData))


  dataTable=as.matrix(dataObj[[assay]]@data)
  sce <- SingleCellExperiment(assays=list(logcounts=dataTable),
                                                    colData=colData)
  if (!is.null(reduction)) { #defined reduction, use reduction in dataObj
    reducedDims(sce) <-
      list(PCA = dataObj@reductions[[reduction]]@cell.embeddings)
    sce <- spatialPreprocess(sce,
                             log.normalize = FALSE,
                             platform = platform,
                             skip.PCA = TRUE)
  } else {
    sce <- spatialPreprocess(sce,
                             log.normalize = FALSE,
                             platform = platform)
  }

  sce <- spatialCluster(sce,
                        q = nCluster,
                        nrep = 10000,
                        burn.in = 500,
                        save.chain = TRUE)
  dataObj@meta.data[[name]] <- colData(sce)$spatial.cluster

  bayesSpaceResultChainMatrix <- rhdf5::h5read(metadata(sce)$chain.h5, "z")
  colnames(bayesSpaceResultChainMatrix) <- colnames(sce)
  dataObj@graphs[[graphNameChain]] <- bayesSpaceResultChainMatrix

  #bayesSpaceResultChainMatrix

  bayesSpaceResultSimilarity <-
    bayesSpaceChainMatrixSimilarity(bayesSpaceResultChainMatrix)
  colnames(bayesSpaceResultSimilarity) <- colnames(dataTable)
  row.names(bayesSpaceResultSimilarity) <- colnames(dataTable)
  bayesSpaceResultDist <- as.dist(1 - bayesSpaceResultSimilarity)

  dataObj@graphs[[graphNameSimilarity]] <- bayesSpaceResultSimilarity
  dataObj@graphs[[graphNameDist]] <- bayesSpaceResultDist

  prune.SNN=1/100
  dataObj <- FindNeighborsByDistance(
    dataObj,
    distMatrix = as.matrix(bayesSpaceResultDist),
    assay = assay,
    prune.SNN = prune.SNN,
    graphs.name = graphs.name,
    nn.name = nn.name
    )
  return(dataObj)

}


#' bayesSpaceChainMatrixSimilarity
#'
#' @param bayesSpaceResultChainMatrix  MCMC chain saved in object
#'
#' @return
#' @export
#'
#' @examples
bayesSpaceChainMatrixSimilarity <- function(bayesSpaceResultChainMatrix) {
  nRep <- nrow(bayesSpaceResultChainMatrix)
  nSample <- ncol(bayesSpaceResultChainMatrix)

  bayesSpaceResChainMatrixSimValue <-
    combn(ncol(bayesSpaceResultChainMatrix), 2,
          function(x) sum(bayesSpaceResultChainMatrix[,x[1]] ==
                            bayesSpaceResultChainMatrix[,x[2]])) / nRep
  bayesSpaceResChainMatrixSimValue[bayesSpaceResChainMatrixSimValue==1] <-
    1 - 1/nRep

  bayesSpaceResSimMatrix <- matrix(NA, ncol = nSample, nrow = nSample)

  #lower.tri and upper.tri are all column order.
  #So lower.tri is correct order upper.tri is incorrect order
  #use lower.tri first and t(), and use lower.tri again
  bayesSpaceResSimMatrix[lower.tri(bayesSpaceResSimMatrix)] <-
    bayesSpaceResChainMatrixSimValue
  bayesSpaceResSimMatrix <- t(bayesSpaceResSimMatrix)

  bayesSpaceResSimMatrix[lower.tri(bayesSpaceResSimMatrix)] <-
    bayesSpaceResChainMatrixSimValue
  diag(bayesSpaceResSimMatrix) <- 1

  return(bayesSpaceResSimMatrix)

}
