#' DoBayesSpaceIntegration
#' @include doBayesSpace.R
#' @include findSNNClusters.R
#' @include similarityIntegration.R
#' @param dataObj A \code{Seurat} object.
#' @param clusterNum number of clusters
#' @param BayesSpaceClusterNum
#' @param assay the two assays used for integration
#'
#' @return
#' @export
#'
#' @examples
#' \dontrun{
#'   object <- DoBayesSpaceIntegration(
#'     object,
#'     clusterNum = clusterNum)
#' }
#'
DoBayesSpaceIntegration <- function(dataObj,
                                   clusterNum = 10,
                                   BayesSpaceClusterNum = clusterNum,
                                   assay = c("SCT","ImageFeature")) {
  #BayesSpace
  dataObj <- doBayesSpace(
    dataObj,
    assay = assay[1],
    nCluster = BayesSpaceClusterNum,
    name = paste0(assay[1],"BayesSpace_Cluster", clusterNum)
  )
  dataObj <- doBayesSpace(
    dataObj,
    assay = assay[2],
    nCluster = BayesSpaceClusterNum,
    name = paste0(assay[2],"BayesSpace_Cluster", clusterNum)
  )
  dataObj <- IntegrationByDistance(
    dataObj,
    distMatrix =
      paste0(assay,"BayesSpace_Cluster",clusterNum,"_Dist"),
    snnMatrix =
      paste0(assay,"BayesSpace_Cluster",clusterNum,"_snn"),
    name = c(
      "BayesSpace_DistIntegrated"
    )
  )

  dataObj <-
    findSNNClusters(
      dataObj,
      nCluster = clusterNum,
      graph.name = "BayesSpace_DistIntegrated_snn",
      nn.name = "BayesSpace_DistIntegrated_Neighbor",
      clusterColumnName = paste0("BayesSpace_DistIntegrated_Cluster",
                                 clusterNum),
      resolutionMin = 0.001,
      resolutionMax = 4,
      RunUMAP = FALSE
    )
  return(dataObj)
}

