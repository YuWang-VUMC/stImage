#' SpatialAwareProcess
#'
#' @param object  A \code{Seurat} object.
#' @param platform platform of ST data
#' @param assay Assay to use in construction of (S)NN
#' @param slot slot of assays used for processing
#' @param SpatialAwareMethod spatial-aware processing method, one of
#' "BayesSpace", "SpatialPCA", and "stLearn"
#' @param clusterNum number of clusters
#' @param pcaDim_g number of PCs when running dimension reduction on gene
#' expression (spatial) matrix, 30 as default
#' @param pcaDim_i number of PCs when running dimension reduction on image
#' feature matrix, 30 as default
#' @param percentCut cutoff value of percentage of spots of their
#' values larger than minimum value for filtering genes after normalization.
#' @param customGenes custom gene sets when performing SpatialPCA, default NULL.
#' @param distanceMethod method for distance calculation, \code{cosine} or
#'\code{euclidean}
#' @param stLearnweights define of normalization weights, based on 1)
#' \code{weights_matrix_all} consider gene level correlation in weight;
#' 2) weights from morphological Similarly and physical distance;
#' 3) physical distance only.
#'
#' @return A \code{Seurat} object
#' @export
#'
#' @examples \dontrun{
#' object <-
#'   SpatialAwareProcess(
#'     object,
#'     patform = "ST",
#'     assay = "SCT",
#'     slot = "data",
#'     SpatialAwareMethod = "stLearn",
#'     clusterNum = 8,
#'     pcaDim_g = 20,
#'     pcaDim_i = 20,
#'     percentCut = 0.2,
#'     customGenes = NULL,
#'     distanceMethod = "cosine",
#'     stLearnweights = "weights_matrix_all"
#'  )
#' }
#'
SpatialAwareProcess <- function(object,
                                platform = c("ST", "Visium"),
                                assay = c("Spatial", "SCT", "ImageFeature"),
                                slot = c("counts", "data", "scale.data"),
                                SpatialAwareMethod = c("SpatialPCA",
                                                       "BayesSpace",
                                                       "stLearn"),
                                clusterNum = NULL,
                                pcaDim_g = 30,
                                pcaDim_i = 30,
                                percentCut = 0.2,
                                customGenes = NULL,
                                distanceMethod = c("cosine", "euclidean"),
                                stLearnweights = c("weights_matrix_all",
                                                   "weights_matrix_pd_md",
                                                   "physical_distance")) {
  SpatialAwareMethod <- match.arg(SpatialAwareMethod)
  if(grepl(assay, "ImageFeature")) {
    pcaDim <- pcaDim_i
  } else {
    pcaDim <- pcaDim_g
  }

  if(SpatialAwareMethod == "BayesSpace") {
    object <- DoBayesSpace(object,
                           platform=platform,
                           assay=assay,
                           nCluster=clusterNum)
  }
  if(SpatialAwareMethod == "SpatialPCA") {
    object <- RunDimReduc(
      object,
      DimReducMethod = "SpatialPCA",
      assay = assay,
      percentCut = percentCut,
      pcaDim = pcaDim,
      customGenes=customGenes
    )
    object <- FindSNNClusters(object,
                              reduction = paste0(assay, "SpatialPCA"),
                              nCluster = clusterNum)
  }
  if(SpatialAwareMethod == "stLearn") {
    object <- DostLearn(object,
                        platform = platform,
                        pcaDim_i = pcaDim_i,
                        pcaDim_g = pcaDim_g,
                        assay = assay,
                        dataSlot = slot,
                        geneDimReducName = paste0(assay, "PCA"),#
                        distanceMethod = distanceMethod,
                        weights=stLearnweights)
    object <- RunDimReduc(object,
                          assay = paste0(assay, "NormalizedByImage"),
                          percentCut = percentCut,
                          pcaDim = pcaDim_g)
    object <- FindSNNClusters(object,
                              dims = 1:pcaDim_g,
                              reduction = paste0(assay, "NormalizedByImagePCA"),
                              nCluster = clusterNum)
  }

  return(object)
}
