#' DoSpectrum
#' @param object A \code{Seurat} object.
#' @param reduction.list list of dimension reduction matrices for integration
#' @param graphs graph slot used for similarity calculation if don't define
#' \code{reduction.list}.
#' @param clusterEachModality cluster each modality if set to TRUE. TRUE as
#' default.
#' @param nCluster the number of clusters
#' @param clusterColumnName column name to store clustering results
#' @param integrated.name integrated name used in graph slot
#' @param nn.name Multimodal neighbor object name in neighbors slot
#' @param reduction.name Name of projected UMAP to store in the query
#' @param reduction.key Value for the projected UMAP key
#' @importFrom Spectrum CNN_kernel
#' @importFrom Spectrum integrate_similarity_matrices
#' @importFrom Spectrum cluster_similarity
#' @importFrom Seurat RunUMAP
#' @return a seurat object
#' @export
#'
DoSpectrum <- function(object,
                       reduction.list = NULL,
                       graphs = NULL,
                       clusterEachModality = TRUE,
                       nCluster = 4,
                       clusterColumnName = paste0("Spectrum_Cluster", nCluster),
                       integrated.name = "SpectrumIntegrated_Similarity",
                       nn.name = "SpectrumNeighbor",
                       reduction.name = paste0(nn.name,'UMAP'),
                       reduction.key = paste0(nn.name,'UMAP_')) {
  if (is.null(reduction.list) & is.null(graphs)) {
    stop("Need define at least one reduction.list or one graphs")
  }
  #similarity
  if (!is.null(reduction.list)) { #defined reduction.list, not defined graphs,
    #get similarity from reduction.list by Spectrum
    similarityList=NULL
    for (i in 1:length(reduction.list)) {
      similarityList[[reduction.list[[i]]]] <-
        CNN_kernel(t(object@reductions[[reduction.list[[i]]]]@cell.embeddings))
    }
  } else if (!is.null(graphs)) {  #defined graphs. get similarity from graphs
    similarityList <- object@graphs[graphs]
  }

  if (clusterEachModality) { #doing single cluster
    for (i in 1:length(similarityList)) {
      clusterResult <-
        Spectrum::cluster_similarity(similarityList[[i]],
                                     k = nCluster,
                                     clusteralg = 'GMM')
      if (!is.null(reduction.list)) {
        clusterColumnNameSingle <-
          paste0(reduction.list[[i]], "Spectrum_Cluster", nCluster)
      } else {
        clusterColumnNameSingle <-
          paste0(graphs[i], "Spectrum_Cluster", nCluster)
      }
      object@meta.data[[clusterColumnNameSingle]] <- clusterResult
    }
  }

  #integration similarity matrix
  sIntegrated <- Spectrum::integrate_similarity_matrices(similarityList)
  SpectrumClusterResult <-
    Spectrum::cluster_similarity(sIntegrated,
                                 k = nCluster,
                                 clusteralg='GMM')
  object@meta.data[[clusterColumnName]] <- SpectrumClusterResult

  object@graphs[[integrated.name]] <- sIntegrated
  #get NeighborObj for umap
  objNeighbor <- DistToNeighborObj(1-sIntegrated,k.param=20)
  object@neighbors[[nn.name]] <- objNeighbor
  object <- RunUMAP(object,
                     nn.name = nn.name,
                     reduction.name = reduction.name,
                     reduction.key = reduction.key)

  return(object)
}
