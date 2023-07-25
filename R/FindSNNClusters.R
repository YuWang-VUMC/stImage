#' FindSNNClusters
#' @param object A \code{Seurat} object.
#' @param reduction dimension reduction method
#' @param graph.name naming parameter for stored (S)NN graph
#' @param nn.name Name of knn output on which to run UMAP
#' @param dims which dimensions to use as input features
#' @param resolution value of the resolution parameter, use a value above
#' (below) 1.0 if you want to obtain a larger (smaller) number of communities.
#' @param assay Assay to use in construction of (S)NN; used only when dims is
#' NULL
#' @param nCluster number of clusters
#' @param k.param Defines k for the k-nearest neighbor algorithm
#' @param clusterColumnName column name in meta data to store the clustering
#' results
#' @param reduction.name name of projected UMAP to store in the query
#' @param reduction.key value for the projected UMAP key
#' @param RunUMAP logical value to define if or not to run RunUMAP
#' @param ... Arguments passed to \code{\link{FindNClusters}}
#' @importFrom Seurat FindNeighbors DefaultAssay RunUMAP FindClusters
#' @importFrom Seurat DefaultAssay<-
#' @return a seurat object
#' @export
#'
#' @examples \dontrun{
#' object <- FindSNNClusters(
#'   object,
#'   dims = 1:20,
#'   reduction = "SCTPCA",
#'   nCluster = clusterNum
#'   )
#' }
#'
FindSNNClusters <- function(
    object,
    reduction = "pca",
    graph.name = paste0(reduction, c("_nn", "_snn")),
    nn.name = NULL,
    dims = 1:30,
    resolution = 0.8,
    assay = DefaultAssay(object),
    nCluster = NULL,
    k.param = 20,
    clusterColumnName = if(!is.null(nCluster))
      paste0(reduction, "_cluster", nCluster) else NULL,
    reduction.name = ifelse(!is.null(reduction),
                            paste0(reduction,'UMAP'), paste0(nn.name,'UMAP')),
    reduction.key = ifelse(!is.null(reduction),
                           paste0(reduction,'UMAP_'), paste0(nn.name,'UMAP_')),
    RunUMAP = TRUE,
    ...) {
  DefaultAssay(object) <- assay
  #dims <- intersect(dims,
  #                  1:ncol(object@reductions[[reduction]]@cell.embeddings))
  if (length(graph.name) == 1 && graph.name %in% names(object@graphs)) {
    #SNN graph.name exist and no need to do FindNeighbors
    message(paste0("Using ", graph.name,
                   " to FindClusters. Skipping FindNeighbors."))
    snnGraphName <- graph.name
  } else {
    dims <- intersect(dims,
                      1:ncol(object@reductions[[reduction]]@cell.embeddings))
    object <- FindNeighbors(object,
                            assay = assay,
                            reduction = reduction,
                            graph.name = graph.name,
                            dims = dims,
                            k.param=k.param)
    snnGraphName <- graph.name[2]
  }
  if (!is.null(nCluster)) { #have target nCluster number
    object <- FindNClusters(object,
                            nCluster = nCluster,
                            resolution = resolution,
                            graph.name = snnGraphName, ...)
  } else {
    object <- FindClusters(object,
                           verbose = FALSE,
                           graph.name = snnGraphName,
                           resolution = resolution)
  }
  if (!is.null(clusterColumnName)) {
    #need to store cluster result in defined clusterColumnName
    object@meta.data[[clusterColumnName]] <-
      object@meta.data[["seurat_clusters"]]
  }

  #RunUMAP
  if (RunUMAP) {
    if (!is.null(nn.name)) { #nn.name defined. Use nn.name to run RunUMAP
      object <- RunUMAP(object,
                        nn.name = nn.name,
                        reduction.name = reduction.name,
                        reduction.key = reduction.key)
    } else { #nn.name NOT defined. Use reduction to run RunUMAP
      dims <- intersect(dims, c(1:ncol(object[[reduction]])))
      object <- RunUMAP(object,
                        reduction = reduction,
                        dims = dims,
                        assay = assay,
                        reduction.name = reduction.name,
                        reduction.key = reduction.key)
    }
  }

  return(object)
}
