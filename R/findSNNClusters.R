#' Title
#'
#' @param object
#' @param reduction
#' @param graph.name
#' @param nn.name
#' @param dims
#' @param resolution
#' @param assay
#' @param nCluster
#' @param clusterColumnName
#' @param reduction.name
#' @param reduction.key
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
findSNNClusters <- function(object,
                            reduction = "pca",
                            graph.name = paste0(reduction, c("_nn", "_snn")),
                            nn.name = NULL,
                            dims = 1:30,
                            resolution=0.8,
                            assay = DefaultAssay(object),
                            nCluster = NULL,
                            clusterColumnName=if (!is.null(nCluster)) paste0(reduction,"_cluster",nCluster) else NULL,
                            reduction.name =ifelse(!is.null(reduction),paste0(reduction,'UMAP'),paste0(nn.name,'UMAP')),
                            reduction.key = ifelse(!is.null(reduction),paste0(reduction,'UMAP_'),paste0(nn.name,'UMAP_')),
                            RunUMAP=TRUE,
                            ...) {
  #dims <- intersect(dims, 1:ncol(object@reductions[[reduction]]@cell.embeddings))
  if (length(graph.name) == 1 && graph.name %in% names(object@graphs)) { #SNN graph.name exist and no need to do FindNeighbors
    message(paste0("Using ",graph.name," to FindClusters. Skipping FindNeighbors."))
    snnGraphName <- graph.name
  } else {
    dims=intersect(dims,1:ncol(object@reductions[[reduction]]@cell.embeddings))
    object <- FindNeighbors(object, reduction = reduction, graph.name = graph.name, dims = dims)
    snnGraphName <- graph.name[2]
  }
  if (!is.null(nCluster)) { #have target nCluster number
    object <- FindNClusters(object, nCluster = nCluster, resolution = resolution, graph.name = snnGraphName, ...)
  } else {
    object <- FindClusters(object, verbose = FALSE, graph.name = snnGraphName, resolution = resolution)
  }
  if (!is.null(clusterColumnName)) { #need to store cluster result in defined clusterColumnName
    object@meta.data[[clusterColumnName]] <- object@meta.data[["seurat_clusters"]]
  }

  #RunUMAP
  if (RunUMAP) {
    if (!is.null(nn.name)) { #nn.name defined. Use nn.name to run RunUMAP
      object <- RunUMAP(object,
                        nn.name = nn.name,
                        reduction.name = reduction.name,
                        reduction.key = reduction.key)
    } else { #nn.name NOT defined. Use reduction to run RunUMAP
      dims=intersect(dims,c(1:ncol(object[[reduction]])))
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
