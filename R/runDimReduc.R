#' Title
#'
#' @param object
#' @param DimReducMethod
#' @param assay
#' @param percentCut
#' @param pcaDim
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
runDimReduc <- function(object,
                        DimReducMethod = "PCA",
                        assay = "SCT",
                        percentCut = 0.05,
                        pcaDim=30,...) {

  DimReducName <- paste0(assay, "", DimReducMethod)
  DimReducKeyName <- paste0(DimReducName, "_")

  if(DimReducMethod == "PCA") {
    #remove less frequent genes
    geneExpressionPercent <- apply(object[[assay]]@counts, 1, function(x) length(which(x>0)) / length(x))
    VariableFeatures(object) <- setdiff(VariableFeatures(object), names(which(geneExpressionPercent <= percentCut)))

    object <- RunPCA(object, assay = assay, npcs = pcaDim, verbose = FALSE,
                     reduction.name = DimReducName,reduction.key = DimReducKeyName)
  } else if(DimReducMethod == "SpatialPCA") {
    rawcount <- as.matrix(object[[assay]]@counts)
    #need to do gene filter here as CreateSeuratObject can't deal with min.cells and min.features correctly to keep all cells with >0 gene
    min.loctions <- as.integer(ncol(rawcount)*percentCut)
    featureExpressionN <- apply(object[[assay]]@counts, 1, function(x) length(which(x>0)))
    rawcountFiltered <- rawcount[which(featureExpressionN > min.loctions), ]
    rawcountFiltered <- rawcountFiltered[ , which(colSums(rawcountFiltered) > 0)]
    location <- as.matrix(object@images$images@coordinates)
    locationFiltered <- location[colnames(rawcountFiltered), ]

    ##Can't use CreateSpatialPCAObject here as it will use SCTransform for normalization ,which is not good for image features
    ##Use @scale.data directly
    SPCAobj <- new(Class = "SpatialPCA",
                   counts = rawcountFiltered,
                   location = locationFiltered,
                  project = object@project.name)
    SPCAobj@normalized_expr <- object[[assay]]@scale.data[intersect(row.names(rawcountFiltered), VariableFeatures(object)), ]
    SPCAobj <- SpatialPCAWorkflow(SPCAobj, SpatialPCnum = pcaDim)
    SPCApcs <- t(SPCAobj@SpatialPCs)
    colnames(SPCApcs) <- paste0(DimReducKeyName, 1:ncol(SPCApcs))
    object <- subset(object, cells=row.names(SPCApcs))
    object[[DimReducName]] <- CreateDimReducObject(embeddings = SPCApcs,
                                                   key = DimReducKeyName,
                                                   assay = assay)
  } else {
    stop(paste0("DimReducMethod has to be PCA or SpatialPCA"))
  }
  return(object)
}
