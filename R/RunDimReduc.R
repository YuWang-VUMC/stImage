#' RunDimReduc
#' @include SpatialPCAWorkflow.R
#' @param object A \code{Seurat} object
#' @param DimReducMethod dimension reduction method. Can be either \code{PCA}
#' or \code{SpatialPCA}
#' @param assay assay used for dimension reduction analysis. i.e., \code{SCT}
#' for gene matrix when preprocessed by \code{SCTransform}
#' @param percentCut cutoff value of percentage of spots of their
#' values larger than minimum value
#' @param pcaDim number of PCs when running dimension reduction
#' @param customGenes vector of custom genes to run dimension reduction step
#' @param ... potentially needed Arguments
#'
#' @return a seurat object
#' @importFrom SpatialPCA CreateSpatialPCAObject
#' @importFrom Seurat VariableFeatures RunPCA CreateDimReducObject
#' @importFrom Seurat VariableFeatures<-
#' @export
#'
RunDimReduc <- function(object,
                        DimReducMethod = "PCA",
                        assay = "SCT",
                        percentCut = 0.05,
                        pcaDim = 30,
                        customGenes = NULL,
                        ...) {

  DimReducName <- paste0(assay, "", DimReducMethod)
  DimReducKeyName <- paste0(DimReducName, "_")

  if(DimReducMethod == "PCA") {
    if (is.null(customGenes)) {
      #remove less frequent genes
      geneExpressionPercent <- apply(object[[assay]]@counts, 1,
                                     function(x) length(which(x>0)) / length(x))
      VariableFeatures(object) <-
        setdiff(VariableFeatures(object),
                names(which(geneExpressionPercent <= percentCut)))
    } else {
      VariableFeatures(object) <-
        intersect(customGenes, row.names(object[[assay]]@counts))
    }

    object <- RunPCA(object,
                     assay = assay,
                     npcs = pcaDim,
                     verbose = FALSE,
                     reduction.name = DimReducName,
                     reduction.key = DimReducKeyName)
  } else if(DimReducMethod == "SpatialPCA") {
    rawcount <- as.matrix(object[[assay]]@counts)
    #need to do gene filter here as CreateSeuratObject can't
    #deal with min.cells and min.features correctly to keep
    #all cells with >0 gene
    min.loctions <- as.integer(ncol(rawcount) * percentCut)
    featureExpressionN <- apply(object[[assay]]@counts, 1,
                                function(x) length(which(x>0)))
    rawcountFiltered <- rawcount[which(featureExpressionN > min.loctions), ]
    rawcountFiltered <- rawcountFiltered[, which(colSums(rawcountFiltered) > 0)]
    location <- as.matrix(object@images[[1]]@coordinates)
    if (("imagecol" %in% colnames(location)) &
        "imagerow" %in% colnames(location) ) { #10X Visum data
      locationFiltered <- location[colnames(rawcountFiltered),
                                   c("imagecol","imagerow")]
    } else { #ST data
      locationFiltered <- location[colnames(rawcountFiltered), ]
    }

    ##Can't use CreateSpatialPCAObject here as it will use
    ##SCTransform for normalization ,which is not good for image features
    ##Use @scale.data directly
    SPCAobj <- new(Class = "SpatialPCA",
                   counts = rawcountFiltered,
                   location = locationFiltered,
                   project = object@project.name)
    if (!is.null(customGenes)) {
      VariableFeatures(object) <-
        intersect(row.names(rawcountFiltered), customGenes)
      message(paste0(length(customGenes)," customGenes defined and ",
                     length(VariableFeatures(object)),
                     " were overlapped with filtered data and used as features
                     for SpatialPCA"))
    }
    SPCAobj@normalized_expr <-
      object[[assay]]@scale.data[intersect(row.names(rawcountFiltered),
                                           VariableFeatures(object)), ]
    SPCAobj <- SpatialPCAWorkflow(SPCAobj, SpatialPCnum = pcaDim)
    SPCApcs <- t(SPCAobj@SpatialPCs)
    colnames(SPCApcs) <- paste0(DimReducKeyName, 1:ncol(SPCApcs))
    object <- subset(object, cells = row.names(SPCApcs))
    object[[DimReducName]] <- CreateDimReducObject(embeddings = SPCApcs,
                                                   key = DimReducKeyName,
                                                   assay = assay)
  } else {
    stop(paste0("DimReducMethod has to be PCA or SpatialPCA"))
  }
  return(object)
}
