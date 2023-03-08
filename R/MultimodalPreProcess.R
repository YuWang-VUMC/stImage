#' MultimodalPreProcess
#' @include runDimReduc.R
#' @param object A \code{Seurat} object prepared by \code{LoadImageFeature} or
#' \code{LoadImageFeatureVisium}
#' @param normalizeMethod the normalization method for gene expression matrix.
#' \code{log} for \code{LogNormalize}, \code{SCT} for \code{SCTransform}.
#' @param pcaDim_s number of PCs when running dimension reduction on gene
#' expression (spatial) matrix, 30 as default
#' @param pcaDim_i number of PCs when running dimension reduction on image
#' feature matrix, 30 as default
#' @param pcaDim_c number of PCs when running dimension reduction on RGB
#' quantile matrix, 6 as default
#' @param DimReducMethod dimension reduction method. Can be either \code{PCA}
#' or \code{SpatialPCA}
#' @param prefiltergenePercentCut cutoff value of percentage of spots of their
#' values larger than minimum value for prefiltering genes.
#' @param prefilterimagePercentCut cutoff value of percentage of spots of their
#' values larger than minimum value for prefiltering image features.
#' @param genePercentCut cutoff value of percentage of spots of their
#' values larger than minimum value for filtering genes after normalization.
#' @param imagePercentCut cutoff value of percentage of spots of their
#' values larger than minimum value for prefiltering image features after
#' normalization.
#' @param customGenes custom gene list provide by user. If is not \code{NULL},
#' this list will be used instead of the most variable genes.
#' @param ... Arguments passed to \code{\link{runDimReduc}}
#'
#' @return  A \code{Seurat} object
#' @importFrom SpatialPCA CreateSpatialPCAObject
#' @export
#'
#' @examples
#' \dontrun{
#' object <-
#'   MultimodalPreProcess(
#'     object,
#'     pcaDim_s = 20,
#'     pcaDim_i = 20,
#'     pcaDim_c = 6,
#'     normalizeMethod = "SCT",
#'     DimReducMethod = "PCA",
#'     prefiltergenePercentCut=0,
#'     prefilterimagePercentCut=0,
#'     genePercentCut = 0.25,
#'     imagePercentCut = 0.2
#'  )
#' }
#'
MultimodalPreProcess <- function(object,
                                 normalizeMethod = c("SCT", "log"),
                                 pcaDim_s = 30,
                                 pcaDim_i = 30,
                                 pcaDim_c = 6,
                                 DimReducMethod = c("PCA", "SpatialPCA"),
                                 prefiltergenePercentCut=0.05,
                                 prefilterimagePercentCut=0.05,
                                 genePercentCut=0.05,
                                 imagePercentCut=0.05,
                                 customGenes=NULL,
                                 ...) {
  normalizeMethod <- match.arg(normalizeMethod)
  DimReducMethod <- match.arg(DimReducMethod)

  if (!is.null(customGenes)) {
    message("## customGenes defined and wil be used for analysis.
            FindVariableFeatures will not be performed.")
  }

  #data filtering at raw data level
  message("## Data filtering at raw level")
  ### gene level
  assay <- "Spatial"
  percentCut <- prefiltergenePercentCut

  minValue <- min(object[[assay]]@counts)
  if (minValue != 0) {
    warning(paste0("Min value in ", assay, " is not equal to 0.
                   Need to confirm data filtering by >=", minValue,
                   " percent is correct!"))
  }
  geneExpressionPercent <- apply(object[[assay]]@counts, 1, function(x)
    length(which(x > minValue)) / length(x))
  geneToKept <- names(which(geneExpressionPercent >= percentCut))
  message(paste0("## ", length(geneExpressionPercent) - length(geneToKept),
                 " features in ", assay, " assay were removed. ",
                 length(geneToKept), " kept."))
  #can only do this because object is new and no information in data
  object[[assay]]@counts <- object[[assay]]@counts[geneToKept,]
  object[[assay]]@data <- object[[assay]]@data[geneToKept,]


  ###ImageFeature
  assay <- "ImageFeature"
  percentCut <- prefilterimagePercentCut

  minValue <- min(object[[assay]]@counts)
  if (minValue != 0) {
    warning(paste0("Min value in ", assay, " is not equal to 0.
                   Need to confirm data filtering by >=", minValue,
                   " percent is correct!"))
  }
  geneExpressionPercent <- apply(object[[assay]]@counts, 1, function(x)
    length(which(x>minValue)) / length(x))
  geneToKept <- names(which(geneExpressionPercent >= percentCut))
  message(paste0("## ", length(geneExpressionPercent) - length(geneToKept),
                 " features in ", assay, " assay were removed. ",
                 length(geneToKept), " kept."))
  #can only do this because object is new and no information in data
  object[[assay]]@counts <- object[[assay]]@counts[geneToKept,]
  object[[assay]]@data <- object[[assay]]@data[geneToKept,]

  #gene processing
  message("## Working on gene expression")
  DefaultAssay(object) <- "Spatial"
  if (normalizeMethod == "SCT") {
    assay <- "SCT"
    object <-
      SCTransform(object,
                  assay = "Spatial",
                  verbose = FALSE,
                  residual.features = customGenes)
    DefaultAssay(object) <- "SCT"
  } else if (normalizeMethod == "log") {
    assay <- "Spatial"
    object <- NormalizeData(object)
    if (!is.null(customGenes)) {
      VariableFeatures(object) <- customGenes
    } else {
      object <- FindVariableFeatures(object)
    }
    object <- ScaleData(object)

  } else {
    stop(paste0("normalizeMethod has to be SCT or log"))
  }
  object <- runDimReduc(
    object,
    DimReducMethod = DimReducMethod,
    assay = assay,
    percentCut = genePercentCut,
    pcaDim = pcaDim_s,
    customGenes=customGenes,
    ...
  )

  #image processing
  message("## Working on image features")
  assay <- "ImageFeature"
  DefaultAssay(object) <- assay
  VariableFeatures(object) <- rownames(object[[assay]])
  object <-
    NormalizeData(object, normalization.method = 'CLR', margin = 1)
  object[[assay]]@data[is.na(object[[assay]]@data)] <- 0
  object <- object %>% ScaleData()

  object <- runDimReduc(
    object,
    DimReducMethod = DimReducMethod,
    assay = assay,
    percentCut = imagePercentCut,
    pcaDim = pcaDim_i,
    ...
  )

  #RGB features processing
  if ("RGB" %in% Seurat::Assays(object)) {
    message("## Working on RGB features")

    assay <- "RGB"
    DefaultAssay(object) <- assay
    VariableFeatures(object) <- row.names(object[[assay]])
    rgb_norm <- 1 - (object@assays$RGB@data / 255)
    object@assays$RGB@data <- as(object = rgb_norm, Class = 'dgCMatrix')
    object <- object %>% ScaleData()

    object <- runDimReduc(
      object,
      DimReducMethod = "PCA",
      assay = assay,
      percentCut = imagePercentCut,
      pcaDim = pcaDim_c,
      ...
    )
  }

  return(object)
}
