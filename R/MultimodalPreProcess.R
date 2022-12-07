#' Title
#'
#' @param object
#' @param normalizeMethod
#' @param pcaDim_s
#' @param pcaDim_i
#' @param pcaDim_c
#' @param DimReducMethod
#' @param genePercentCut
#' @param imagePercentCut
#' @param ...
#'
#' @return
#' @importFrom SpatialPCA CreateSpatialPCAObject
#' @export
#'
#' @examples
MultimodalPreProcess <- function(object,
                                 normalizeMethod = c("SCT", "log"),
                                 pcaDim_s = 30,
                                 pcaDim_i = 30,
                                 pcaDim_c = 6,
                                 DimReducMethod = c("PCA", "SpatialPCA"),
                                 genePercentCut=0.05,
                                 customGenes=NULL,
                                 imagePercentCut=0.05,
                                 ...) {
  normalizeMethod <- match.arg(normalizeMethod)
  DimReducMethod <- match.arg(DimReducMethod)

  if (!is.null(customGenes)) {
    message("## customGenes defined and wil be used for analysis. FindVariableFeatures will not be performed.")
  }

  #data filtering at raw data level
  message("## Data filtering at raw level")
  ### gene level
  assay="Spatial"
  percentCut=genePercentCut

  minValue=min(object[[assay]]@counts)
  if (minValue!=0) {
    warning(paste0("Min value in ",assay," is not equal to 0. Need to confirm data filtering by >=",minValue,
    " percent is correct!"))
  }
  geneExpressionPercent <- apply(object[[assay]]@counts, 1, function(x) length(which(x>minValue)) / length(x))
  geneToKept=names(which(geneExpressionPercent>=percentCut))
  message(paste0("## ",length(geneExpressionPercent)-length(geneToKept)," features in ",assay," assay were removed. ",length(geneToKept)," kept."))
  #can only do this because object is new and no information in data
  object[[assay]]@counts=object[[assay]]@counts[geneToKept,]
  object[[assay]]@data=object[[assay]]@data[geneToKept,]


  ###ImageFeature
  assay="ImageFeature"
  percentCut=imagePercentCut

  minValue=min(object[[assay]]@counts)
  if (minValue!=0) {
    warning(paste0("Min value in ",assay," is not equal to 0. Need to confirm data filtering by >=",minValue,
                   " percent is correct!"))
  }
  geneExpressionPercent <- apply(object[[assay]]@counts, 1, function(x) length(which(x>minValue)) / length(x))
  geneToKept=names(which(geneExpressionPercent>=percentCut))
  message(paste0("## ",length(geneExpressionPercent)-length(geneToKept)," features in ",assay," assay were removed. ",length(geneToKept)," kept."))
  #can only do this because object is new and no information in data
  object[[assay]]@counts=object[[assay]]@counts[geneToKept,]
  object[[assay]]@data=object[[assay]]@data[geneToKept,]

  #gene processing
  message("## Working on gene expression")
  DefaultAssay(object) <- "Spatial"
  if (normalizeMethod == "SCT") {
    assay <- "SCT"
    object <-
      SCTransform(object, assay = "Spatial", verbose = FALSE,residual.features=customGenes)
    DefaultAssay(object) <- "SCT"
  } else if (normalizeMethod == "log") {
    assay <- "Spatial"
    # object <- NormalizeData(object) %>%
    #   FindVariableFeatures() %>%
    #   ScaleData()
    object <- NormalizeData(object)
    if (!is.null(customGenes)) {
      VariableFeatures(object)=customGenes
    } else {
      object=FindVariableFeatures(object)
    }
    object=ScaleData(object)

  } else {
    stop(paste0("normalizeMethod has to be SCT or log"))
  }
  object = runDimReduc(
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
    DimReducMethod = "PCA",
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
    object = object %>% ScaleData()

    object = runDimReduc(
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
