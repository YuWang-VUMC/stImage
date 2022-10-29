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
                                 imagePercentCut=0.05,
                                 ...) {
  normalizeMethod <- match.arg(normalizeMethod)
  DimReducMethod <- match.arg(DimReducMethod)

  #gene level
  message("## Working on gene expression")
  DefaultAssay(object) <- "Spatial"
  if (normalizeMethod == "SCT") {
    assay <- "SCT"
    object <-
      SCTransform(object, assay = "Spatial", verbose = FALSE)
    DefaultAssay(object) <- "SCT"
  } else if (normalizeMethod == "log") {
    assay <- "Spatial"
    object <- NormalizeData(object) %>%
      FindVariableFeatures() %>%
      ScaleData()
  } else {
    stop(paste0("normalizeMethod has to be SCT or log"))
  }
  object = runDimReduc(
    object,
    DimReducMethod = DimReducMethod,
    assay = assay,
    percentCut = genePercentCut,
    pcaDim = pcaDim_s,
    ...
  )

  #image level
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

  #RGB features
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
