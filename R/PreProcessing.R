#' MultiModalPreProcessing
#' @include RunDimReduc.R
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
#' @param ... Arguments passed to \code{\link{RunDimReduc}}
#'
#' @return  A \code{Seurat} object
#' @importFrom SpatialPCA CreateSpatialPCAObject
#' @importFrom Seurat SCTransform DefaultAssay NormalizeData ScaleData
#' @importFrom Seurat VariableFeatures<-
#' @importFrom Seurat FindVariableFeatures VariableFeatures DefaultAssay<-
#' @importFrom dplyr %>%
#' @export
#'
#' @examples \dontrun{
#' object <-
#'   PreProcessing(
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
PreProcessing <- function(object,
                          normalizeMethod = c("SCT", "log"),
                          pcaDim_s = 30,
                          pcaDim_i = 30,
                          pcaDim_c = 6,
                          DimReducMethod = "PCA",
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
  DefaultAssay(object) <- assay
  percentCut <- prefiltergenePercentCut
  assayData_counts <- GetAssayData(object = object, assay = assay, slot = "counts")
  minValue <- min(assayData_counts)
  if (minValue != 0) {
    warning(paste0("Min value in ", assay, " is not equal to 0.
                   Need to confirm data filtering by >=", minValue,
                   " percent is correct!"))
  }
  geneExpressionPercent <- apply(assayData_counts, 1, function(x)
    length(which(x > minValue)) / length(x))
  geneToKept <- names(which(geneExpressionPercent >= percentCut))
  message(paste0("## ", length(geneExpressionPercent) - length(geneToKept),
                 " features in ", assay, " assay were removed. ",
                 length(geneToKept), " kept."))
  #can only do this because object is new and no information in data
  newdata_counts <- assayData_counts[geneToKept,]
  #object <- SetAssayData(object = object, new.data = newdata_counts)
  new_object <- CreateSeuratObject(newdata_counts, project = object@project.name, assay = assay, )

  ###ImageFeature
  assay <- "ImageFeature"
  DefaultAssay(object) <- assay
  percentCut <- prefilterimagePercentCut
  assayData_counts <- GetAssayData(object = object, assay = assay, slot = "counts")
  minValue <- min(assayData_counts)
  if (minValue != 0) {
    warning(paste0("Min value in ", assay, " is not equal to 0.
                   Need to confirm data filtering by >=", minValue,
                   " percent is correct!"))
  }
  geneExpressionPercent <- apply(assayData_counts, 1, function(x)
    length(which(x>minValue)) / length(x))
  geneToKept <- names(which(geneExpressionPercent >= percentCut))
  message(paste0("## ", length(geneExpressionPercent) - length(geneToKept),
                 " features in ", assay, " assay were removed. ",
                 length(geneToKept), " kept."))
  #can only do this because object is new and no information in data
  newdata_counts <- assayData_counts[geneToKept,]
  new_img_obj <- CreateAssayObject(newdata_counts)
  #object <- SetAssayData(object = object, slot = "counts", new.data = newdata_counts)
  new_object[[assay]] <- new_img_obj

  if("RGB" %in% names(object@assays)){
    assay <- "RGB"
    new_rgb_obj <- CreateAssayObject(object[["RGB"]]$counts, key = "rgb_")
    new_object@assays[["RGB"]] <- new_rgb_obj
  }
  new_object@images <- object@images

  #gene processing
  message("## Working on gene expression")
  DefaultAssay(new_object) <- "Spatial"
  if (normalizeMethod == "SCT") {
    assay <- "SCT"
    new_object <-
      SCTransform(new_object,
                  assay = "Spatial",
                  new.assay.name = assay,
                  verbose = FALSE,
                  residual.features = customGenes)
    DefaultAssay(new_object) <- "SCT"
  } else if (normalizeMethod == "log") {
    assay <- "Spatial"
    new_object <- NormalizeData(new_object)
    if (!is.null(customGenes)) {
      VariableFeatures(new_object) <- customGenes
    } else {
      new_object <- FindVariableFeatures(new_object)
    }
    new_object <- ScaleData(new_object)

  } else {
    stop(paste0("normalizeMethod has to be SCT or log"))
  }
  new_object <- RunDimReduc(
    new_object,
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
  DefaultAssay(new_object) <- assay
  VariableFeatures(new_object) <- rownames(new_object[[assay]])
  new_object <-
    NormalizeData(new_object, normalization.method = 'CLR', margin = 1)
  assayData_data <- GetAssayData(object = new_object, assay = assay, slot = "data")
  assayData_data[is.na(assayData_data)] <- 0
  new_object <- SetAssayData(object = new_object, slot = "data", new.data = assayData_data)
  new_object <- new_object %>% ScaleData()

  new_object <- RunDimReduc(
    new_object,
    DimReducMethod = DimReducMethod,
    assay = assay,
    percentCut = imagePercentCut,
    pcaDim = pcaDim_i,
    ...
  )

  #RGB features processing
  if ("RGB" %in% Seurat::Assays(new_object)) {
    message("## Working on RGB features")

    assay <- "RGB"
    DefaultAssay(new_object) <- assay
    VariableFeatures(new_object) <- row.names(new_object[[assay]])
    assayData_data <- GetAssayData(object = new_object, assay = assay, slot = "data")
    rgb_norm <- 1 - (assayData_data / 255)
    new_object <- SetAssayData(object = new_object, slot = "data", new.data = as(object = rgb_norm, Class = 'dgCMatrix'))
    new_object <- new_object %>% ScaleData()

    new_object <- RunDimReduc(
      new_object,
      DimReducMethod = "PCA",
      assay = assay,
      percentCut = imagePercentCut,
      pcaDim = pcaDim_c,
      ...
    )
  }

  return(new_object)
}
