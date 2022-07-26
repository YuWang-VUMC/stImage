#load 10X VisiumV1 data
LoadImageFeatureVisium <- function(dataDir, filename = "filtered_feature_bc_matrix.h5",
                                imageFeatures, assay = "Spatial", slice = "slice1",
                                filter.matrix = TRUE, to.upper = FALSE, image = NULL, ...) {
  suppressMessages(require(Seurat, warn.conflicts = F, quietly = T))

  object <- Load10X_Spatial(dataDir, filename = filename, assay = assay,
                            slice = slice, filter.matrix = filter.matrix,
                            to.upper = to.upper, image = image)
  #removing constant features
  constantColumns <- which(apply(imageFeatures, 2, function(x) sd(x, na.rm=TRUE)) == 0)
  if (length(constantColumns)>0) {
    imageFeatures <- imageFeatures[, -constantColumns]
    warning(paste0(length(constantColumns)," image features have SD=0 and were removed"))
  }
  imageFeaturesObj <- CreateAssayObject(counts = t(imageFeatures))
  object[["ImageFeature"]] <- imageFeaturesObj
  message("Image Feature and expression matrix has been successfully loaded!")
  return(object)
}

#load data from raw matrix files
LoadImageFeature <- function(countTable, imageFeatures, positionTable,
                             project = NULL, ...) {
  imageFeaturesObj <- CreateAssayObject(counts = t(imageFeatures))
  object <- CreateSeuratObject(countTable, project = project, assay = "Spatial")
  object[["ImageFeature"]] <- imageFeaturesObj
  object[['images']] <- new(Class = 'SlideSeq', assay = "Spatial", coordinates = positionTable)
  message("Image Feature and expression matrix has been successfully loaded!")
  return(object)
}
