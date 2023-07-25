#' extract image features
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Functions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#'
#' @param imgFile the histology image file
#' @param positionTable spots coordinate information
#' @param patchSize the diameter of patches. Suggest 256 for ST platform.
#' For Visium platform, the spot diameter in \code{scalefactors_json.json} is
#' suggested.
#' @param scaleFactor A scaling factor that converts pixel positions in the
#' original, full-resolution image to pixel positions in the image file used
#' here. 1 if used the raw image.
#' @param saveRawImage logical value to define wether save raw image to file.
#' @param savePatchOnImage path and file name of plotting patches on the image.
#'
#' @return A list of two matrices
#' @importFrom keras application_vgg16
#' @importFrom keras image_load
#' @importFrom keras image_to_array
#' @importFrom keras imagenet_preprocess_input
#' @importFrom plotrix draw.circle
#' @importFrom dplyr %>%
#' @importFrom grDevices dev.off as.raster png
#' @examples \dontrun{
#' # ST platform
#' imageFeatures.list <-
#'   ExtractFeatures(
#'     imgFile = rawimg,
#'     positionTable = positionTable,
#'     patchSize = patchsize,
#'     rawImage = T,
#'     savePatchOnImage = paste0("sampleName_rawpatch_", patchsize, ".png")
#'   )
#' imageFeatures <- imageFeatures.list[["ImageFeature"]]
#' RGBquantile <- imageFeatures.list[["RGBquantile"]]
#' }
#'
ExtractFeatures <- function(imgFile,
                             positionTable,
                             patchSize = NULL,
                             scaleFactor = 1,
                             saveRawImage = FALSE,
                             savePatchOnImage = NULL) {

  patchRadius <- as.integer(patchSize/2)
  model <- application_vgg16(weights = 'imagenet',
                             include_top = FALSE,
                             pooling='avg',
                             input_shape = c(patchSize, patchSize, 3))

  img <- image_load(imgFile) %>% image_to_array()
  message(paste0("Loaded image file ", imgFile, " with size: ",
                 paste(dim(img), collapse="X")))

  positionTable <- round(positionTable * scaleFactor)
  positionTable <- positionTable + 1 #positionTable 0 based?
  message(paste0("Loaded position table with ", nrow(positionTable),
                 " spots position"))

  imgPatchAll <- array(dim=c(nrow(positionTable), patchSize, patchSize, 3))
  for (i in 1:nrow(positionTable)) {
    imgPatchOne <- img[(positionTable[i,1] - patchRadius):
                         (positionTable[i,1] - patchRadius + patchSize - 1),
                    (positionTable[i,2] - patchRadius):
                      (positionTable[i,2] - patchRadius + patchSize - 1),]
    imgPatchAll[i,,,] <- imgPatchOne
  }
  message(paste0("Loaded ",nrow(positionTable)," patches with size ",patchSize))

  #RGB
  imgPatch_quantile <- array(dim=c(nrow(positionTable),9,3))
  imgPatch_quantile_mtx <- matrix(NA, nrow = nrow(positionTable), ncol = 27)

  for (i in 1:nrow(positionTable)) {
    for(j in 1:3) {
      imgPatch_q1 <- quantile(imgPatchAll[i,,,j],
                              probs = seq(0.1, 0.9, by = 0.1))
      imgPatch_quantile[i,,j] <- imgPatch_q1
    }
    imgPatch_quantile_mtx[i,] <- as.vector(imgPatch_quantile[i,,])
  }
  rownames(imgPatch_quantile_mtx) <- row.names(positionTable)
  colnames(imgPatch_quantile_mtx) <- paste0(rep(c("R", "G", "B"),
                                                each = 9), seq(10, 90, 10))


  imgPatchAllProcessed <- imagenet_preprocess_input(imgPatchAll)
  imgPatchFeaturesAll <- model %>% predict(imgPatchAllProcessed)

  row.names(imgPatchFeaturesAll) <- row.names(positionTable)
  colnames(imgPatchFeaturesAll) <- paste0("ImageFeature",
                                          1:ncol(imgPatchFeaturesAll))

  if (!is.null(savePatchOnImage)) {
    savePatchOnImageWidth <- dim(img)[1]
    savePatchOnImageHeight <- dim(img)[2]
    if(!saveRawImage){
      saveImageMaxSize <- 2000
      if (any(savePatchOnImageWidth > saveImageMaxSize |
              savePatchOnImageHeight > saveImageMaxSize)) {
        saveImageSizeFactor <- max(c(savePatchOnImageWidth/saveImageMaxSize,
                                     savePatchOnImageHeight/saveImageMaxSize))
        savePatchOnImageWidth <- savePatchOnImageWidth/saveImageSizeFactor
        savePatchOnImageHeight <- savePatchOnImageHeight/saveImageSizeFactor
      }
    }

    png(savePatchOnImage,
        width=savePatchOnImageWidth, height=savePatchOnImageHeight)
    (img/max(img)) %>% as.raster() %>% plot()

    for (i in 1:nrow(positionTable)) {
      plotrix::draw.circle(x=positionTable[i,2],
                           y=dim(img)[1] - positionTable[i,1],
                           radius=patchSize/2, lwd = 2)
    }
    dev.off()
  }
  image_list <- list("ImageFeature" = imgPatchFeaturesAll,
                  "RGBquantile" = imgPatch_quantile_mtx)
  return(image_list)
}

#' extract image features for Visium dataset
#'
#' @param spatialDir the folder contains spatial related files. It's
#' \code{spatial} if data were generated by Space Ranger.
#' @param patchSize the diameter of patches. the spot diameter in
#' \code{scalefactors_json.json} is used if not defined.
#' @param rawImage If or not the raw image is used. FALSE by default.
#' @param imgFile the image file used here if \code{rawImage} is TRUE.
#' @param ... Arguments passed to \code{\link{ExtractFeatures}}
#'
#' @return  A list of two matrices
#' @importFrom jsonlite fromJSON
#' @export
#'
#' @examples \dontrun{
#' # Visium platform
#' imageFeatures.list <-
#'   ExtractFeaturesVisium(
#'     spatial_path,
#'     rawImage = T,
#'     imgFile = rawimg,
#'     savePatchOnImage = NULL)
#' imageFeatures <- imageFeatures.list[["ImageFeature"]]
#' RGBquantile <- imageFeatures.list[["RGBquantile"]]
#' }
ExtractFeaturesVisium <- function(spatialDir,
                                  patchSize=NULL,
                                  rawImage = FALSE,
                                  imgFile=NULL, ...) {
  if(file.exists(paste0(spatialDir,"/tissue_positions_list.csv"))){
    positionTableFile <- paste0(spatialDir,"/tissue_positions_list.csv")
  } else if(file.exists(paste0(spatialDir,"/tissue_positions.csv"))){
    positionTableFile <- paste0(spatialDir,"/tissue_positions.csv")
  } else {
    stop("Please check that if tissue_positions.csv or
         tissue_positions_list.csv exist under spatial folder!\n")
  }
  first_line <- readLines(positionTableFile,n=1)
  if(grepl("^barcode", first_line)){
    positionTable <- read.csv(positionTableFile, header=TRUE, row.names=1)
  }
  else{
    positionTable <- read.csv(positionTableFile, header=FALSE, row.names=1)
  }
  positionTable <- positionTable[which(positionTable[,1]==1),]
  positionTable <- positionTable[,c(4:5)] #only need ImageY and ImageX
  colnames(positionTable) <- c("row", "col")
  scaleFactorFile <- paste0(spatialDir,"/scalefactors_json.json")
  if(rawImage) {
    imgFile=imgFile
    scaleFactor=1
    if(is.null(patchSize)){
      patchSize=fromJSON(scaleFactorFile)$spot_diameter_fullres
      message(paste0("Didn't define patchsize. Set the patch size to ",
                     patchSize, " based on spot_diameter_fullres in ",
                     scaleFactorFile))
    }
  } else {
    imgFile <- paste0(spatialDir,"/tissue_hires_image.png")
    scaleFactor <- fromJSON(scaleFactorFile)$tissue_hires_scalef
    if(is.null(patchSize)){
      patchSize=fromJSON(scaleFactorFile)$spot_diameter_fullres
      patchSize=patchSize*scaleFactor
      message(paste0("Didn't define patchsize. Set the patch size to ",
                     patchSize,
                     " based on spot_diameter_fullres",
                     " X tissue_hires_scalef in ", scaleFactorFile))
    }
  }
  image_list <- ExtractFeatures(imgFile,
                                positionTable,
                                patchSize=patchSize,
                                scaleFactor=scaleFactor,
                                ...)

  return(image_list)
}




