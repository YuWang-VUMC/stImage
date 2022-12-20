#' extract image features
#'
#' @param imgFile
#' @param positionTable
#' @param patchSize
#' @param scaleFactor
#' @param rawImage
#' @param savePatchOnImage
#'
#' @return A matrix
#' @importFrom keras application_vgg16
#' @importFrom keras image_load
#' @importFrom keras imagenet_preprocess_input
#' @importFrom plotrix draw.circle
#' @export
#'
#' @examples
extract_features <- function(imgFile,
                             positionTable,
                             patchSize=NULL,
                             scaleFactor=1,
                             saveRawImage = FALSE,
                             savePatchOnImage=NULL) {
  #suppressMessages(require(tensorflow, warn.conflicts = F, quietly = T))
  #suppressMessages(require(keras, warn.conflicts = F, quietly = T))
  #suppressMessages(require(plotrix, warn.conflicts = F, quietly = T))

  patchRadius <- as.integer(patchSize/2)
  model <- application_vgg16(weights = 'imagenet', include_top = FALSE, pooling='avg', input_shape = c(patchSize,patchSize,3))

  img <- image_load(imgFile) %>% image_to_array()
  message(paste0("Loaded image file ", imgFile, " with size: ", paste(dim(img),collapse="X")))

  positionTable <- round(positionTable*scaleFactor)
  positionTable <- positionTable + 1 #positionTable 0 based?
  message(paste0("Loaded position table with ",nrow(positionTable)," spots position"))

  imgPatchAll <- array(dim=c(nrow(positionTable),patchSize,patchSize,3))
  for (i in 1:nrow(positionTable)) {
    imgPatchOne <- img[(positionTable[i,1]-patchRadius):(positionTable[i,1]-patchRadius+patchSize-1),
                    (positionTable[i,2]-patchRadius):(positionTable[i,2]-patchRadius+patchSize-1),]
    imgPatchAll[i,,,] <- imgPatchOne
  }
  message(paste0("Loaded ",nrow(positionTable)," patches with size ",patchSize))

  #RGB
  imgPatch_quantile <- array(dim=c(nrow(positionTable),9,3))
  imgPatch_quantile_mtx <- matrix(NA, nrow = nrow(positionTable), ncol = 27)

  for (i in 1:nrow(positionTable)) {
    for(j in 1:3) {
      imgPatch_q1 <- quantile(imgPatchAll[i,,,j], probs = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9))
      imgPatch_quantile[i,,j] <- imgPatch_q1
    }
    imgPatch_quantile_mtx[i,] <- as.vector(imgPatch_quantile[i,,])
  }
  rownames(imgPatch_quantile_mtx) <- row.names(positionTable)
  colnames(imgPatch_quantile_mtx) <- paste0(rep(c("R", "G", "B"), each = 9), seq(10, 90, 10))


  imgPatchAllProcessed <- imagenet_preprocess_input(imgPatchAll)
  imgPatchFeaturesAll <- model %>% predict(imgPatchAllProcessed)

  row.names(imgPatchFeaturesAll) <- row.names(positionTable)
  colnames(imgPatchFeaturesAll) <- paste0("ImageFeature",1:ncol(imgPatchFeaturesAll))

  if (!is.null(savePatchOnImage)) {
    savePatchOnImageWidth <- dim(img)[1]
    savePatchOnImageHeight <- dim(img)[2]
    if(!saveRawImage){
      saveImageMaxSize <- 2000
      if (any(savePatchOnImageWidth > saveImageMaxSize | savePatchOnImageHeight > saveImageMaxSize)) {
        saveImageSizeFactor <- max(c(savePatchOnImageWidth/saveImageMaxSize,savePatchOnImageHeight/saveImageMaxSize))
        savePatchOnImageWidth <- savePatchOnImageWidth/saveImageSizeFactor
        savePatchOnImageHeight <- savePatchOnImageHeight/saveImageSizeFactor
      }
    }

    png(savePatchOnImage, width=savePatchOnImageWidth, height=savePatchOnImageHeight)
    (img/max(img)) %>% as.raster() %>% plot()

    for (i in 1:nrow(positionTable)) {
      plotrix::draw.circle(x=positionTable[i,2], y=dim(img)[1]-positionTable[i,1], radius=patchSize/2, lwd = 2)
    }
    dev.off()
  }
  if_list <- list("ImageFeature" = imgPatchFeaturesAll, "RGBquantile" = imgPatch_quantile_mtx)
  return(if_list)
}

#' extract image features for Visium dataset
#'
#' @param spatialDir
#' @param patchSize
#' @param rawImage
#' @param imgFile
#' @param ...
#'
#' @return A matrix
#' @importFrom jsonlite fromJSON
#' @export
#'
#' @examples
extract_features_Visium <- function(spatialDir, patchSize=NULL, rawImage = FALSE, imgFile=NULL, ...) {
  suppressMessages(require(jsonlite, warn.conflicts = F, quietly = T))

  positionTableFile <- paste0(spatialDir,"/tissue_positions_list.csv")
  positionTable <- read.csv(positionTableFile,header=FALSE,row.names=1)
  positionTable <- positionTable[which(positionTable[,1]==1),] #spots in tissue only
  positionTable <- positionTable[,c(4:5)] #only need ImageY and ImageX

  scaleFactorFile <- paste0(spatialDir,"/scalefactors_json.json")
  if(rawImage) {
    imgFile=imgFile
    scaleFactor=1
    if(is.null(patchSize)){
      patchSize=fromJSON(scaleFactorFile)$spot_diameter_fullres
    }
  } else {
    imgFile <- paste0(spatialDir,"/tissue_hires_image.png")
    scaleFactor <- fromJSON(scaleFactorFile)$tissue_hires_scalef
  }
  imageFeatures <- extract_features(imgFile, positionTable, patchSize=patchSize, scaleFactor=scaleFactor,  ...)

  return(imageFeatures)
}




