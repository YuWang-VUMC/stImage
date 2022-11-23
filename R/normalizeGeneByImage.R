
#' Title
#'
#' @param spotsCoordinate
#' @param n
#'
#' @return
#' @export
#'
#' @examples
FindCoordinateNeighbor <- function(spotsCoordinate, n=6) {

  spotToNeighborList <- list()
  for (i in 1 : nrow(spotsCoordinate)) {
    distI <- (spotsCoordinate[,1] - spotsCoordinate[i,1])^2 +
      (spotsCoordinate[,2] - spotsCoordinate[i,2])^2
    neighborIndI <- which(rank(distI) <= (n+1) & rank(distI) != 1) #not self
    spotToNeighborList[[i]] <- row.names(spotsCoordinate)[neighborIndI]
  }
  names(spotToNeighborList) <- row.names(spotsCoordinate)
  return(spotToNeighborList)
}


#' Title
#'
#' @param geneExp
#' @param spotToNeighborList
#' @param spotsImageSimilarity
#'
#' @return
#' @export
#'
#' @examples
ExtractGeneExpByNeighbor <- function(geneExp,
                                     spotToNeighborList,
                                     spotsImageSimilarity) {
  geneExpByNeighbor <- matrix(0, ncol = ncol(geneExp), nrow = nrow(geneExp))
  colnames(geneExpByNeighbor) <- colnames(geneExp)
  row.names(geneExpByNeighbor) <- row.names(geneExp)

  for (spotName in colnames(geneExpByNeighbor)) {
    #spotName="AAACAAGTATCTCCCA-1"
    spotNeighborNames <- spotToNeighborList[[spotName]]

    spotNeighborSimilarity <- spotsImageSimilarity[spotName, spotNeighborNames]
    spotNeighborSimilaritySum <- sum(spotNeighborSimilarity)

    if (spotNeighborSimilaritySum > 0) {
      spotNeighborSimilarityPropotion <- spotNeighborSimilarity / spotNeighborSimilaritySum
      spotNeighborGeneExpBySimilarity <- sweep(geneExp[, spotNeighborNames], 2, spotNeighborSimilarityPropotion, "*")
      geneExpByNeighbor[, spotName] <- rowSums(as.matrix(spotNeighborGeneExpBySimilarity))
    }
  }

  #geneExpByNeighbor==0, keep geneExp. Otherwise take average of geneExp and geneExpByNeighbor
  geneExpNormlizedByNeighbor <- (geneExp * (2 - sign(abs(geneExpByNeighbor))) + geneExpByNeighbor) / 2
  geneExpNormlizedByNeighbor[is.na(geneExpNormlizedByNeighbor)] <- 0
  return(geneExpNormlizedByNeighbor)

}


#' Title
#'
#' @param dataObj
#' @param ipcaObj
#' @param nPCs
#' @param distanceMethod
#' @param spotsCoordinate
#' @param coordinateNeighborN
#' @param geneExp
#'
#' @return
#' @export
#'
#' @examples
NormalizeGeneByImage <- function(dataObj,
                                 ipcaObj = NULL,
                                 imageDimReducName = "ImageFeaturePCA",
                                 pcaDim_i = 30,
                                 distanceMethod = "cosine",
                                 spotsCoordinate = NULL,
                                 coordinateNeighborN = 4,
                                 geneExp = NULL,
                                 geneDimReducName = NULL,
                                 assay = "SCT",
                                 dataSlot = c("data", "scale.data")
                                 ) {

  #distance based on Image PCs
  message("Similarity by image features...")
  if (is.null(ipcaObj)) {
    imageFeaturePCs <- t(dataObj[[imageDimReducName]]@cell.embeddings[ ,1:pcaDim_i])
  } else {
    imageFeaturePCs <- t(ipcaObj@cell.embeddings[ , 1:pcaDim_i])
  }
  #distanceMethod=c("euclidean")
  if (distanceMethod == "cosine") {
    spotsImageSimilarity <- lsa::cosine(imageFeaturePCs)
  } else if (distanceMethod == "euclidean") {
    spotsImageSimilarity <- 1 - dist(imageFeaturePCs)
  }
  spotsImageSimilarity[spotsImageSimilarity < 0] <- 0

  #Location distance of spots
  message("find coordinate neighbor...")
  if (is.null(spotsCoordinate)) {
    spotsCoordinate <- dataObj@images[[1]]@coordinates
  }
  spotToNeighborList <- FindCoordinateNeighbor(spotsCoordinate, n = coordinateNeighborN)

  #Normalization
  message("Normalization...")
  if (is.null(geneExp)) { #use geneExp from object
    if (!is.null(geneDimReducName)) { #Normalization PCA
      geneExp <- t(dataObj[[geneDimReducName]]@cell.embeddings)
    } else { #Normalization gene
      dataSlot <- match.arg(dataSlot)
      geneExp <- GetAssayData(dataObj, slot = dataSlot, assay = assay)
    }
  }
  #this is slow
  geneExpByNeighbor <- ExtractGeneExpByNeighbor(geneExp,
                                                spotToNeighborList,
                                                spotsImageSimilarity)

  #return(geneExpByNeighbor)

  if (!is.null(geneDimReducName)) { #Need update PCA
    newGeneDimReducName <- paste0(geneDimReducName, "NormalizedByImage")
    dataObj[[newGeneDimReducName]] <- CreateDimReducObject(embeddings = t(geneExpByNeighbor), key = "PC_", assay = assay)
    message(paste0("Smoothed gene PCs were added to ", newGeneDimReducName))

  } else { #Need update gene expression
    geneExpByNeighborAssay <- dataObj[[assay]]
    geneExpByNeighborAssay <- SetAssayData(object = geneExpByNeighborAssay,
                                           slot = dataSlot,
                                           new.data = geneExpByNeighbor)
    if (dataSlot=="data") { #need update scale.data
      geneExpByNeighborAssay <- FindVariableFeatures(geneExpByNeighborAssay)
      geneExpByNeighborAssay <- geneExpByNeighborAssay %>% ScaleData()
    }
    newAssayName <- paste0(assay,"NormalizedByImage")
    dataObj[[newAssayName]] <- geneExpByNeighborAssay

    DefaultAssay(dataObj) <- newAssayName
    message(paste("Smoothed gene expression was stored and set default assay to ", newAssayName))
  }

  return(dataObj)
}


