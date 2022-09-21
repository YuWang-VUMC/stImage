
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
  geneExpByNeighbor <- matrix(NA, ncol = ncol(geneExp), nrow = nrow(geneExp))
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
      geneExpByNeighbor[, spotName] <- rowSums(spotNeighborGeneExpBySimilarity)
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
                                 ipcaObj=NULL,
                                 nPCs=30,
                                 distanceMethod="cosine",
                                 spotsCoordinate=NULL,
                                 coordinateNeighborN=4,
                                 geneExp=NULL) {

  #distance based on Image PCs
  message("Similarity by image features...")
  if (is.null(ipcaObj)) {
    imageFeaturePCs <- t(dataObj@reductions$ipca@cell.embeddings[ ,1 : nPCs])
  } else {
    imageFeaturePCs <- t(ipcaObj@cell.embeddings[ ,1 : nPCs])
  }
  #distanceMethod=c("euclidean")
  if (distanceMethod == "cosine") {
    spotsImageSimilarity <- lsa::cosine(imageFeaturePCs)
  } else if (distanceMethod == "euclidean") {
    spotsImageSimilarity <- 1 - dist(imageFeaturePCs)
  }
  spotsImageSimilarity[spotsImageSimilarity<0] <- 0

  #Location distance of spots
  message("find_coordinate_neighbor...")
  if (is.null(spotsCoordinate)) {
    spotsCoordinate <- dataObj@images$images@coordinates
  }
  spotToNeighborList <- find_coordinate_neighbor(spotsCoordinate, n = coordinateNeighborN)

  #Normalization
  message("Normalization...")
  if (is.null(geneExp)) {
    geneExp <- dataObj@assays$SCT@data
  }
  #this is slow
  geneExpByNeighbor <- extract_geneExp_ByNeighbor(geneExp,
                                                  spotToNeighborList,
                                                  spotsImageSimilarity)

  return(geneExpByNeighbor)
}


