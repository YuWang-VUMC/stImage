#' FindCoordinateNeighbor
#'
#' @param spotsCoordinate coordinate information of spots
#' @param n number of neighbor spots
#'
#' @return a list containing neighbors of spots
#' @export
#'
FindCoordinateNeighbor <- function(spotsCoordinate,
                                   n=6) {

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


#' ExtractGeneExpByNeighbor
#'
#' @param geneExp gene expression matrix
#' @param spotToNeighborList neighbor spot list returned by
#' \code{FindCoordinateNeighbor}
#' @param spotsImageSimilarity similarity matrix between spots calculated with
#' the neighbor image features
#'
#' @return a gene expression matrix
#' @export
#'
ExtractGeneExpByNeighbor <- function(geneExp,
                                     spotToNeighborList,
                                     spotsImageSimilarity) {
  geneExpByNeighbor <- matrix(0, ncol = ncol(geneExp), nrow = nrow(geneExp))
  colnames(geneExpByNeighbor) <- colnames(geneExp)
  row.names(geneExpByNeighbor) <- row.names(geneExp)

  for (spotName in colnames(geneExpByNeighbor)) {
    spotNeighborNames <- spotToNeighborList[[spotName]]

    spotNeighborSimilarity <- spotsImageSimilarity[spotName, spotNeighborNames]
    spotNeighborSimilaritySum <- sum(spotNeighborSimilarity)

    if (spotNeighborSimilaritySum > 0) {
      spotNeighborSimilarityPropotion <-
        spotNeighborSimilarity / spotNeighborSimilaritySum
      spotNeighborGeneExpBySimilarity <-
        sweep(geneExp[, spotNeighborNames], 2,
              spotNeighborSimilarityPropotion, "*")
      geneExpByNeighbor[, spotName] <-
        rowSums(as.matrix(spotNeighborGeneExpBySimilarity))
    }
  }

  #geneExpByNeighbor==0, keep geneExp.
  #Otherwise take average of geneExp and geneExpByNeighbor
  geneExpNormlizedByNeighbor <-
    (geneExp * (2 - sign(abs(geneExpByNeighbor))) + geneExpByNeighbor) / 2
  geneExpNormlizedByNeighbor[is.na(geneExpNormlizedByNeighbor)] <- 0

  return(geneExpNormlizedByNeighbor)
}


#' dostLearn
#' @param dataObj A \code{Seurat} object
#' @param ipcaObj user provided image feature PCA object
#' @param imageDimReducName image feature dimension reduction matrix in object
#' @param pcaDim_i number of PCs of image feature dimension reduction matrix
#' @param geneDimReducName gene expression dimension reduction matrix in object
#' @param pcaDim_g number of PCs of gene expression dimension reduction matrix
#' @param distanceMethod method for distance calculation, \code{cosine} or
#'\code{euclidean}
#' @param spotsCoordinate spot coordinate information, extract from object if
#' not provided
#' @param geneExp user provided gene expression matrix. i.e. Spark genes
#' @param platform platform information. either "ST" or "Visium"
#' @param normalizePC logical value to turn on normalization of PC matrix.
#' FALSE as default
#' @param assay assay of gene expression
#' @param dataSlot slot of assay of gene expression
#' @param weights define of normalization weights, based on 1)
#' \code{weights_matrix_all} consider gene level correlation in weight;
#' 2) weights from morphological Similarly and physicial distance;
#' 3) physical distance only.
#' @importFrom lsa cosine
#' @importFrom Rfast2 dcora
#' @importFrom Seurat GetAssayData CreateDimReducObject SetAssayData ScaleData
#' @importFrom Seurat FindVariableFeatures DefaultAssay DefaultAssay<-
#' @importFrom dplyr %>%
#' @return a seurat object
#' @export
#'
#' @examples \dontrun{
#' object <- dostLearn(
#'   object,
#'   pcaDim_i = 20,
#'   pcaDim_g = 20,
#'   geneDimReducName = "SCTPCA",
#'   assay = "SCT",
#'   platform = "ST",
#'   distanceMethod = "cosine",
#'   dataSlot = c("data"),
#'   weights="weights_matrix_all")
#' }
#'
dostLearn <- function(
    dataObj,
    ipcaObj = NULL,
    imageDimReducName = "ImageFeaturePCA",
    pcaDim_i = 30,
    geneDimReducName = NULL,
    pcaDim_g = 30,
    distanceMethod = c("cosine","euclidean"),
    spotsCoordinate = NULL,
    platform = ifelse("platform" %in% names(dataObj@misc),
                      dataObj@misc$platform,
                      "Visium"),
    geneExp = NULL,
    normalizePC = FALSE,
    assay = "SCT",
    dataSlot = c("counts","data", "scale.data"),
    weights = c("weights_matrix_all",
              "weights_matrix_pd_md",
              "physical_distance")
) {

  weights <- match.arg(weights)

  #Location distance of spots
  message("find coordinate neighbor ...")
  if (is.null(spotsCoordinate)) {
    spotsCoordinate <- dataObj@images[[1]]@coordinates
    if (("imagecol" %in% colnames(spotsCoordinate)) &
        "imagerow" %in% colnames(spotsCoordinate) ) { #10X Visium data
      #this is too large and slow, still use row and col
      #spotsCoordinate <- spotsCoordinate[,c("imagerow","imagecol")]
      spotsCoordinate <- spotsCoordinate[,c("row","col")]
      platform <- "Visium"
    } else if (("row" %in% colnames(spotsCoordinate)) &
               "col" %in% colnames(spotsCoordinate) ) { #ST
      spotsCoordinate <- spotsCoordinate[,c("row","col")]
      platform <- "ST"
    } else if (("x" %in% colnames(spotsCoordinate)) &
               "y" %in% colnames(spotsCoordinate)) { #don't know platform
      spotsCoordinate <- spotsCoordinate[,c("y","x")]
    } else {
      stop("Can't find imagerow/imagecol or row/col or x/y in spotsCoordinate")
    }
  }

  coordinateNeighborN <- 6
  if (platform == "ST") {
    coordinateNeighborN <- 4
  }
  spotToNeighborList <-
    FindCoordinateNeighbor(spotsCoordinate, n = coordinateNeighborN)

  if (weights=="physical_distance") {
    #only need normlization based on physical_distance
    message("Only use physical_distance for smoothing ...")
    spotsImageSimilarity <-
      matrix(1, ncol=length(spotToNeighborList),
             nrow=length(spotToNeighborList))
    row.names(spotsImageSimilarity) <- names(spotToNeighborList)
    colnames(spotsImageSimilarity) <- names(spotToNeighborList)
  } else { #need normalization based on image or gene similarity
    #distance based on Image PCs
    message("Similarity by image feature PCs ...")
    if (is.null(ipcaObj)) {
      imageFeaturePCs <-
        t(dataObj[[imageDimReducName]]@cell.embeddings[, 1:pcaDim_i])
    } else {
      imageFeaturePCs <- t(ipcaObj@cell.embeddings[, 1:pcaDim_i])
    }

    #distanceMethod <- c("euclidean")
    distanceMethod <- match.arg(distanceMethod)
    if (distanceMethod == "cosine") {
      spotsImageSimilarity <- lsa::cosine(imageFeaturePCs)
    } else if (distanceMethod == "euclidean") {
      spotsImageSimilarity <- 1 - dist(imageFeaturePCs)
    }
    spotsImageSimilarity[spotsImageSimilarity < 0] <- 0

    if (weights == "weights_matrix_all") {
      #consider gene level correlation in weight, see stLearn code
      message("Similarity by gene PCs ...")
      if (is.null(geneDimReducName)) {
        geneDimReducName="SCTPCA"
        warning(paste0("selected weights_matrix_all but geneDimReducName not
                       defined. Use SCTPCA as geneDimReducName"))
      }
      genePCs <- t(dataObj[[geneDimReducName]]@cell.embeddings[, 1:pcaDim_g])
      spotsGeneSimilarity <- Rfast2::dcora(genePCs)
    } else {
      spotsGeneSimilarity <- 1
    }
    spotsImageSimilarity <- spotsImageSimilarity * spotsGeneSimilarity
  }

  #Normalization
  message("Normalization...")
  if (is.null(geneExp)) { #use geneExp from object
    if (normalizePC) { #Normalization PCA
      if (is.null(geneDimReducName)) {
        geneDimReducName="SCTPCA"
        warning(paste0("selected normalizePC but geneDimReducName not defined.
                       Use SCTPCA as geneDimReducName"))
      }
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

  if (dataSlot=="counts") {
    geneExpByNeighbor <- round(geneExpByNeighbor)
    geneExpByNeighbor <- as(geneExpByNeighbor, "dgCMatrix")
  }

  if (normalizePC) { #Need update PCA
    newGeneDimReducName <- paste0(geneDimReducName, "NormalizedByImage")
    dataObj[[newGeneDimReducName]] <-
      CreateDimReducObject(embeddings = t(geneExpByNeighbor),
                           key = "PC_",
                           assay = assay)
    message(paste0("Smoothed gene PCs were added to ", newGeneDimReducName))

  } else { #Need update gene expression
    geneExpByNeighborAssay <- dataObj[[assay]]
    geneExpByNeighborAssay <-
      SetAssayData(object = geneExpByNeighborAssay,
                   slot = dataSlot,
                   new.data = as.matrix(geneExpByNeighbor))
    if (dataSlot=="data") { #need update scale.data and count
      #scale.data slot
      geneExpByNeighborAssay <- FindVariableFeatures(geneExpByNeighborAssay)
      geneExpByNeighborAssay <- geneExpByNeighborAssay %>% ScaleData()
      #counts slot
      geneExpByNeighborAssay <-
        SetAssayData(object = geneExpByNeighborAssay,
        slot = "counts",
        new.data = round(exp(as.matrix(geneExpByNeighbor)) - 1))
    } else if (dataSlot=="counts") {  #need update scale.data and data
      #data slot
      geneExpByNeighborAssay <-
        SetAssayData(object = geneExpByNeighborAssay,
                     slot = "data",
                     new.data = log1p(as.matrix(geneExpByNeighbor)))
      #scale.data slot
      geneExpByNeighborAssay <- FindVariableFeatures(geneExpByNeighborAssay)
      geneExpByNeighborAssay <- geneExpByNeighborAssay %>% ScaleData()
    } else { #dataSlot=="scale.data"
      warning("dataSlot=='scale.data',
              data and counts were NOT normalized by Image")
    }

    newAssayName <- paste0(assay,"NormalizedByImage")
    dataObj[[newAssayName]] <- geneExpByNeighborAssay

    DefaultAssay(dataObj) <- newAssayName
    message(paste0("Gene expression in slot ", dataSlot,
                   " was normlized, stored, and set default assay to ",
                   newAssayName))
  }

  return(dataObj)
}

