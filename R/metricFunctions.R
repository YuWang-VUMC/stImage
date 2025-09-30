

#Average Spot distance to neighbor, comparing with average spot distance to random spots
#to indicate if the assumption holds (spots are similar to their neighbors)
#large (positive) means spots are similar to their neighbors than random spots. Should use spatial method
#small (negative) means spots are NOT similar to their neighbors than random spots. Should NOT use spatial method
#' Title
#'
#' @param x result table from extractNeighborDist
#' @param varName1 var for Neighbor distance. default is Neighbor.mean
#' @param varName2 var for distance to compare with. default is Random.mean
#' @param main
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
AverageSimilarityPlot <- function(x,
                          varName1 = "Neighbor.mean",
                          varName2 = "Random.mean",
                          main=paste0("Difference between ", varName1," and ", varName2),
                          ...) {
  print("Propotion of positive and negative difference")
  print(table(sign(x[,varName1] - x[,varName2])) / nrow(x))
  
  # hist(x[,varName1]-x[,varName2],main=main,xlab="",...)
  # abline(v=0)
  
  x$PlotValue <- x[,varName1] - x[,varName2]
  x$NeighborVsRandomSpots <- ifelse(x$PlotValue>0, "No", "Yes")
  
  p <- ggplot(x, aes(x = PlotValue, fill = NeighborVsRandomSpots)) +
    geom_histogram()+
    xlab("Neighbors vs. Random Spots") +
    theme_bw()
  return(p)
}



#extract distance between Neighbors, and other spots
#' Title
#'
#' @param dataObj
#' @param neighborsK number of neighbors to consider
#'
#' @return
#' @export
#'
#' @examples
extractNeighborDist <- function(dataObj,
                                neighborsK = 10,
                                sampleN=20000,
                                pcaN =1:20,
                                reduction="SCTPCA") {
  if (ncol(dataObj)>sampleN) {
    h=sample(1:ncol(dataObj),sampleN)
  } else {
    h=1:ncol(dataObj)
  }
  
  #correlation of neighbors
  #source("~/source/stImage/R/normalizeGeneByImage.R")
  if (class(dataObj@images[[1]])=="FOV") {
    spotsCoordinate <- dataObj@images[[1]]@boundaries$centroids@coords
    row.names(spotsCoordinate)=row.names((dataObj@meta.data))
  } else {
    spotsCoordinate <- dataObj@images[[1]]@coordinates
  }
  if (("imagecol" %in% colnames(spotsCoordinate)) & "imagerow" %in% colnames(spotsCoordinate) ) { #10X Visium data
    #this is too large and slow, still use row and col
    #spotsCoordinate <- spotsCoordinate[,c("imagerow","imagecol")]
    spotsCoordinate <- spotsCoordinate[,c("row","col")]
    platform <- "Visium"
  } else if (("row" %in% colnames(spotsCoordinate)) & "col" %in% colnames(spotsCoordinate) ) { #ST
    spotsCoordinate <- spotsCoordinate[,c("row","col")]
    platform <- "ST"
  } else if (("x" %in% colnames(spotsCoordinate)) & "y" %in% colnames(spotsCoordinate)) { #other platforms
    spotsCoordinate <- spotsCoordinate[,c("y","x")]
  } else if (("coor_x" %in% colnames(spotsCoordinate)) & "coor_y" %in% colnames(spotsCoordinate)) { #other platforms
    spotsCoordinate <- spotsCoordinate[,c("coor_y","coor_x")]
  } else {
    stop("Can't find imagerow/imagecol or row/col or x/y in spotsCoordinate")
  }
  
  FindCoordinateNeighbor= function(spotsCoordinate,index=1 : nrow(spotsCoordinate), n=6) {
    spotToNeighborList <- list()
    for (i in 1:length(index)) {
      distI <- (spotsCoordinate[,1] - spotsCoordinate[index[i],1])^2 +
        (spotsCoordinate[,2] - spotsCoordinate[index[i],2])^2
      neighborIndI <- which(rank(distI) <= (n+1) & rank(distI) != 1) #not self
      spotToNeighborList[[i]] <- row.names(spotsCoordinate)[neighborIndI]
    }
    names(spotToNeighborList) <- row.names(spotsCoordinate)[index]
    return(spotToNeighborList)
  }
  spotToNeighbors <- FindCoordinateNeighbor(spotsCoordinate, index=h,n = neighborsK)
  includedSpotNames=na.omit(unique(c(unname(unlist(spotToNeighbors)),names(spotToNeighbors))))
  
  #distance
  distByGene <- as.matrix(dist(dataObj@reductions[[reduction]]@cell.embeddings[includedSpotNames,pcaN]))
  
  #number of neighbor spots in 20 nearest spots
  
  distToNeighborsOutAll <- NULL
  distToAllSpotsOutAll <- NULL
  distToRandomSpotsOutAll <- NULL
  for (j in 1:length(spotToNeighbors)) {
    distToNeighbors <- distByGene[names(spotToNeighbors)[j], spotToNeighbors[[j]]]
    distToNeighborsOut <- c(min(distToNeighbors), max(distToNeighbors),
                            median(distToNeighbors), mean(distToNeighbors),
                            var(distToNeighbors), sd(distToNeighbors),
                            quantile(distToNeighbors, c(0.025, 0.05, 0.1, 0.25, 0.75, 0.95)))
    distToNeighborsOutAll <- rbind(distToNeighborsOutAll, distToNeighborsOut)
    colnames(distToNeighborsOutAll) <- paste0("Neighbor.", c("min","max","median","mean","var","sd","Q2.5","Q5","Q10","Q25","Q75","Q95"))
    
    #distToAllSpots=distByGene[names(spotToNeighbors)[j],]
    distToAllSpots <- distByGene[names(spotToNeighbors)[j], setdiff(colnames(distByGene), spotToNeighbors[[j]])]
    distToAllSpots <- distToAllSpots[distToAllSpots != 0]
    distToAllSpotsOut <- c(min(distToAllSpots), max(distToAllSpots),
                           median(distToAllSpots), mean(distToAllSpots),
                           var(distToAllSpots), sd(distToAllSpots),
                           quantile(distToAllSpots,c(0.025, 0.05, 0.1, 0.25, 0.75, 0.95)),
                           mean(sort(distToAllSpots)[1:neighborsK]))
    distToAllSpotsOutAll <- rbind(distToAllSpotsOutAll, distToAllSpotsOut)
    colnames(distToAllSpotsOutAll) <- paste0("All.", c("min","max","median","mean","var","sd","Q2.5","Q5","Q10","Q25","Q75","Q95","Nearest.mean"))
    
    distToRandomSpots <- sample(distToAllSpots,neighborsK)
    distToRandomSpots <- distToRandomSpots[distToRandomSpots != 0]
    distToRandomSpotsOut <- c(min(distToRandomSpots), max(distToRandomSpots),
                              median(distToRandomSpots), mean(distToRandomSpots),
                              var(distToRandomSpots), sd(distToRandomSpots),
                              quantile(distToRandomSpots, c(0.025, 0.05, 0.1, 0.25, 0.75, 0.95)))
    distToRandomSpotsOutAll <- rbind(distToRandomSpotsOutAll, distToRandomSpotsOut)
    colnames(distToRandomSpotsOutAll) <- paste0("Random.", c("min","max","median","mean","var","sd","Q2.5","Q5","Q10","Q25","Q75","Q95"))
    
  }
  row.names(distToNeighborsOutAll) <- names(spotToNeighbors)
  return(data.frame(distToNeighborsOutAll, distToAllSpotsOutAll, distToRandomSpotsOutAll))
}

#plot distance consistency between two modality embeddings
#' Title
#'
#' @param dataObj
#' @param i number of PCs to use
#' @param sampleN number of samples in resampling
#' @param normalize normalize distance by min distance, True or False
#' @param diffCutoff cutoff value to highlight spot pairs 
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
ModalityConcordancePlot <- function(dataObj,
                                i=1:20,
                                sampleN=20000,
                                normalize=TRUE,
                                diffCutoff=5,
                                geneReduction="SCTPCA",
                                ImageReduction="ImageFeaturePCA",
                                #highlight.cutoff=0.3,
                                ...) {
  #distance of spots
  temp1 <- dist(dataObj@reductions[[ImageReduction]]@cell.embeddings[,i])
  temp2 <- dist(dataObj@reductions[[geneReduction]]@cell.embeddings[,i])
  j <- sample(1:length(as.vector(temp1)), sampleN)
  if (normalize) {
    temp1ToPlot <- as.vector(temp1)[j] / min(as.vector(temp1)[j])
    temp2ToPlot <- as.vector(temp2)[j] / min(as.vector(temp2)[j])
    #temp1ToPlot=(as.vector(temp1)[j]-min(as.vector(temp1)[j]))/(range(as.vector(temp1))[2]-range(as.vector(temp1))[1])
    #temp2ToPlot=(as.vector(temp2)[j]-min(as.vector(temp2)[j]))/(range(as.vector(temp2))[2]-range(as.vector(temp2))[1])
  } else {
    temp1ToPlot <- as.vector(temp1)[j]
    temp2ToPlot <- as.vector(temp2)[j]
  }
  #plot(temp1ToPlot,temp2ToPlot,xlab="Distance between spots (Image)",ylab="Distance between spots (Gene)",pch=16,...)
#  dataForPlot <- data.frame(Image = temp1ToPlot,
#                            Gene = temp2ToPlot)
  dataForPlot <- data.frame(Image = temp1ToPlot,
                            Gene = temp2ToPlot
                            #Issue1=ifelse(abs(temp1ToPlot-temp2ToPlot)>=highlight.cutoff,"Yes","No"),
                            #Issue2=ifelse(abs(log(temp1ToPlot/temp2ToPlot))>=abs(log(highlight.cutoff)),"Yes","No")
                            )
  p <- ggplot(dataForPlot,aes(x=Image,y=Gene)) +
    geom_point( alpha = 0.6) +
    #geom_point(aes(colour=Issue1), alpha = 0.3) +
    stat_density_2d(aes(fill = after_stat(level)), geom = "polygon", alpha = 0.6) +
    #geom_density_2d_filled(contour_var = "density",alpha = 0.3)+
    xlab("Image Distance") +
    ylab("Gene Distance") +
    theme_bw()
  
  return(p)
  #print(cor(as.vector(temp1)[j],as.vector(temp2)[j]))
}

#Max Spot distance to neighbor, comparing with min spot distance to all non-neighbor spots
#to indicate if the some spots are more similar to their non-neighbors
#large (positive) means spot not similar to neighbor spots are similar to non-neighbor spots. Should NOT use spatial method
#small (negative) means spot not similar to neighbor spots are also NOT similar to non-neighbor spots. Should use spatial method
#' Title
#'
#' @param x result table from extractNeighborDist
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
HighSimilarityNonNeighborPlot <- function(x, var.x="Neighbor.min",var.x.ref="All.Q5", 
                                          var.y="Neighbor.median",var.y.ref="All.Q5",
                                          var.x.highlight=0.5,
                                          var.y.highlight=0.5, ...) {
  temp1=x[, c("Neighbor.max")] / x[,c("All.Q25")]
  temp2=x[,var.y]/x[,var.y.ref]
  temp3=x[,var.x]/x[,var.x.ref]
  
  temp1 <- log2(temp1)
  temp2 <- log2(temp2)
  temp3 <- log2(temp3)
  #plot(temp3,temp1)
  
  dataForPlot <- data.frame(temp1, temp2, temp3)
  dataForPlot$MoreSimilarityToNonNeighbor <- ifelse(temp3 > var.x.highlight| temp2 > var.y.highlight, "SimilarityToNonNeighbor", "Uncertain")
  #dataForPlot$PlotValueCategory <- ifelse(temp3 > 0, "No", "Yes")
  print(table(dataForPlot$MoreSimilarityToNonNeighbor)/nrow(dataForPlot))
  p <- ggplot(dataForPlot, aes(x = temp3, y = temp2, colour = MoreSimilarityToNonNeighbor)) +
  #p <- ggplot(dataForPlot, aes(x = temp3, y = temp2)) +
    geom_point() +
    xlab("Most Similar Neighbor vs. Most Similar Non-Neighbors") +
    ylab("Median Neighbors vs. Most Similar Non-Neighbors") +
    theme_bw()
  return(p)
  
  abline(v = 0)
  abline(v = c(median(temp1), median(temp2), median(temp3)), lty = 2)
  # print("Median Difference:")
  # print(median(temp1)-median(temp2))
  #print(mean(temp1)-mean(temp2))
  print("Median value of Distance to Neighbor(Max)")
  print(median(temp1))
  
  print("Median value of Distance to Neighbor(Mean)")
  print(median(temp2))
  
  print("Median value of Distance to Neighbor(Min)")
  print(median(temp3))
  #abline(h=1)
}
