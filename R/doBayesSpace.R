
library(BayesSpace)

doBayesSpace=function(dataObj,platform=c("ST","Visium"),
                      assay=c("SCT","Spatial","ImageFeature"),
                      name=paste0(assay,"BayesSpace_Cluster",nCluster),
                      graphs.name=paste0(name,c("_nn","_snn")),
                      graphNameSimilarity=paste0(name,c("_Similarity")),
                      graphNameDist=paste0(name,c("_Dist")),
                      graphNameChain=paste0(name,c("_ChainMatrix")),
                      nn.name=paste0(name,c("_Neighbor")),
                      clusterColumnName=paste0(name),
                      reduction=NULL,
                      nCluster=4) {

  #sce <- as.SingleCellExperiment(dataObj)

  #assay=match.arg(assay)

  #platform=match.arg(platform)
  colData <- dataObj@images[[1]]@coordinates
  if (("imagecol" %in% colnames(colData)) & "imagerow" %in% colnames(colData) ) { #10X Visum data
    colData <- colData[,c("row","col")]   #still use row, not imagerow for BayesSpace
    platform="Visium"
  } else { #ST data
    platform="ST"
  }
  colnames(colData)=gsub("y","row",colnames(colData))
  colnames(colData)=gsub("x","col",colnames(colData))


  dataTable=as.matrix(dataObj[[assay]]@data)
  sce <- SingleCellExperiment::SingleCellExperiment(assays=list(logcounts=dataTable),
                                                    colData=colData)
  if (!is.null(reduction)) { #defined reduction, use reduction in dataObj
    reducedDims(sce) <- list(PCA=dataObj@reductions[[reduction]]@cell.embeddings)
    sce <- spatialPreprocess(sce, log.normalize=FALSE,platform=platform,skip.PCA = TRUE)
  } else {
    sce <- spatialPreprocess(sce, log.normalize=FALSE,platform=platform)
  }


  sce <- spatialCluster(sce, q=nCluster,
                        #nrep=1000, burn.in=100,
                        nrep=10000, burn.in=500,
                        save.chain=TRUE)
  #dataObj@meta.data[[paste0(assay,"BayesSpace_cluster",nCluster)]]=colData(sce)$spatial.cluster
  dataObj@meta.data[[clusterColumnName]]=colData(sce)$spatial.cluster

  bayesSpaceResultChainMatrix = rhdf5::h5read(metadata(sce)$chain.h5,"z")
  colnames(bayesSpaceResultChainMatrix)=colnames(sce)
  #bayesSpaceResultChainMatrix=apply(bayesSpaceResultChainMatrix, 2, as.character)
  #dataObj@graphs[[paste0(assay,"BayesSpace_cluster",nCluster,"_ChainMatrix")]]=bayesSpaceResultChainMatrix
  dataObj@graphs[[graphNameChain]]=bayesSpaceResultChainMatrix

  bayesSpaceResultChainMatrix

  #bayesSpaceResultDist=cluster::daisy(t(bayesSpaceResultChainMatrix),metric ="gower")

  bayesSpaceResultSimilarity=bayesSpaceChainMatrixSimilarity(bayesSpaceResultChainMatrix)
  colnames(bayesSpaceResultSimilarity)=colnames(dataTable)
  row.names(bayesSpaceResultSimilarity)=colnames(dataTable)
  bayesSpaceResultDist=as.dist(1-bayesSpaceResultSimilarity)

  #graphNameSimilarity=paste0(assay,"BayesSpace_Similarity_Cluster",nCluster)
  #graphNameDist=paste0(assay,"BayesSpace_dist_Cluster",nCluster)
  dataObj@graphs[[graphNameSimilarity]]=bayesSpaceResultSimilarity
  dataObj@graphs[[graphNameDist]]=bayesSpaceResultDist

  #browser()
  prune.SNN=1/100
  dataObj=FindNeighborsByDistance(dataObj,
                                  distMatrix=as.matrix(bayesSpaceResultDist),
                                  assay=assay,
                                  prune.SNN=prune.SNN,
                                  graphs.name=graphs.name,
                                  nn.name=nn.name
                                  )

  # #browser()
  # prune.SNN=1/15
  # bayesSpaceResultNeighborsList=FindNeighbors(bayesSpaceResultDist,
  #                                             prune.SNN =prune.SNN,
  #                                             force.recalc=TRUE)
  # #row.names(bayesSpaceResultNeighborsList[[1]])
  # #browser()
  # graphNameNN=paste0(assay,"BayesSpace_nn_Cluster",nCluster)
  # graphNameSNN=paste0(assay,"BayesSpace_snn_Cluster",nCluster)
  # DefaultAssay(bayesSpaceResultNeighborsList[["nn"]])=assay
  # DefaultAssay(bayesSpaceResultNeighborsList[["snn"]])=assay
  # dataObj[[graphNameNN]]=bayesSpaceResultNeighborsList[["nn"]]
  # dataObj[[graphNameSNN]]=bayesSpaceResultNeighborsList[["snn"]]


  # neighborName=paste0(assay,"BayesSpace_Neightbor_Cluster",nCluster)
  # dataObj@neighbors[[neighborName]]=FindNeighbors(bayesSpaceResultDist,
  #                                                 return.neighbor=TRUE,
  #                                                 prune.SNN =prune.SNN)


#  dataObj <- FindNClusters(dataObj, nCluster = nCluster, graph.name = graphNameSNN,resolutionMax = 4)
  return(dataObj)

}


bayesSpaceChainMatrixSimilarity=function(bayesSpaceResultChainMatrix) {
  nRep=nrow(bayesSpaceResultChainMatrix)
  nSample=ncol(bayesSpaceResultChainMatrix)

  bayesSpaceResultChainMatrixSimilarityValue=combn(ncol(bayesSpaceResultChainMatrix), 2, function(x) sum(bayesSpaceResultChainMatrix[,x[1]]==bayesSpaceResultChainMatrix[,x[2]]))/nRep
  bayesSpaceResultChainMatrixSimilarityValue[bayesSpaceResultChainMatrixSimilarityValue==1]=1-1/nRep

  bayesSpaceResultSimilarityMatrix=matrix(NA,ncol=nSample,nrow=nSample)

  #lower.tri and upper.tri are all column order. So lower.tri is correct order upper.tri is incorrect order
  #use lower.tri first and t(), and use lower.tri again
  bayesSpaceResultSimilarityMatrix[lower.tri(bayesSpaceResultSimilarityMatrix)]=bayesSpaceResultChainMatrixSimilarityValue
  #bayesSpaceResultSimilarityMatrix[upper.tri(bayesSpaceResultSimilarityMatrix)]=bayesSpaceResultChainMatrixSimilarityValue
  bayesSpaceResultSimilarityMatrix=t(bayesSpaceResultSimilarityMatrix)

  bayesSpaceResultSimilarityMatrix[lower.tri(bayesSpaceResultSimilarityMatrix)]=bayesSpaceResultChainMatrixSimilarityValue
  diag(bayesSpaceResultSimilarityMatrix)=1

  return(bayesSpaceResultSimilarityMatrix)

}


