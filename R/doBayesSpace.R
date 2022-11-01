
library(BayesSpace)

doBayesSpace=function(dataObj,platform="ST",assay=c("SCT","Spatial","ImageFeature"),nCluster=4) {

  #sce <- as.SingleCellExperiment(dataObj)

  assay=match.arg(assay)

  colData <- dataObj@images$images@coordinates
  colnames(colData)=gsub("y","row",colnames(colData))
  colnames(colData)=gsub("x","col",colnames(colData))

  dataTable=as.matrix(dataObj[[assay]]@data)
  sce <- SingleCellExperiment::SingleCellExperiment(assays=list(logcounts=dataTable),
                                                    colData=colData)
  sce <- spatialPreprocess(sce, log.normalize=FALSE,platform=platform)

  sce <- spatialCluster(sce, q=nCluster,
                        nrep=1000, burn.in=100,
                        save.chain=TRUE)
  dataObj@meta.data[[paste0(assay,"BayesSpace_cluster",nCluster)]]=colData(sce)$spatial.cluster

  bayesSpaceResultChainMatrix = rhdf5::h5read(metadata(sce)$chain.h5,"z")
  colnames(bayesSpaceResultChainMatrix)=colnames(sce)
  #bayesSpaceResultChainMatrix=apply(bayesSpaceResultChainMatrix, 2, as.character)
  dataObj@graphs[[paste0("BayesSpace_cluster",nCluster,"_ChainMatrix")]]=bayesSpaceResultChainMatrix

  #bayesSpaceResultDist=cluster::daisy(t(bayesSpaceResultChainMatrix),metric ="gower")

  bayesSpaceResultSimilarity=bayesSpaceChainMatrixSimilarity(bayesSpaceResultChainMatrix)
  colnames(bayesSpaceResultSimilarity)=colnames(dataTable)
  row.names(bayesSpaceResultSimilarity)=colnames(dataTable)
  bayesSpaceResultDist=1-bayesSpaceResultSimilarity

  graphNameSimilarity=paste0(assay,"BayesSpace_Similarity_Cluster",nCluster)
  graphNameDist=paste0(assay,"BayesSpace_dist_Cluster",nCluster)
  dataObj@graphs[[graphNameSimilarity]]=bayesSpaceResultSimilarity
  dataObj@graphs[[graphNameDist]]=bayesSpaceResultDist


  bayesSpaceResultNeighborsList=FindNeighbors(bayesSpaceResultDist)
  #row.names(bayesSpaceResultNeighborsList[[1]])

  graphNameNN=paste0(assay,"BayesSpace_nn_Cluster",nCluster)
  graphNameSNN=paste0(assay,"BayesSpace_snn_Cluster",nCluster)
  dataObj@graphs[[graphNameNN]]=bayesSpaceResultNeighborsList[["nn"]]
  dataObj@graphs[[graphNameSNN]]=bayesSpaceResultNeighborsList[["snn"]]

#  dataObj <- FindNClusters(dataObj, nCluster = nCluster, graph.name = graphNameSNN,resolutionMax = 4)
  return(dataObj)

}


bayesSpaceChainMatrixSimilarity=function(bayesSpaceResultChainMatrix) {
  nRep=nrow(bayesSpaceResultChainMatrix)
  nSample=ncol(bayesSpaceResultChainMatrix)

  bayesSpaceResultChainMatrixSimilarityValue=combn(ncol(bayesSpaceResultChainMatrix), 2, function(x) sum(bayesSpaceResultChainMatrix[,x[1]]==bayesSpaceResultChainMatrix[,x[2]]))/nRep
  bayesSpaceResultChainMatrixSimilarityValue[bayesSpaceResultChainMatrixSimilarityValue==1]=1-1/nRep

  bayesSpaceResultSimilarityMatrix=matrix(NA,ncol=nSample,nrow=nSample)
  bayesSpaceResultSimilarityMatrix[lower.tri(bayesSpaceResultSimilarityMatrix)]=bayesSpaceResultChainMatrixSimilarityValue
  bayesSpaceResultSimilarityMatrix[upper.tri(bayesSpaceResultSimilarityMatrix)]=bayesSpaceResultChainMatrixSimilarityValue
  diag(bayesSpaceResultSimilarityMatrix)=1

  return(bayesSpaceResultSimilarityMatrix)

}


