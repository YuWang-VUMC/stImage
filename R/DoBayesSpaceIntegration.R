
DoBayesSpaceIntegration = function(dataObj,
                                   clusterNum=10, BayesSpaceClusterNum = clusterNum,
                                   assay=c("SCT","ImageFeature")
                           ) {
  #BayesSpace
  dataObj = doBayesSpace(
    dataObj,
    assay = assay[1],
    nCluster = BayesSpaceClusterNum,
    name=paste0(assay[1],"BayesSpace")
  )
  #reduction="SCTPCA")  #Not good using these PCs. Need use BayesSpace PCs and variable genes
  dataObj = doBayesSpace(
    dataObj,
    assay = assay[2],
    nCluster = BayesSpaceClusterNum,
    name=paste0(assay[2],"BayesSpace")
  )
  #reduction="ImageFeaturePCA")
  dataObj = IntegrationByDistance(
    dataObj,
    distMatrix =
      paste0(assay,"BayesSpace_Dist"),
      # c("SCTBayesSpace_Dist",
      #   "ImageFeatureBayesSpace_Dist"),
    snnMatrix =
      paste0(assay,"BayesSpace_snn"),
      # c("SCTBayesSpace_snn",
      #   "ImageFeatureBayesSpace_snn"),
    name = c(
      "BayesSpace_DistIntegrated"
    )
  )

  dataObj <-
    findSNNClusters(
      dataObj,
      nCluster = clusterNum,
      graph.name = "BayesSpace_DistIntegrated_snn",
      nn.name = "BayesSpace_DistIntegrated_Neighbor",
      clusterColumnName = paste0("BayesSpace_DistIntegrated_Cluster",clusterNum),
      resolutionMin = 0.01,
      resolutionMax = 4,
      RunUMAP = FALSE
    )
  return(dataObj)
}

