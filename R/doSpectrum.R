
doSpectrum=function(dataObj,
                    reduction.list=NULL,graphs=NULL,
                    clusterEachModality=TRUE,
                    nCluster=4,
                    clusterColumnName=paste0("Spectrum_Cluster",nCluster),
                    integrated.name="SpectrumIntegrated_Similarity",
                    nn.name="SpectrumNeighbor",
                    reduction.name = paste0(nn.name,'UMAP'),
                    reduction.key = paste0(nn.name,'UMAP_')) {
  if (is.null(reduction.list) & is.null(graphs)) {
    stop("Need define at least one reduction.list or one graphs")
  }
  #similarity
  if (!is.null(reduction.list)) { #defined reduction.list, not defined graphs. get similarity from reduction.list by Spectrum
    similarityList=NULL
    for (i in 1:length(reduction.list)) {
      similarityList[[reduction.list[[i]]]]= Spectrum::CNN_kernel(t(dataObj@reductions[[reduction.list[[i]]]]@cell.embeddings))
    }
  } else if (!is.null(graphs)) {  #defined graphs. get similarity from graphs
    similarityList=dataObj@graphs[graphs]
  }

  if (clusterEachModality) { #doing single cluster
    for (i in 1:length(similarityList)) {
      clusterResult=Spectrum::cluster_similarity(similarityList[[i]],k=nCluster,clusteralg='GMM')
      if (!is.null(reduction.list)) {
        clusterColumnNameSingle=paste0(reduction.list[[i]],"Spectrum_Cluster",nCluster)
      } else {
        clusterColumnNameSingle=paste0(graphs[i],"Spectrum_Cluster",nCluster)
      }
      dataObj@meta.data[[clusterColumnNameSingle]]=clusterResult
    }
  }

  #integration similarity matrix
  sIntegrated <- Spectrum::integrate_similarity_matrices(similarityList)
  SpectrumClusterResult <- Spectrum::cluster_similarity(sIntegrated,k=nCluster,clusteralg='GMM')
  dataObj@meta.data[[clusterColumnName]]=SpectrumClusterResult

  dataObj@graphs[[integrated.name]]=sIntegrated
  #get NeighborObj for umap
  objNeighbor=distToNeighborObj(1-sIntegrated,k.param=20)
  dataObj@neighbors[[nn.name]]=objNeighbor
  dataObj <- RunUMAP(dataObj,
                     nn.name = nn.name,
                     reduction.name = reduction.name,
                     reduction.key = reduction.key)

  return(dataObj)
}
