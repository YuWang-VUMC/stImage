
doSpectrum=function(dataObj,reduction.list=NULL,graphs=NULL,nCluster=4,clusterColumnName=paste0("Spectrum_cluster",nCluster)) {
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

  #integration similarity matrix
  sIntegrated <- Spectrum::integrate_similarity_matrices(similarityList)
  SpectrumClusterResult <- Spectrum::cluster_similarity(sIntegrated,k=nCluster,clusteralg='GMM')
  object@meta.data[[clusterColumnName]]=SpectrumClusterResult

  return(object)
}
