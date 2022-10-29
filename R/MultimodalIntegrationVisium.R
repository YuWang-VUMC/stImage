#' integrating and processing for Visium dataset
#'
#' @param object
#' @param pcaDim_s
#' @param pcaDim_i
#' @param cia.nf
#' @param IntNMF.nf
#' @param MultimodalMethod
#' @param WnnImageWeightFactor
#' @param reduction.list
#' @param resolution
#' @param nCluster
#' @param clusterColumnName
#' @param ...
#' @param normalizeMethod
#'
#' @return
#' @importFrom omicade4 mcia
#' @importFrom IntNMF nmf.mnnals
#' @export
#'
#' @examples
MultiModalIntegrationVisium <- function(object,
                                        normalizeMethod = c("SCT", "log"),
                                        reduction.list = list("SCTPCA", "ImageFeaturePCA"),
                                        pcaDim_s = 30,
                                        pcaDim_i = 30,
                                        resolution = 0.8,
                                        nCluster = NULL,
                                        clusterColumnName = if (!is.null(nCluster)) paste0("wnn_cluster",nCluster) else NULL,
                                        cia.nf = 20,
                                        IntNMF.nf = 20,
                                        genePercentCut=0.05,
                                        imagePercentCut=0.05,
                                        MultimodalMethod = c("WNN", "MCIA", "IntNMF", "tICA","Spectrum"),
                                        WnnImageWeightFactor = 1,
                                        ...) {
  normalizeMethod <- match.arg(normalizeMethod)
  MultimodalMethod=match.arg(MultimodalMethod)
  message("Integration by ",paste(MultimodalMethod,collapse=";"))

  pcaDim_s <- min(pcaDim_s, ncol(object@reductions[[reduction.list[[1]]]]@cell.embeddings))
  pcaDim_i <- min(pcaDim_i, ncol(object@reductions[[reduction.list[[2]]]]@cell.embeddings))

  if ("WNN" %in% MultimodalMethod) {
    # WNN
    object <- FindMultiModalNeighbors(
      object,
      reduction.list = reduction.list,
      dims.list = list(1:pcaDim_s, 1:pcaDim_i),
      knn.graph.name = "wknn",
      snn.graph.name = "wsnn",
      weighted.nn.name = "weighted.nn",
      modality.weight.name = c("SCT.weight", "ImageFeature.weight"),
      return.intermediate = TRUE
    )

    if (WnnImageWeightFactor != 1) {
      modalityWeightObj <- (object@misc$modality.weight)
      modalityWeightObj@modality.weight.list[[reduction.list[[2]]]] <-
        modalityWeightObj@modality.weight.list[[reduction.list[[2]]]] * WnnImageWeightFactor
      object <- FindMultiModalNeighbors(
        object,
        reduction.list = reduction.list,
        dims.list = list(1:pcaDim_s, 1:pcaDim_i),
        knn.graph.name = "wknn",
        snn.graph.name = "wsnn",
        weighted.nn.name = "weighted.nn",
        modality.weight.name = c("SCT.weight", "ImageFeature.weight"),
        modality.weight = modalityWeightObj
      )
    }
    object <- findSNNClusters(object,
                              reduction  = reduction.list[[1]],
                              graph.name = "wsnn",
                              algorithm = 3,
                              nCluster = nCluster,
                              resolution = resolution,
                              clusterColumnName = clusterColumnName,
                              nn.name = "weighted.nn",
                              ...)
    #object <- FindClusters(object, graph.name = "wsnn", algorithm = 3, resolution = resolution, verbose = FALSE)
    #object <- RunUMAP(object, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
  }

  #Spectrum
  if ("Spectrum" %in% MultimodalMethod) {
    s1 <- Spectrum::CNN_kernel(t(dataObj@reductions[[reduction.list[[1]]]]@cell.embeddings))
    s2 <- Spectrum::CNN_kernel(t(dataObj@reductions[[reduction.list[[2]]]]@cell.embeddings))
    sIntegrated <- Spectrum::integrate_similarity_matrices(list(s1,s2))
    SpectrumClusterResult <- Spectrum::cluster_similarity(sIntegrated,k=nCluster,clusteralg='GMM')
    clusterColumnName=paste0("Spectrum_cluster",nCluster)
    object@meta.data[[clusterColumnName]]=SpectrumClusterResult
  }

  list <- list()
  if (normalizeMethod == "SCT") {
    sct.data <- object@assays$SCT@data
    geneExpressionPercent <- apply(sct.data, 1, function(x) length(which(x>0)) / length(x))
    filterFeatures <- setdiff(rownames(sct.data), names(which(geneExpressionPercent <= genePercentCut)))
    sct.data <- sct.data[filterFeatures, ]

    list[["SCT"]] <- as.matrix(sct.data)
  } else if (normalizeMethod == "log") {
    spatial.data <- object@assays$Spatial@data
    geneExpressionPercent <- apply(spatial.data, 1, function(x) length(which(x>0)) / length(x))
    filterFeatures <- setdiff(rownames(spatial.data), names(which(geneExpressionPercent <= genePercentCut)))
    spatial.data <- spatial.data[filterFeatures, ]

    list[["Spatial"]] <- as.matrix(spatial.data)
  }
  if.data <- object@assays$ImageFeature@data
  imagePercent <- apply(if.data, 1, function(x) length(which(x>0)) / length(x))
  filterFeatures <- setdiff(rownames(if.data), names(which(imagePercent <= imagePercentCut)))
  if.data <- if.data[filterFeatures, ]

  list[["ImageFeature"]] <- if.data

  if ("RGB" %in% names(object@assays)) {
    list[["RGB"]] <- as.matrix(object@assays$RGB@data)
  }
  list_pos <- list()

  for (j in 1:length(list)) {
    if (min(list[[j]]) < 0) {
      list_pos[[j]] <- list[[j]] + abs(min(list[[j]]))
    } else {
      list_pos[[j]] <- list[[j]]
    }
    list_pos[[j]] <- list_pos[[j]] / max(list_pos[[j]])
  }
  if ("MCIA" %in% MultimodalMethod) {
    require(omicade4)
    DefaultAssay(object) <- "SCT"
    mcoin <- mcia(list, cia.nf = cia.nf)
    factor_mcia <- as.matrix(mcoin$mcoa$SynVar)
    rownames(factor_mcia) <-
      gsub("\\.", "\\-", rownames(factor_mcia))
    metagenes_mcia <-
      as.matrix(mcoin$mcoa$axis[1:dim(object)[1],])
    object[["mcia"]] <-
      CreateDimReducObject(embeddings = factor_mcia,
                           key = "mcia_",
                           assay = DefaultAssay(object))
    object <-
      FindNeighbors(
        object,
        reduction = "mcia",
        graph.name = c("mcia_nn", "mcia_snn"),
        dims = 1:cia.nf
      )
    object <- findSNNClusters(object,
                              reduction  = "mcia",
                              graph.name = "mcia_snn",
                              algorithm = 3,
                              nCluster = nCluster,
                              resolution = resolution,
                              clusterColumnName = paste0("mcia_snn_cluster",nCluster),
                              ...)
  }
  if ("IntNMF" %in% MultimodalMethod) {
    require(IntNMF)
    factorizations_intnmf <-
      nmf.mnnals(dat = lapply(list_pos, function(x)
        t(x)), k = IntNMF.nf)
    factors_intNMF <- as.matrix(factorizations_intnmf$W)
    colnames(factors_intNMF) <- paste0("intnmf_", 1:IntNMF.nf)
    DefaultAssay(object) <- "SCT"
    object[["intnmf"]] <-
      CreateDimReducObject(embeddings = factors_intNMF,
                           key = "intnmf_",
                           assay = DefaultAssay(object))
    object <-
      FindNeighbors(
        object,
        reduction = "intnmf",
        graph.name = c("intnmf_nn", "intnmf_snn"),
        dims = 1:IntNMF.nf
      )
    object <- findSNNClusters(object,
                              reduction  = "intnmf",
                              graph.name = "intnmf_snn",
                              algorithm = 3,
                              nCluster = nCluster,
                              resolution = resolution,
                              clusterColumnName = paste0("intnmf_snn_cluster",nCluster),
                              ...)
  }
  if ("tICA" %in% MultimodalMethod) {
    require(tensorBSS)
    omics_tensor <- list()
    for (j in 1:length(list)) {
      omics_tensor[[j]] <- cor(list[[j]], method = "spearman")
      if (any(is.na(omics_tensor[[j]]))) {
        omics_tensor[[j]][is.na(omics_tensor[[j]])] <- 0
      }
    }
    S <-
      vector(length = dim(list[[1]])[2] * dim(list[[1]])[2] * length(list))
    dim(S) <-
      c(length(list), dim(list[[1]])[2], dim(list[[1]])[2])
    for (j in 1:length(list)) {
      S[j, ,] <- t(omics_tensor[[j]])
    }
    tICA <- DoTICA(S, cia.nf, method = "FOBI")
    factor_tica <- as.matrix(tICA$signals)
    colnames(factor_tica) <- paste0("tica_", 1:cia.nf)
    rownames(factor_tica) <- Cells(object)
    object[["tica"]] <-
      CreateDimReducObject(embeddings = factor_tica,
                           key = "tica_",
                           assay = DefaultAssay(object))
    object <-
      FindNeighbors(
        object,
        reduction = "tica",
        graph.name = c("tica_nn", "tica_snn"),
        dims = 1:cia.nf
      )
    object <- findSNNClusters(object,
                              reduction  = "tica",
                              graph.name = "tica_snn",
                              algorithm = 3,
                              nCluster = nCluster,
                              resolution = resolution,
                              clusterColumnName = paste0("tica_snn_cluster",nCluster),
                              ...)
  }
  return(object)
}
