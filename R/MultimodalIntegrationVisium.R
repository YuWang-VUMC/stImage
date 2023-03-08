#' integrating and processing for Visium dataset
#' @inheritParams Seurat
#' @param object A \code{Seurat} object.
#' @param normalizeMethod the normalization method for gene expression matrix.
#' \code{log} for \code{LogNormalize}, \code{SCT} for \code{SCTransform}.
#' @param reduction.list list of dimension reduction matrices for integration
#' @param pcaDim_s number of PCs when running dimension reduction on gene
#' expression (spatial) matrix, 30 as default
#' @param pcaDim_i number of PCs when running dimension reduction on image
#' feature matrix, 30 as default
#' @param resolution resolution of SNN clustering
#' @param nCluster number of clusters
#' @param clusterColumnName column name in meta data to store the clustering
#' results
#' @param cia.nf An integer indicating the number of kept axes in \code{mcia}
#' function
#' @param IntNMF.nf Number of clusters in \code{nmf.mnnals} function
#' @param MultimodalMethod integration methods used for integrating
#' multi-modalites, can be multiple. Note: if include \code{IntNMF}, this step
#' take extremely long time.
#' @param WnnImageWeightFactor numeric value of setting weight of image feature
#' when do WNN integration. 1 as default
#' @param genePercentCut cutoff value of percentage of spots of their
#' values larger than minimum value for filtering genes after normalization.
#' @param imagePercentCut cutoff value of percentage of spots of their
#' values larger than minimum value for prefiltering image features after
#' normalization.
#' @param k.nn the number of multimodal neighbors to compute. 20 by default
#' @param knn.range The number of approximate neighbors to compute
#' @param ... Arguments passed to \code{\link{findSNNClusters}}
#'
#' @return
#' @importFrom omicade4 mcia
#' @importFrom IntNMF nmf.mnnals
#' @importFrom tensorBSS  tensorCentering tPCA tensorTransform tFOBI tJADE
#' @export
#'
#' @examples
#' \dontrun{
#' object <-
#'   MultiModalIntegrationVisium(
#'     object,
#'     MultimodalMethod = "WNN",
#'     pcaDim_s = 20,
#'     pcaDim_i = 20,
#'     reduction.list = list("SCTPCA", "ImageFeaturePCA"),
#'     nCluster = clusterNum,
#'     resolutionMax = 3
#'     )
#' }
#'
MultiModalIntegrationVisium <- function(
    object,
    normalizeMethod = c("SCT", "log"),
    reduction.list = list("SCTPCA", "ImageFeaturePCA"),
    pcaDim_s = 30,
    pcaDim_i = 30,
    resolution = 0.8,
    nCluster = NULL,
    clusterColumnName = if(!is.null(nCluster))
      paste0("wnn_cluster",nCluster) else NULL,
    cia.nf = 20,
    IntNMF.nf = 20,
    genePercentCut=0.05,
    imagePercentCut=0.05,
    MultimodalMethod = c("WNN", "MCIA", "IntNMF", "tICA","Spectrum"),
    WnnImageWeightFactor = 1,
    k.nn = 20,
    knn.range = 200,
    ...) {
  normalizeMethod <- match.arg(normalizeMethod)
  MultimodalMethod <- match.arg(MultimodalMethod,several.ok = TRUE)
  message("Integration by ", paste(MultimodalMethod, collapse = ";"))

  pcaDim_s <- min(pcaDim_s,
                ncol(object@reductions[[reduction.list[[1]]]]@cell.embeddings))
  pcaDim_i <- min(pcaDim_i,
                ncol(object@reductions[[reduction.list[[2]]]]@cell.embeddings))

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
        modalityWeightObj@modality.weight.list[[reduction.list[[2]]]] *
          WnnImageWeightFactor
      object <- FindMultiModalNeighbors(
        object,
        reduction.list = reduction.list,
        dims.list = list(1:pcaDim_s, 1:pcaDim_i),
        knn.graph.name = "wknn",
        snn.graph.name = "wsnn",
        weighted.nn.name = "weighted.nn",
        modality.weight.name = c("SCT.weight", "ImageFeature.weight"),
        k.nn = k.nn,
        knn.range = knn.range,
        modality.weight = modalityWeightObj
      )
    }
    object <- findSNNClusters(object,
                              reduction  = NULL,
                              graph.name = "wsnn",
                              algorithm = 3,
                              nCluster = nCluster,
                              resolution = resolution,
                              clusterColumnName = clusterColumnName,
                              nn.name = "weighted.nn",
                              ...)
  }

  #Spectrum
  if ("Spectrum" %in% MultimodalMethod) {
    object <- doSpectrum(object,
                         reduction.list=reduction.list,
                         nCluster=nCluster)
  }

  list <- list()
  if (normalizeMethod == "SCT") {
    sct.data <- object@assays$SCT@data
    geneExpressionPercent <-
      apply(sct.data, 1, function(x) length(which(x>0)) / length(x))
    filterFeatures <-
      setdiff(rownames(sct.data),
              names(which(geneExpressionPercent <= genePercentCut)))
    sct.data <- sct.data[filterFeatures, ]

    list[["SCT"]] <- as.matrix(sct.data)
  } else if (normalizeMethod == "log") {
    spatial.data <- object@assays$Spatial@data
    geneExpressionPercent <-
      apply(spatial.data, 1, function(x) length(which(x>0)) / length(x))
    filterFeatures <-
      setdiff(rownames(spatial.data),
              names(which(geneExpressionPercent <= genePercentCut)))
    spatial.data <- spatial.data[filterFeatures, ]

    list[["Spatial"]] <- as.matrix(spatial.data)
  }
  if.data <- object@assays$ImageFeature@data
  imagePercent <- apply(if.data, 1, function(x) length(which(x>0)) / length(x))
  filterFeatures <-
    setdiff(rownames(if.data), names(which(imagePercent <= imagePercentCut)))
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

  #MCIA
  if ("MCIA" %in% MultimodalMethod) {
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
    object <- FindNeighbors(
      object,
      reduction = "mcia",
      graph.name = c("mcia_nn", "mcia_snn"),
      dims = 1:cia.nf
    )
    object <- findSNNClusters(
      object,
      reduction  = "mcia",
      graph.name = "mcia_snn",
      algorithm = 3,
      nCluster = nCluster,
      resolution = resolution,
      clusterColumnName = paste0("mcia_snn_cluster", nCluster),
      ...)
  }
  #IntNMF
  if ("IntNMF" %in% MultimodalMethod) {
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
    object <- findSNNClusters(
      object,
      reduction  = "intnmf",
      graph.name = "intnmf_snn",
      algorithm = 3,
      nCluster = nCluster,
      resolution = resolution,
      clusterColumnName = paste0("intnmf_snn_cluster",nCluster),
      ...)
  }
  #tICA
  if ("tICA" %in% MultimodalMethod) {
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
    object <- findSNNClusters(
      object,
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
