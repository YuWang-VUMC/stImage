#' integrating and processing for Visium dataset
#'
#' @param object
#' @param normalizeMethod
#' @param ImageFeaturePCs
#' @param pcaDim_s
#' @param pcaDim_i
#' @param geneResolution
#' @param imageFeatureResolution
#' @param imageAndGeneResolution
#' @param cia.nf
#' @param IntNMF.nf
#' @param genePercentCut
#' @param imagePercentCut
#' @param MultimodalMethod
#' @param WnnImageWeightFactor
#'
#' @return
#' @importFrom omicade4 mcia
#' @importFrom IntNMF nmf.mnnals
#' @export
#'
#' @examples
MultimodalImageWorkflowVisium <- function(object, normalizeMethod = "SCT", ImageFeaturePCs = NULL,
                                          pcaDim_s=30, pcaDim_i=30,
                                          geneResolution=0.8, imageFeatureResolution=geneResolution,
                                          imageAndGeneResolution=geneResolution, cia.nf = 20, IntNMF.nf = 20,
                                          genePercentCut=0.1, imagePercentCut=0.3, MultimodalMethod = c("WNN", "MCIA", "IntNMF", "tICA"),
                                          WnnImageWeightFactor=1){
  # analyze spatial matrix
  DefaultAssay(object) <- "Spatial"
  if(normalizeMethod == "SCT") {
    assay <- "SCT"
    object <- SCTransform(object, assay = "Spatial", verbose = FALSE)
    DefaultAssay(object) <- "SCT"
  } else if(normalizeMethod == "log") {
    assay <- "Spatial"
    object <- NormalizeData(object) %>%
      FindVariableFeatures() %>%
      ScaleData()
  }

  #remove less frequent genes
  geneExpressionPercent <- apply(object[[assay]]@counts,1,function(x) length(which(x>0))/length(x))
  VariableFeatures(object) <- setdiff(VariableFeatures(object),names(which(geneExpressionPercent<=genePercentCut)))

  object <- RunPCA(object, assay = assay, npcs = pcaDim_s, verbose = FALSE)
  object <- FindNeighbors(object, reduction = "pca", dims = 1:pcaDim_s)
  object <- FindClusters(object, verbose = FALSE,graph.name = paste0(assay,"_snn"),resolution = geneResolution)
  object <- RunUMAP(object, reduction = "pca", dims = 1:pcaDim_s,
                    assay = assay, reduction.name = 'st.umap', reduction.key = 'stUMAP_')
  # object <- FindSpatiallyVariableFeatures(object, assay = "SCT",
  #                                         features = VariableFeatures(object),
  #                                         selection.method = "markvariogram")
  if(is.null(ImageFeaturePCs)) {
    # analyze imagefeature matrix
    DefaultAssay(object) <- "ImageFeature"
    VariableFeatures(object) <- rownames(object[["ImageFeature"]])
    object <- NormalizeData(object, normalization.method = 'CLR', margin = 1)
    object@assays$ImageFeature@data[is.na(object@assays$ImageFeature@data)] <- 0
    #remove less frequent image features
    imageExpressionPercent <- apply(object[["ImageFeature"]]@counts,1,function(x) length(which(x>0))/length(x))
    VariableFeatures(object) <- setdiff(VariableFeatures(object),names(which(imageExpressionPercent<=imagePercentCut)))


    object <- object %>%
      ScaleData() %>%
      RunPCA(reduction.name = 'ipca', npcs = pcaDim_i, reduction.key = "iPC_")
    object <- FindNeighbors(object, reduction = "ipca", dims = 1:pcaDim_i)
    object <- FindClusters(object, verbose = FALSE,graph.name = "ImageFeature_snn",resolution = imageFeatureResolution)
    object <- RunUMAP(object, reduction = 'ipca', dims = 1:pcaDim_i,
                      assay = 'ImageFeature', reduction.name = 'if.umap', reduction.key = 'ifUMAP_')
  } else {
    DefaultAssay(object) <- "ImageFeature"
    imagefeaturepcs <- ImageFeaturePCs
    object[["ipca"]] <- CreateDimReducObject(embeddings = imagefeaturepcs, key = "iPC_", assay = DefaultAssay(object))
    object <- FindNeighbors(object, reduction = "ipca", dims = 1:pcaDim_i)
    object <- FindClusters(object, verbose = FALSE,graph.name = "ImageFeature_snn",resolution = imageFeatureResolution)
    object <- RunUMAP(object, reduction = 'ipca', dims = 1:pcaDim_i,
                      assay = 'ImageFeature', reduction.name = 'if.umap', reduction.key = 'ifUMAP_')
  }

  if("WNN" %in% MultimodalMethod) {
    # WNN
    object <- FindMultiModalNeighbors(
      object, reduction.list = list("pca", "ipca"),
      dims.list = list(1:pcaDim_s, 1:pcaDim_i),
      modality.weight.name = c("SCT.weight","ImageFeature.weight"),
      return.intermediate = TRUE)

    if (WnnImageWeightFactor!=1) {
      modalityWeightObj=(object@misc$modality.weight)
      modalityWeightObj@modality.weight.list[["ipca"]]=modalityWeightObj@modality.weight.list[["ipca"]]*WnnImageWeightFactor
      object <- FindMultiModalNeighbors(
        object, reduction.list = list("pca", "ipca"),
        dims.list = list(1:pcaDim_s, 1:pcaDim_i), modality.weight.name = c("SCT.weight","ImageFeature.weight"),
        modality.weight=modalityWeightObj)
    }
    object <- FindClusters(object, graph.name = "wsnn", algorithm = 3,resolution = imageAndGeneResolution, verbose = FALSE)
    object <- RunUMAP(object, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
  }
  list <- list()
  if(normalizeMethod == "SCT") {
    list[["SCT"]] <- as.matrix(object@assays$SCT@data)
  } else if (normalizeMethod == "log") {
    list[["Spatial"]] <- as.matrix(object@assays$SCT@data)
  }
  list[["ImageFeature"]] <- object@assays$ImageFeature@data
  list_pos <- list()

  for(j in 1:length(list)) {
    if(min(list[[j]]) < 0) {
      list_pos[[j]] <- list[[j]]+abs(min(list[[j]]))
    } else {
      list_pos[[j]] <- list[[j]]
    }
    list_pos[[j]] <- list_pos[[j]]/max(list_pos[[j]])
  }
  if("MCIA" %in% MultimodalMethod) {
    require(omicade4)
    DefaultAssay(object) <- "SCT"
    mcoin <- mcia(list,cia.nf=cia.nf)
    factor_mcia <- as.matrix(mcoin$mcoa$SynVar)
    rownames(factor_mcia) <- gsub("\\.", "\\-", rownames(factor_mcia))
    metagenes_mcia <- as.matrix(mcoin$mcoa$axis[1:dim(object)[1],])
    object[["mcia"]] <- CreateDimReducObject(embeddings = factor_mcia, key = "mcia_", assay = DefaultAssay(object))
    object <- FindNeighbors(object, reduction = "mcia", graph.name = "mcia_snn", dims = 1:cia.nf)
    object <- FindClusters(object, verbose = FALSE, graph.name = "mcia_snn", resolution = imageAndGeneResolution)
  }
  if("IntNMF" %in% MultimodalMethod) {
    require(IntNMF)
    factorizations_intnmf <- nmf.mnnals(dat=lapply(list_pos, function(x) t(x)), k=IntNMF.nf)
    factors_intNMF <- as.matrix(factorizations_intnmf$W)
    colnames(factors_intNMF) <- paste0("intnmf_", 1:IntNMF.nf)
    DefaultAssay(object) <- "SCT"
    object[["intnmf"]] <- CreateDimReducObject(embeddings = factors_intNMF, key = "intnmf_", assay = DefaultAssay(object))
    object <- FindNeighbors(object, reduction = "intnmf", graph.name = c("intnmf_nn", "intnmf_snn"), dims = 1:IntNMF.nf)
    object <- FindClusters(object, verbose = FALSE, graph.name = "intnmf_snn", resolution = imageAndGeneResolution)
  }
  if("tICA" %in% MultimodalMethod){
    require(tensorBSS)
    omics_tensor <- list()
    for(j in 1:length(list)){
      omics_tensor[[j]] <- cor(list[[j]], method = "spearman")
      if(any(is.na(omics_tensor[[j]]))){
        omics_tensor[[j]][is.na(omics_tensor[[j]])] <- 0
      }
    }
    S <- vector(length = dim(list[[1]])[2]*dim(list[[1]])[2]*length(list))
    dim(S) <- c(length(list), dim(list[[1]])[2], dim(list[[1]])[2])
    for(j in 1:length(list)){
      S[j,,] <- t(omics_tensor[[j]])
    }
    tICA <- DoTICA(S,cia.nf,method="FOBI")
    factor_tica <- as.matrix(tICA$signals)
    colnames(factor_tica) <- paste0("tica_", 1:cia.nf)
    rownames(factor_tica) <- Cells(object)
    object[["tica"]] <- CreateDimReducObject(embeddings = factor_tica, key = "tica_", assay = DefaultAssay(object))
    object <- FindNeighbors(object, reduction = "tica", graph.name = "tica_snn", dims = 1:cia.nf)
    object <- FindClusters(object, verbose = FALSE, graph.name = "tica_snn", resolution = imageAndGeneResolution)
  }
  return(object)
}
