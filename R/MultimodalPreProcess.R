#' Title
#'
#' @param object
#' @param GeneFeaturePCobj
#' @param ImageFeaturePCobj
#' @param normalizeMethod
#' @param pcaDim_s
#' @param pcaDim_i
#' @param genePercentCut
#' @param imagePercentCut
#' @param geneResolution
#' @param imageFeatureResolution
#'
#' @return
#' @export
#'
#' @examples
MultimodalPreProcess <- function(object, GeneFeaturePCobj = NULL, ImageFeaturePCobj = NULL,
                                 normalizeMethod = c("SCT", "log"),
                                 pcaDim_s=30, pcaDim_i=30,
                                 genePercentCut=0.1, imagePercentCut=0.3,
                                 geneResolution=0.8, imageFeatureResolution=geneResolution){
  if(is.null(GeneFeaturePCobj)) {
    # analyze spatial matrix
    DefaultAssay(object) <- "Spatial"
    #rawcount <- object[[assay]]@counts
    #positionTable <- object@images$images@coordinates
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
    object <- FindClusters(object, verbose = FALSE, graph.name = paste0(assay,"_snn"), resolution = geneResolution)
    object <- RunUMAP(object, reduction = "pca", dims = 1:pcaDim_s,
                      assay = assay, reduction.name = 'st.umap', reduction.key = 'stUMAP_')
  } else {
    if(normalizeMethod == "SCT") {
      assay <- "SCT"
      DefaultAssay(object) <- "SCT"
    } else if(normalizeMethod == "log") {
      assay <- "Spatial"
      DefaultAssay(object) <- "Spatial"
    }
    #object[["pca"]] <- CreateDimReducObject(embeddings = GeneFeaturePCs, key = "PC_", assay = DefaultAssay(object))
    object[["pca"]] <- GeneFeaturePCobj
    object <- FindNeighbors(object, reduction = "pca", dims = 1:pcaDim_s)
    object <- FindClusters(object, verbose = FALSE, graph.name = paste0(assay,"_snn"), resolution = geneResolution)
    object <- RunUMAP(object, reduction = 'pca', dims = 1:pcaDim_s,
                      assay = assay, reduction.name = 'st.umap', reduction.key = 'stUMAP_')
  }

  if(is.null(ImageFeaturePCobj)) {
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
    #object[["ipca"]] <- CreateDimReducObject(embeddings = ImageFeaturePCs, key = "iPC_", assay = DefaultAssay(object))
    object[["ipca"]] <- ImageFeaturePCobj
    object <- FindNeighbors(object, reduction = "ipca", dims = 1:pcaDim_i)
    object <- FindClusters(object, verbose = FALSE,graph.name = "ImageFeature_snn",resolution = imageFeatureResolution)
    object <- RunUMAP(object, reduction = 'ipca', dims = 1:pcaDim_i,
                      assay = 'ImageFeature', reduction.name = 'if.umap', reduction.key = 'ifUMAP_')
  }
}

#Example of preparing DimReducObject of SpatialPCA results

#if(pcaMethod == "SpatialPCA"){
#  suppressMessages(require(SpatialPCA, warn.conflicts = F, quietly = T))

#    spcaObj <- CreateSpatialPCAObject(counts=rawcount, location=positionTable,
#                                 project = "SpatialPCA", gene.type="spatial",
#                                 sparkversion="spark", gene.number=3000,
#                                 customGenelist=NULL,min.loctions = 20, min.features=20)
#
#    spcaObj <- SpatialPCA_buildKernel(spcaObj, kerneltype="gaussian", bandwidthtype="SJ")
#    spcaObj <- SpatialPCA_EstimateLoading(spcaObj, fast=FALSE, SpatialPCnum=20)
#    spcaObj <- SpatialPCA_SpatialPCs(spcaObj, fast=FALSE)

#    object[["pca"]] <- CreateDimReducObject(embeddings = t(SpatialPCs@SpatialPCs),
#                                            loadings = SpatialPCs@W,
#                                            key = "PC_", assay = DefaultAssay(object))
#}
