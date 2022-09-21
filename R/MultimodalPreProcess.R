#' Title
#'
#' @param object
#' @param GeneFeaturePCs
#' @param ImageFeaturePCs
#' @param normalizeMethod
#' @param pcaDim_s
#' @param pcaDim_i
#' @param genePercentCut
#' @param imagePercentCut
#' @param geneResolution
#' @param imageFeatureResolution
#' @param DimReducMethod
#' @param sparkversion
#' @param numCores_spark
#' @param gene.number
#' @param customGenelist
#' @param min.loctions
#' @param min.features
#'
#' @return
#' @importFrom SpatialPCA CreateSpatialPCAObject
#' @export
#'
#' @examples
MultimodalPreProcess <- function(object,
                                 GeneFeaturePCs = NULL,
                                 ImageFeaturePCs = NULL,
                                 normalizeMethod = c("SCT", "log"),
                                 pcaDim_s = 30,
                                 pcaDim_i = 30,
                                 DimReducMethod = c("PCA", "SpatialPCA"),
                                 genePercentCut = 0.1,
                                 imagePercentCut = 0.3,
                                 geneResolution = 0.8,
                                 imageFeatureResolution = geneResolution,
                                 sparkversion = "spark",
                                 numCores_spark = 1,
                                 gene.number = 3000,
                                 customGenelist = NULL,
                                 min.loctions = 30,
                                 min.features = 30) {
  location <- as.matrix(object@images$images@coordinates)
  if(is.null(GeneFeaturePCs)) {
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
    if(DimReducMethod == "PCA") {
      #remove less frequent genes
      geneExpressionPercent <- apply(object[[assay]]@counts, 1, function(x) length(which(x>0)) / length(x))
      VariableFeatures(object) <- setdiff(VariableFeatures(object), names(which(geneExpressionPercent <= genePercentCut)))

      object <- RunPCA(object, assay = assay, npcs = pcaDim_s, verbose = FALSE)
      object <- FindNeighbors(object, reduction = "pca", dims = 1:pcaDim_s)
      object <- FindClusters(object, verbose = FALSE, graph.name = paste0(assay,"_snn"), resolution = geneResolution)
      object <- RunUMAP(object,
                        reduction = "pca",
                        dims = 1:pcaDim_s,
                        assay = assay,
                        reduction.name = 'st.umap',
                        reduction.key = 'stUMAP_')
    } else if(DimReducMethod == "SpatialPCA") {
      rawcount <- object@assays[["Spatial"]]@counts
      SPCAobj <- CreateSpatialPCAObject(counts = rawcount,
                                        location = location,
                                        project = object@project.name,
                                        gene.type = "spatial",
                                        sparkversion = sparkversion,
                                        gene.number = gene.number,
                                        customGenelist = customGenelist,
                                        min.loctions = min.loctions,
                                        min.features = min.features)

      SPCAobj <- SpatialPCAWorkflow(SPCAobj)
      SPCApcs <- t(SPCAobj@SpatialPCs)
      colnames(SPCApcs) <- paste0("PC",1:ncol(SPCApcs))
      object[["GeneSpatialPCA"]] <- CreateDimReducObject(embeddings = SPCApcs, key = "GeneSpatial_", assay = DefaultAssay(object))
      object <- FindNeighbors(object, reduction = "GeneSpatialPCA", graph.name = c("st_spca_nn", "st_spca_snn"), dims = 1:ncol(SPCApcs))
      object <- FindClusters(object, verbose = FALSE, graph.name = "st_spca_snn", resolution = geneResolution)
      object <- RunUMAP(object,
                        reduction = 'GeneSpatialPCA',
                        dims = 1:ncol(SPCApcs),
                        assay = DefaultAssay(object),
                        reduction.name = 'st.spca.umap',
                        reduction.key = 'stspcaUMAP_')
    }
  } else {
    if(normalizeMethod == "SCT") {
      assay <- "SCT"
      DefaultAssay(object) <- "SCT"
    } else if(normalizeMethod == "log") {
      assay <- "Spatial"
      DefaultAssay(object) <- "Spatial"
    }
    if(DimReducMethod == "PCA") {
      object[["pca"]] <- CreateDimReducObject(embeddings = GeneFeaturePCs, key = "PC_", assay = DefaultAssay(object))
      object <- FindNeighbors(object, reduction = "pca", dims = 1:pcaDim_s)
      object <- FindClusters(object, verbose = FALSE, graph.name = paste0(assay,"_snn"), resolution = geneResolution)
      object <- RunUMAP(object, reduction = 'pca', dims = 1:pcaDim_s,
                        assay = assay, reduction.name = 'st.umap', reduction.key = 'stUMAP_')
    } else if(DimReducMethod == "SpatialPCA") {
      object[["GeneSpatialPCA"]] <- CreateDimReducObject(embeddings = GeneFeaturePCs, key = "GeneSpatial_", assay = DefaultAssay(object))
      object <- FindNeighbors(object, reduction = "GeneSpatialPCA", graph.name = c("st_spca_nn", "st_spca_snn"), dims = 1:ncol(GeneFeaturePCs))
      object <- FindClusters(object, verbose = FALSE, graph.name = "st_spca_snn", resolution = geneResolution)
      object <- RunUMAP(object,
                        reduction = 'GeneSpatialPCA',
                        dims = 1:ncol(GeneFeaturePCs),
                        assay = DefaultAssay(object),
                        reduction.name = 'st.spca.umap',
                        reduction.key = 'stspcaUMAP_')
    }
  }
  if(is.null(ImageFeaturePCs)) {
    # analyze imagefeature matrix
    DefaultAssay(object) <- "ImageFeature"
    VariableFeatures(object) <- rownames(object[["ImageFeature"]])
    object <- NormalizeData(object, normalization.method = 'CLR', margin = 1)
    object@assays$ImageFeature@data[is.na(object@assays$ImageFeature@data)] <- 0
    #remove less frequent image features
    imageExpressionPercent <- apply(object[["ImageFeature"]]@counts, 1, function(x) length(which(x>0)) / length(x))
    VariableFeatures(object) <- setdiff(VariableFeatures(object), names(which(imageExpressionPercent <= imagePercentCut)))
    if(DimReducMethod == "PCA") {
      object <- object %>%
      ScaleData() %>%
      RunPCA(assay = DefaultAssay(object), reduction.name = 'ipca', npcs = pcaDim_i, reduction.key = "iPC_")
      object <- FindNeighbors(object, reduction = "ipca", dims = 1:pcaDim_i)
      object <- FindClusters(object, verbose = FALSE, graph.name = "ImageFeature_snn", resolution = imageFeatureResolution)
      object <- RunUMAP(object, reduction = 'ipca',
                        dims = 1:pcaDim_i,
                        assay = DefaultAssay(object),
                        reduction.name = 'if.umap',
                        reduction.key = 'ifUMAP_')
    } else if(DimReducMethod == "SpatialPCA") {
      imageFeaturesNormlized <- object@assays[["ImageFeature"]]@data
      IFSPCAobj <- CreateSpatialPCAObject(counts = imageFeaturesNormlized,
                                          location = location,
                                          project = object@project.name,
                                          gene.type="custom",
                                          customGenelist = row.names(imageFeaturesNormlized),
                                          min.loctions = min.loctions,
                                          min.features = min.features)

      IFSPCAobj <- SpatialPCAWorkflow(IFSPCAobj)
      IFSPCApcs <- t(IFSPCAobj@SpatialPCs)
      colnames(IFSPCApcs) <- paste0("PC",1:ncol(IFSPCApcs))
      object[["ImageSpatialPCA"]] <- CreateDimReducObject(embeddings = IFSPCApcs, key = "ImageSpatial_", assay = DefaultAssay(object))
      object <- FindNeighbors(object, reduction = "ImageSpatialPCA", graph.name = c("if_spca_nn", "if_spca_snn"), dims = 1:ncol(IFSPCApcs))
      object <- FindClusters(object, verbose = FALSE, graph.name = "if_spca_snn", resolution = imageFeatureResolution)
      object <- RunUMAP(object,
                        reduction = 'ImageSpatialPCA',
                        dims = 1:ncol(IFSPCApcs),
                        assay = DefaultAssay(object),
                        reduction.name = 'if.spca.umap',
                        reduction.key = 'ifspcaUMAP_')
    }
  } else {
    DefaultAssay(object) <- "ImageFeature"
    if(DimReducMethod == "PCA") {
      object[["ipca"]] <- CreateDimReducObject(embeddings = ImageFeaturePCs, key = "iPC_", assay = DefaultAssay(object))
      object <- FindNeighbors(object, reduction = "ipca", dims = 1:ncol(ImageFeaturePCs))
      object <- FindClusters(object, verbose = FALSE,graph.name = "ImageFeature_snn",resolution = imageFeatureResolution)
      object <- RunUMAP(object,
                        reduction = 'ipca',
                        dims = 1:ncol(ImageFeaturePCs),
                        assay = DefaultAssay(object),
                        reduction.name = 'if.umap',
                        reduction.key = 'ifUMAP_')
    } else if(DimReducMethod == "SpatialPCA"){
      object[["ImageSpatialPCA"]] <- CreateDimReducObject(embeddings = ImageFeaturePCs, key = "ImageSpatial_", assay = DefaultAssay(object))
      object <- FindNeighbors(object, reduction = "ImageSpatialPCA", graph.name = c("if_spca_nn", "if_spca_snn"), dims = 1:ncol(ImageFeaturePCs))
      object <- FindClusters(object, verbose = FALSE, graph.name = "if_spca_snn", resolution = imageFeatureResolution)
      object <- RunUMAP(object,
                        reduction = 'ImageSpatialPCA',
                        dims = 1:ncol(ImageFeaturePCs),
                        assay = DefaultAssay(object),
                        reduction.name = 'if.spca.umap',
                        reduction.key = 'ifspcaUMAP_')
    }

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
