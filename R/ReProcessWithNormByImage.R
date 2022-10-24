#' Title
#'
#' @param object
#' @param normalizeMethod
#' @param renormbyimage
#' @param pcaDim_s
#' @param pcaDim_i
#' @param ipcaObj
#' @param distanceMethod
#' @param spotsCoordinate
#' @param coordinateNeighborN
#' @param geneExp
#' @param DimReducMethod
#' @param genePercentCut
#' @param imagePercentCut
#' @param geneResolution
#' @param imageFeatureResolution
#' @param sparkversion
#' @param numCores_spark
#' @param gene.number
#' @param customGenelist
#' @param min.loctions
#' @param min.features
#'
#' @return
#' @export
#'
#' @examples
ReProcessWithNormByImage <- function(object,
                                        normalizeMethod = c("SCT", "log"),
                                        renormbyimage = c("GeneExp", "GeneExpPCs"),
                                        pcaDim_s = 30,
                                        pcaDim_i = 30,
                                        ipcaObj = NULL,
                                        distanceMethod = "cosine",
                                        spotsCoordinate = NULL,
                                        coordinateNeighborN = 6,
                                        geneExp = NULL,
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
  if(renormbyimage == "GeneExp") {
    if(normalizeMethod == "SCT") {
      DefaultAssay(object) <- "SCT"
      assay <- "SCT"
      renorm_gene_exp <- normalizeGeneByImage(object,
                                              ipcaObj = ipcaObj,
                                              nPCs = pcaDim_s,
                                              distanceMethod = distanceMethod,
                                              coordinateNeighborN = coordinateNeighborN)
      object@assays$SCT@data <- renorm_gene_exp
      object <- ScaleData(object)
    } else if(normalizeMethod == "log") {
      DefaultAssay(object) <- "Spatial"
      assay <- "Spatial"
      renorm_gene_exp <- normalizeGeneByImage(object,
                                              ipcaObj = ipcaObj,
                                              nPCs = pcaDim_s,
                                              distanceMethod = distanceMethod,
                                              geneExp = object@assays$Spatial@data,
                                              coordinateNeighborN = coordinateNeighborN)
      object@assays$Spatial@data <- renorm_gene_exp
      object <- ScaleData()
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
  } else if(renormbyimage == "GeneExpPCs") {
    if(normalizeMethod == "SCT") {
      assay <- "SCT"
      DefaultAssay(object) <- "SCT"
    } else if(normalizeMethod == "log") {
      assay <- "Spatial"
      DefaultAssay(object) <- "Spatial"
    }

    genePCs <- t(object@reductions$pca@cell.embeddings[,1:pcaDim_s])
    norm_genePCs_by_image <- normalizeGeneByImage(object,
                                                  ipcaObj = ipcaObj,
                                                  nPCs = pcaDim_s,
                                                  distanceMethod = distanceMethod,
                                                  coordinateNeighborN = coordinateNeighborN,
                                                  geneExp = genePCs)


    object[[DimReducMethod]]=CreateDimReducObject(embeddings = t(norm_genePCs_by_image), key = "PC_", assay = DefaultAssay(object))
    object=findSNNClusters(object,reduction = DimReducMethod,
                           graph.name = paste0(DimReducMethod,c("_nn", "_snn")),
                           dims = 1:pcaDim_s,
                           resolution=geneResolution,
                           assay = DefaultAssay(object),
                           reduction.name = paste0(DimReducMethod,c("_umap")),
                           reduction.key = paste0(DimReducMethod,c('UMAP_')))



    if(DimReducMethod == "PCA") {
      object[["pca"]] <- CreateDimReducObject(embeddings = t(norm_genePCs_by_image), key = "PC_", assay = DefaultAssay(object))
      object <- FindNeighbors(object, reduction = "pca", dims = 1:pcaDim_s)
      object <- FindClusters(object, verbose = FALSE, graph.name = paste0(assay,"_snn"), resolution = geneResolution)
      object <- RunUMAP(object, reduction = 'pca', dims = 1:pcaDim_s,
                        assay = assay, reduction.name = 'st.umap', reduction.key = 'stUMAP_')
    } else if(DimReducMethod == "SpatialPCA") {
      object[["GeneSpatialPCA"]] <- CreateDimReducObject(embeddings = t(norm_genePCs_by_image), key = "GeneSpatial_", assay = DefaultAssay(object))
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
}
