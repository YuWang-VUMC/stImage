#' gg_color_hue
#'
#' @param n number of colors
#' @importFrom grDevices hcl
#' @return a color palette
#' @export
#'
#' @examples \dontrun{
#'   gg_color_hue(8)
#' }
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 95)[1:n]
}

#' st_voronoi_point
#'
#' @param points sf data frame of points
#' @importFrom sf st_combine st_voronoi st_collection_extract
#' @importFrom sf st_geometry_type st_geometry st_intersects
#' @return a sf data frame
#' @export
#'
#' @examples \dontrun{
#'   location <- GetTissueCoordinates(object)
#'   p <- st_as_sf(data.frame(x=location$y, y=dplyr::desc(location$x),
#'   A=1:nrow(location)), coords=1:2)
#'   v <- st_voronoi_point(p)
#' }
st_voronoi_point <- function(points){
  ## points must be POINT geometry
  # check for point geometry and execute if true
  if(!all(st_geometry_type(points) == "POINT")){
    stop("Input not  POINT geometries")
  }
  g = st_combine(st_geometry(points)) # make multipoint
  v = st_voronoi(g)
  v = st_collection_extract(v)
  return(v[unlist(st_intersects(points, v))])
}

#' AllSpatialDimPlot
#' @param object a seurat object
#' @param cluster.frame a logical value to define if show the cluster frame
#' @param cluster.highlight the cluster id to show if \code{cluster.frame} is
#' TRUE
#' @param col a color palette
#' @param groups the colname of clusters in metadata of object
#' @param ... Arguments passed to \code{\link{ggplot}}
#' @importFrom sf st_as_sf st_convex_hull st_union
#' @importFrom sf st_intersection st_cast
#' @importFrom ggplot2 ggplot geom_sf geom_point aes_string scale_color_manual
#' @importFrom ggplot2 theme_void guide_legend guides ggtitle theme element_text
#' @importFrom Seurat GetTissueCoordinates Cells
#' @return a list of ggplot objects
#' @export
#'
#' @examples \dontrun{
#'   plots <- AllSpatialDimPlot(object,
#'   cluster.frame = F,
#'   col = colors)
#' }
AllSpatialDimPlot <- function(object,
                              cluster.frame = TRUE,
                              cluster.highlight = NULL,
                              col = NULL,
                              groups=NULL,
                              ...) {

  location <- GetTissueCoordinates(object)[,c(1,2)]
  if (identical(colnames(location), c("imagerow","imagecol"))) {
    #Visum, change to Y and X to match codes in other parts
    colnames(location)=c("x","y")
  }

  if (is.null(groups)) { #groups to plot NOT defined. Plot all groups in data
    geneClusterName <-
      grep("SCTPCA_cluster",
           colnames(object@meta.data),
           value=TRUE)
    renormgeneClusterName <-
      grep("SCTNormalizedByImagePCA_cluster",
           colnames(object@meta.data),
           value=TRUE)
    geneBayesSpaceClusterName <-
      grep("SCTBayesSpace_Cluster",
           colnames(object@meta.data),
           value=TRUE)
    ImageFeatureBayesSpaceClusterName <-
      grep("ImageFeatureBayesSpace_Cluster",
           colnames(object@meta.data),
           value=TRUE)
    BayesSpacewnnClusterName <-
      grep("BayesSpace_DistIntegrated_Cluster",
           colnames(object@meta.data),
           value=TRUE)
    geneSpectrumClusterName <-
      grep("SCTPCASpectrum_Cluster",
           colnames(object@meta.data),
           value=TRUE)
    ImageFeatureSpectrumClusterName <-
      grep("ImageFeaturePCASpectrum_Cluster",
           colnames(object@meta.data),
           value=TRUE)
    RGBSpectrumClusterName <-
      grep("RGBPCASpectrum_Cluster",
           colnames(object@meta.data),
           value=TRUE)
    SpectrumClusterName <-
      grep("Spectrum_Cluster",
           colnames(object@meta.data),
           value=TRUE)
    geneSPCASpectrumClusterName <-
      grep("SCTSpatialPCASpectrum_Cluster",
           colnames(object@meta.data),
           value=TRUE)
    ImageFeatureSPCASpectrumClusterName <-
      grep("ImageFeatureSpatialPCASpectrum_Cluster",
           colnames(object@meta.data),
           value=TRUE)
    geneSPCAClusterName <-
      grep("SCTSpatialPCA_cluster",
           colnames(object@meta.data),
           value=TRUE)
    imageClusterName <-
      grep("ImageFeaturePCA_cluster",
           colnames(object@meta.data),
           value=TRUE)
    imageFeatureSPCAClusterName <-
      grep("ImageFeatureSpatialPCA_cluster",
           colnames(object@meta.data),
           value=TRUE)
    RGBClusterName <-
      grep("RGBPCA_cluster",
           colnames(object@meta.data),
           value=TRUE)
    wnnClusterName <-
      grep("wnn_cluster",
           colnames(object@meta.data),
           value=TRUE)
    mciaClusterName <-
      grep("mcia_snn_cluster",
           colnames(object@meta.data),
           value=TRUE)
    intnmfClusterName <-
      grep("intnmf_snn_cluster",
           colnames(object@meta.data),
           value=TRUE)
    ticaClusterName <-
      grep("tica_snn_cluster",
           colnames(object@meta.data),
           value=TRUE)

    groups <- c(geneClusterName, renormgeneClusterName,
                geneBayesSpaceClusterName, ImageFeatureBayesSpaceClusterName,
                BayesSpacewnnClusterName, geneSpectrumClusterName,
                ImageFeatureSpectrumClusterName, RGBSpectrumClusterName,
                SpectrumClusterName, geneSPCASpectrumClusterName,
                ImageFeatureSPCASpectrumClusterName, geneSPCAClusterName,
                imageClusterName, imageFeatureSPCAClusterName,
                RGBClusterName, wnnClusterName,
                mciaClusterName, intnmfClusterName, ticaClusterName)
    groups <- groups[!is.na(groups)]
  } else { #groups to plot defined.
    groups <- intersect(groups, colnames(object@meta.data))
  }

  plots <- list()
  for (i in 1:length(groups)) {

    cl <- setNames(object@meta.data[, groups[i]],
                   rownames(object@meta.data[groups[i]]))
    location$cluster <- as.factor(cl)
    # arrange colors by Palo if not set by user

    #pal <- gg_color_hue(length(unique(cl)))
    #palopal <- Palo(location, cl, pal, init_iter = 100,
    #                refine_iter = 200,
    #                early_stop = 50)
#
    #if(is.null(col)){
    #  col <- palopal
    #}
#
    if(cluster.frame){
      if(is.null(cluster.highlight)) {
        stop("Must provide the cluster ids when setting cluster.frame == True!")
      } else {
        # get cluster frame with sf
        p <- st_as_sf(data.frame(x=location$y, y=dplyr::desc(location$x),
                                 A=1:nrow(location)), coords=1:2)
        v <- st_voronoi_point(p)
        hull <- st_convex_hull(st_union(p))

        vcut <- st_intersection(st_cast(v), hull)
        sf <- st_as_sf(vcut, celltype = cl)
        sf_0 <- sf[sf$celltype %in% as.character(cluster.highlight),]

        plots[[groups[i]]] <- ggplot(data = location) +
          geom_sf(data = st_union(sf_0),
                  color = col[as.character(cluster.highlight)]) +
          geom_point(aes_string(x = location$y,
                                y = dplyr::desc(location$x),
                                col = location$cluster), size = 0.5) +
          scale_color_manual(values = col) +
          ggtitle(groups[i]) +
          guides(color = guide_legend(title = groups[i])) +
          theme_void() +
          theme(legend.position="none",
                plot.title = element_text(hjust = 0.5))
      }
    } else {
      plots[[groups[i]]] <- ggplot(data = location) +
        #geom_sf(data = st_union(sf_0),
        #  color = col[as.character(cluster.highlight)]) +
        geom_point(aes_string(x = location$y,
                              y = dplyr::desc(location$x),
                              col = location$cluster), size = 0.5) +
        scale_color_manual(values = col) +
        ggtitle(groups[i]) +
        guides(color = guide_legend(title = groups[i])) +
        theme_void() +
        theme(legend.position="none",
              plot.title = element_text(hjust = 0.5))
    }
  }
  return(plots)
}

#' AllSpatialFeaturePlot
#'
#' @param object a seurat object
#' @param normalizeMethod parameter for choosing the assay for gene expression
#' visualization. Can be either "SCT" or "log"
#' @param integreation.method the integration method for the final clustering.
#' related to the cluster frame if \code{cluster.frame} is TRUE
#' @param features a vector to define features want to plot
#' @param cluster.frame a logical value to define if show the cluster frame
#' @param cluster.highlight the cluster id to show if \code{cluster.frame} is
#' TRUE
#' @param col a color palette
#' @param ... Arguments passed to \code{\link{ggplot}}
#' @importFrom sf st_as_sf st_convex_hull st_union
#' @importFrom sf st_intersection st_cast
#' @importFrom ggplot2 ggplot geom_sf geom_point aes_string scale_color_manual
#' @importFrom ggplot2 scale_color_gradient2 scale_alpha
#' @importFrom ggplot2 theme_void guide_legend guides ggtitle theme element_text
#' @importFrom Seurat GetTissueCoordinates ScaleData
#' @return a list of ggplot objects
#' @export
#'
#' @examples \dontrun{
#'   plots <- AllSpatialFeaturePlot(object,
#'   normalizeMethod = "SCT",
#'   integreation.method = "wnn",
#'   features = c("TM4SF1", "S100A4"),
#'   cluster.frame = TRUE,
#'   col = colors,
#'   cluster.highlight = 1)
#' }
AllSpatialFeaturePlot <- function(
    object,
    normalizeMethod = c("SCT", "log"),
    slot="data",
    cl_id=NULL,
    pt.size=0.5,
    line.size=0.5,
    integreation.method = c("wnn", "mcia", "intnmf", "tica"),
    features = NULL,
    cluster.frame = TRUE,
    cluster.highlight = NULL,
    col = NULL,
    ...) {

  if(is.null(features)){
    stop("Must provide the features you want to plot!")
  } else {
    if (normalizeMethod == "SCT") {
      #exp <- object@assays$SCT@data
      assay <- "SCT"
      #exp <- slot(object@assays$SCT,x)
      exp <- GetAssayData(object = object, assay = assay, slot = slot)
    } else if(normalizeMethod == "log"){
      #exp <- object@assays$Spatial@data
      #exp <- slot(object@assays$Spatial,x)
      exp <- GetAssayData(object = object, assay = "Spatial", slot = slot)
    } else {
      stop("Please tell which slots of assays will be used to generate feature
           plots!")
    }
    exp <- as.matrix(exp)
    scale_exp <- ScaleData(exp, features = row.names(exp))

    location <- GetTissueCoordinates(object)[,c(1,2)]
    if (identical(colnames(location), c("imagerow","imagecol"))) {
      #Visum, change to Y and X to match codes in other parts
      colnames(location)=c("x","y")
    }

    wnnClusterName <-
      grep("wnn_cluster",
           colnames(object@meta.data),
           value=TRUE)
    mciaClusterName <-
      grep("mcia_snn_cluster",
           colnames(object@meta.data),
           value=TRUE)
    intnmfClusterName <-
      grep("intnmf_snn_cluster",
           colnames(object@meta.data),
           value=TRUE)
    ticaClusterName <-
      grep("tica_snn_cluster",
           colnames(object@meta.data),
           value=TRUE)

    if (is.null(cl_id)) {
      #cl_id="SCTPCA_cluster6_renamed"
      if(integreation.method == "wnn") {
        cl_id <- wnnClusterName
      } else if(integreation.method == "mcia") {
        cl_id <- mciaClusterName
      } else if(integreation.method == "intnmf") {
        cl_id <- intnmfClusterName
      } else if(integreation.method == "tica") {
        cl_id <- ticaClusterName
      }
    }



    plots <- list()
    for (i in 1:length(features)) {
      cl <- setNames(object@meta.data[,cl_id],
                     rownames(object@meta.data[cl_id]))
      scale_features <- t(scale_exp[features[i], , drop = F])
      scale_features <- cbind(location, scale_features)
      scale_features$alpha <-
        (scale_features[, features[i]] +
           abs(min(scale_features[, features[i]]))) /
        (abs(min(scale_features[, features[i]])) +
           max(scale_features[, features[i]]))
      ## arrange colors by Palo if not set by user
      #pal <- gg_color_hue(length(unique(cl)))
      #palopal <- Palo(location, cl, pal)
      #if(is.null(col)){
      #  col <- palopal
      #}

      if(cluster.frame){
        if(is.null(cluster.highlight)) {
          stop("Must provide the cluster ids when setting cluster.frame ==
               True!")
        } else {
          # get cluster frame with sf
          p <- st_as_sf(data.frame(x=location$y,
                                   y=dplyr::desc(location$x),
                                   A=1:nrow(location)), coords=1:2)
          v <- st_voronoi_point(p)
          hull <- st_convex_hull(st_union(p))

          vcut <- st_intersection(st_cast(v), hull)
          sf <- st_as_sf(vcut, celltype = cl)
          sf_0 <- sf[sf$celltype %in% as.character(cluster.highlight),]

          plots[[features[i]]] <- ggplot(scale_features) +
            geom_sf(data = st_union(sf_0),
                    color = "black", size = line.size,linewidth =line.size, fill = NA) +
            geom_point(aes_string(x = scale_features$y,
                                  y = dplyr::desc(scale_features$x),
                                  col = features[i],
                                  fill = features[i],
                                  alpha = "alpha"),
                       size = pt.size) +
            #scale_color_viridis(option="viridis") +
            ggtitle(features[i]) +
            #scale_color_gradient2(palette = "PiYG") +
            scale_color_gradient2(low = 'green', mid = 'white', high = 'red') +
            scale_alpha(range = c(0.2, 1)) +
            guides(color = guide_legend(title = features[i])) +
            theme_void() +
            theme(legend.position="none",
                  plot.title = element_text(hjust = 0.5))
        }
      } else {
        plots[[features[i]]] <- ggplot(scale_features) +
          #geom_sf(data = st_union(sf_0), color = col[cluster.highlight]) +
          geom_point(aes_string(x = scale_features$y,
                                y = dplyr::desc(scale_features$x),
                                col = features[i],
                                fill = features[i],
                                alpha = "alpha"),
                     size = pt.size) +
          #scale_color_viridis(option="viridis") +
          scale_color_gradient2(low = 'green', mid = 'white', high = 'red') +
          scale_alpha(range = c(0.2, 1)) +
          ggtitle(features[i]) +
          guides(color = guide_legend(title = features[i])) +
          theme_void() +
          theme(legend.position="none",
                plot.title = element_text(hjust = 0.5))
      }
    }
  }
  return(plots)
}
