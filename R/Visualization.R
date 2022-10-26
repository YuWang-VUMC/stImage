require(Palo)
require(sp)
require(sf)
require(RColorBrewer)
require(viridis)

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 95)[1:n]
}

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

AllSpatialDimPlot <- function(object,
                              cluster.frame = TRUE,
                              cluster.highlight = NULL,
                              col = NULL,
                              ...) {

  location <- GetTissueCoordinates(object)

  geneClusterName=grep("SCTPCA_cluster",colnames(object@meta.data),value=TRUE)
  renormgeneClusterName=grep("SCTNormalizedByImagePCA_cluster",colnames(object@meta.data),value=TRUE)
  geneSPCAClusterName=grep("SCTSpatialPCA_cluster",colnames(object@meta.data),value=TRUE)
  imageClusterName=grep("ImageFeaturePCA_cluster",colnames(object@meta.data),value=TRUE)
  RGBClusterName=grep("RGBPCA_cluster",colnames(object@meta.data),value=TRUE)
  wnnClusterName=grep("wnn_cluster",colnames(object@meta.data),value=TRUE)
  mciaClusterName=grep("mcia_snn_cluster",colnames(object@meta.data),value=TRUE)
  intnmfClusterName=grep("intnmf_snn_cluster",colnames(object@meta.data),value=TRUE)
  ticaClusterName=grep("tica_snn_cluster",colnames(object@meta.data),value=TRUE)

  groups <- c(geneClusterName, renormgeneClusterName, geneSPCAClusterName, imageClusterName, RGBClusterName, wnnClusterName, mciaClusterName, intnmfClusterName, ticaClusterName)
  groups <- groups[!is.na(groups)]

  plots <- list()
  for (i in 1:length(groups)) {
    cl <- setNames(object@meta.data[,groups[i]], rownames(object@meta.data[groups[i]]))
    location$cluster <- cl
    # arrange colors by Palo if not set by user
    pal <- gg_color_hue(length(unique(cl)))
    palopal <- Palo(location, cl, pal)
    if(is.null(col)){
      col <- palopal
    }

    if(cluster.frame){
      if(is.null(cluster.highlight)) {
        stop("Must provide the cluster ids when setting cluster.frame == True!")
      } else {
        # get cluster frame with sf
        p <- st_as_sf(data.frame(x=location$y, y=desc(location$x), A=1:nrow(location)), coords=1:2)
        v <- st_voronoi_point(p)
        hull <- st_convex_hull(st_union(p))

        vcut <- st_intersection(st_cast(v), hull)
        sf <- st_as_sf(vcut, celltype = cl)
        sf_0 <- sf[sf$celltype %in% as.character(cluster.highlight),]

        plots[[groups[i]]] <- ggplot(data = location) +
          geom_sf(data = st_union(sf_0), color = col[as.character(cluster.highlight)]) +
          geom_point(aes(x = y, y = desc(x), col = cluster), size = 1) +
          scale_color_manual(values = col) +
          ggtitle(groups[i]) +
          guides(color = guide_legend(title = groups[i])) +
          theme_void() +
          theme(legend.position="none",
                plot.title = element_text(hjust = 0.5))
      }
    } else {
      plots[[groups[i]]] <- ggplot(data = location) +
        #geom_sf(data = st_union(sf_0), color = col[as.character(cluster.highlight)]) +
        geom_point(aes(x = y, y = desc(x), col = cluster), size = 1) +
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

features <- c("TPT1", "IGHA1")

AllSpatialFeaturePlot <- function(object,
                                  normalizeMethod = c("SCT", "log"),
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
      exp <- object@assays$SCT@data
    } else if(normalizeMethod == "log"){
      exp <- object@assays$Spatial@data
    } else {
      stop("Please tell which slots of assays will be used to generate feature plots!")
    }
    exp <- as.matrix(exp)
    pseudo_v <- min(exp[exp > 0]) / 10^6
    scale_exp <- ScaleData(exp, features = row.names(exp))

    location <- GetTissueCoordinates(object)

    wnnClusterName=grep("wnn_cluster",colnames(object@meta.data),value=TRUE)
    mciaClusterName=grep("mcia_snn_cluster",colnames(object@meta.data),value=TRUE)
    intnmfClusterName=grep("intnmf_snn_cluster",colnames(object@meta.data),value=TRUE)
    ticaClusterName=grep("tica_snn_cluster",colnames(object@meta.data),value=TRUE)

    if(integreation.method == "wnn") {
      cl <- wnnClusterName
    } else if(integreation.method == "mcia") {
      cl <- mciaClusterName
    } else if(integreation.method == "intnmf") {
      cl <- intnmfClusterName
    } else if(integreation.method == "tica") {
      cl <- ticaClusterName
    }

    plots <- list()
    for (i in 1:length(features)) {

      scale_features <- t(scale_exp[features[i], , drop = F])
      scale_features <- cbind(location, scale_features)
      # arrange colors by Palo if not set by user
      pal <- gg_color_hue(length(unique(cl)))
      palopal <- Palo(location, cl, pal)
      if(is.null(col)){
        col <- palopal
      }

      if(cluster.frame){
        if(is.null(cluster.highlight)) {
          stop("Must provide the cluster ids when setting cluster.frame == True!")
        } else {
          # get cluster frame with sf
          p <- st_as_sf(data.frame(x=location$y, y=desc(location$x), A=1:nrow(location)), coords=1:2)
          v <- st_voronoi_point(p)
          hull <- st_convex_hull(st_union(p))

          vcut <- st_intersection(st_cast(v), hull)
          sf <- st_as_sf(vcut, celltype = cl)
          sf_0 <- sf[sf$celltype %in% as.character(cluster.highlight),]

          plots[[features[i]]] <- ggplot(scale_features) +
            geom_sf(data = st_union(sf_0), color = col[as.character(cluster.highlight)]) +
            geom_point(aes_string(x = scale_features$y, y = desc(scale_features$x), col = features[i]), size = 1) +
            scale_color_viridis(option="viridis") +
            ggtitle(features[i]) +
            #scale_color_brewer(palette = "PiYG") +
            #scale_color_gradient2(low = 'blue', mid = 'white', high = 'red') +
            guides(color = guide_legend(title = features[i])) +
            theme_void() +
            theme(legend.position="none",
                  plot.title = element_text(hjust = 0.5))
        }
      } else {
        plots[[features[i]]] <- ggplot(scale_features) +
          #geom_sf(data = st_union(sf_0), color = col[as.character(cluster.highlight)]) +
          geom_point(aes_string(x = scale_features$y, y = desc(scale_features$x), col = features[i]), size = 1) +
          scale_color_viridis(option="viridis") +
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
