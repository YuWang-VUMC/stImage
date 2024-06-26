% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SimilarityIntegration.R
\name{FindNeighborsByDistance}
\alias{FindNeighborsByDistance}
\title{FindNeighborsByDistance}
\usage{
FindNeighborsByDistance(
  dataObj = NULL,
  distMatrix,
  prune.SNN = 1/15,
  assay = NULL,
  graphs.name = c("nn", "snn"),
  nn.name = "neighbor",
  reduction.name = paste0(nn.name, "UMAP"),
  reduction.key = paste0(nn.name, "UMAP_"),
  k.param = 20
)
}
\arguments{
\item{dataObj}{A \code{Seurat} object.}

\item{distMatrix}{distance matrix in graphs slot}

\item{prune.SNN}{Sets the cutoff for acceptable Jaccard index when
computing the neighborhood overlap for the SNN construction. Any edges with
values less than or equal to this will be set to 0 and removed from the SNN
graph. Essentially sets the stringency of pruning (0 --- no pruning, 1 ---
prune everything).}

\item{assay}{assay name}

\item{graphs.name}{graph names: nn and snn}

\item{nn.name}{Name of knn output on which to run UMAP}

\item{reduction.name}{Name of projected UMAP to store in the query}

\item{reduction.key}{Value for the projected UMAP key}

\item{k.param}{number of k}
}
\value{
a seurat object
}
\description{
FindNeighborsByDistance
}
