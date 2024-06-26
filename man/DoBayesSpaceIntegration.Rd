% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DoBayesSpaceIntegration.R
\name{DoBayesSpaceIntegration}
\alias{DoBayesSpaceIntegration}
\title{DoBayesSpaceIntegration}
\usage{
DoBayesSpaceIntegration(
  dataObj,
  clusterNum = 10,
  BayesSpaceClusterNum = clusterNum,
  assay = c("SCT", "ImageFeature")
)
}
\arguments{
\item{dataObj}{A \code{Seurat} object.}

\item{clusterNum}{number of clusters}

\item{BayesSpaceClusterNum}{number of clusters by BayesSpace method}

\item{assay}{the two assays used for integration}
}
\value{
a seurat object
}
\description{
DoBayesSpaceIntegration
}
\examples{
\dontrun{
  object <- DoBayesSpaceIntegration(
    object,
    clusterNum = clusterNum)
}

}
