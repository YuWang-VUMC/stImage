#' Title
#'
#' @param SpatialPCAObj
#' @param SpatialPCnum
#' @param kerneltype
#' @param bandwidthtype
#' @param fast
#'
#' @return
#' @export
#'
#' @examples
SpatialPCAWorkflow <- function(SpatialPCAObj,
                               SpatialPCnum = 30,
                               kerneltype = "gaussian",
                               bandwidthtype = "SJ",
                               fast = FALSE) {
  SpatialPCAObj <- SpatialPCA_buildKernel(SpatialPCAObj,
                                          kerneltype = kerneltype,
                                          bandwidthtype = bandwidthtype)
  SpatialPCAObj <- SpatialPCA_EstimateLoading(SpatialPCAObj,
                                              fast = fast,
                                              SpatialPCnum = SpatialPCnum)
  SpatialPCAObj <- SpatialPCA_SpatialPCs(SpatialPCAObj,
                                         fast = fast)

  return(SpatialPCAObj)
}