#' SpatialPCAWorkflow
#' @inheritParams SpatialPCA
#' @param SpatialPCAObj A \code{SpatialPCA} object
#' @param SpatialPCnum number of PCs when running SpatialPCA
#' @param kerneltype The type of kernel to be used, either "gaussian" for
#' gaussian kernel, or "cauchy" for cauchy kernel, or "quadratic" for rational
#' quadratic kernel.
#' @param bandwidthtype The type of bandwidth to be used in Gaussian kernel,
#' "SJ" for Sheather & Jones (1991) method (usually used in small sample size
#' datasets), "Silverman" for Silverman's ‘rule of thumb’ method (1986)(usually
#' used in large sample size datasets).
#' @param fast Select "TRUE" to accrelerate the algorithm by performing
#' low-rank approximation on the kernel matrix, otherwise "FALSE" for
#' calculation without low-rank approximation on the kernel matrix.
#'
#' @return
#' @importFrom SpatialPCA SpatialPCA_buildKernel
#' @importFrom SpatialPCA SpatialPCA_EstimateLoading
#' @importFrom SpatialPCA SpatialPCA_SpatialPCs
#' @export
#'
#' @examples
#' \dontrun{SPCAobj <- SpatialPCAWorkflow(SPCAobj, SpatialPCnum = pcaDim)}
#'
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
