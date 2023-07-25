#' DoSpark
#' @param object A \code{Seurat} object
#' @param n.core number of cores for parallel computing
#' @param genePercentCut The percentage of cells that are expressed for analysis
#' @param min_total_counts The minimum counts for each cell for filtering
#' @param method set automatically based on the sample size: if < 1000 spots,
#' \code{method} == "spark", otherwise "sparkx".
#'
#' @return a matrix
#' @importFrom SPARK CreateSPARKObject spark.vc spark.test sparkx
#' @importFrom Seurat Cells GetTissueCoordinates
#' @export
#'
DoSpark <- function(object,
                    n.core = 1,
                    genePercentCut = 0.05,
                    min_total_counts = 10,
                    method = ifelse(length(Cells(object)) < 1000, "spark",
                                    "sparkx")) {
  cell_n <- length(Cells(object))
  counts <- object@assays$Spatial@counts
  location <- GetTissueCoordinates(object)
  if(method=="spark"){
    spark <- CreateSPARKObject(counts = counts,
                               location = location,
                               percentage = genePercentCut,
                               min_total_counts = min_total_counts)
    spark@lib_size <- apply(spark@counts, 2, sum)
    spark <- spark.vc(spark,
                      covariates = NULL,
                      lib_size = spark@lib_size,
                      num_core = n.core,
                      fit.model="gaussian",
                      verbose = F)
    spark <- spark.test(spark,
                        check_positive = T,
                        verbose = F)
    result <- spark@res_mtest[,c("combined_pvalue","adjusted_pvalue")]
  } else if (method=="sparkx") {
    sparkX <- sparkx(counts,
                     location,
                     numCores = n.core,
                     option = "mixture")
    result <- sparkX$res_mtest
  } else {
    stop("method should be spark or sparkx")
  }
  return(result)
}
