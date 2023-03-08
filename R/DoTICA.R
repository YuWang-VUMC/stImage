#### tensorial ICA: tWFOBI and tWJADE
# inherited from:
# https://github.com/cantinilab/momix-notebook/blob/master/scripts/tICA.R
#' tICA
#'
#' @param data this is the data object, which contains the multi-way data in
#' two different formats. The "A" entry of data (data$A) gives the array or
#' data-tensor format, whereas the "L" entry of data (data$L) gives the data in
#' list format. In the former case, and in our specific applications, the first
#' mode defines the tissue, cell or data-type, the second mode defines the
#' samples and the third mode the features (e.g. CpGs or genes). In the latter
#' case, each list entry corresponds to the cell/tissue or data-type and
#' consists of the data-matrix which rows representing features and columns
#' representing samples.
#' @param dim this is a vector which contains the number of significant
#' components of each data matrix to search for, and is typically obtained by
#' applying RMT to each separate data/tissue-type matrix (i.e. to the
#' individual entries of data$L above).
#' @param method either \code{FOBI} or \code{JADE}
#'
#' @return
#' @importFrom tensorBSS tensorCentering tPCA tensorTransform tFOBI tJADE
#' @export
#'
#' @examples
#'
DoTICA <- function(data,
                   dim,
                   method = c("FOBI", "JADE")){
  nt <- dim(data)[1]
  ng <- dim(data)[3]
  dim.v <- c(nt, max(dim))
  cdata.a <- tensorCentering(data)
  tpca.o <- tPCA(cdata.a, d=dim.v)
  pdata.a <- tensorTransform(cdata.a, t(tpca.o$U[[2]]), 2) ## whiten the data


  if(method == "FOBI"){
    tica.o <- tFOBI(pdata.a)
  } else {
    tica.o <- tJADE(pdata.a)
  }

  projS.lm <- list()
  for(t in 1:nt){
    projS.lm[[t]] <- tica.o$S[t,,]
  }
  signals <- tpca.o$U[[2]] %*% tica.o$W[[2]] ## back to original space
  return(list(projS = projS.lm,
              S = tica.o$S,
              U = tica.o$W,
              nt = nt,
              ng = ng,
              signals = signals))

}
