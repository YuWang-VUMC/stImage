#### tensorial ICA: tWFOBI and tWJADE
#' tICA
#'
#' @param data
#' @param dim
#' @param method
#'
#' @return
#' @importFrom tensorBSS tensorCentering tPCA
#' tensorTransform tFOBI tJADE
#' @export
#'
#' @examples
DoTICA <- function(data,dim,method=c("FOBI","JADE")){
  nt <- dim(data)[1]
  ng <- dim(data)[3]
  dim.v <- c(nt,max(dim))
  cdata.a <- tensorCentering(data)
  tpca.o <- tPCA(cdata.a,d=dim.v)
  pdata.a <- tensorTransform(cdata.a,t(tpca.o$U[[2]]),2) ## whiten the data


  if(method=="FOBI"){
    tica.o <- tFOBI(pdata.a)
  } else {
    tica.o <- tJADE(pdata.a)
  }

  projS.lm <- list()
  for(t in 1:nt){
    projS.lm[[t]] <- tica.o$S[t,,]
  }
  signals <- tpca.o$U[[2]] %*% tica.o$W[[2]] ## back to original space
  return(list(projS=projS.lm,S=tica.o$S,U=tica.o$W,nt=nt,ng=ng,signals=signals))

}
