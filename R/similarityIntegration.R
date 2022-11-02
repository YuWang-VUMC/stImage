
distToNeighborObj=function(distMatrix,k.param=20) {

  #distMatrix=as.matrix(bayesSpaceResultDist)
  n.cells=nrow(distMatrix)
  knn.mat <- matrix(data = 0, ncol = k.param, nrow = n.cells)
  knd.mat <- knn.mat
  for (i in 1:n.cells) {
    knn.mat[i, ] <- order(distMatrix[i, ])[1:k.param]
    knd.mat[i, ] <- distMatrix[i, knn.mat[i, ]]
  }
  # nn.ranked <- knn.mat[, 1:k.param]
  # nn.idx=t(apply(nnMatrix,1,function(x) which(x>0)))
  # nn.dist=matrix(NA,ncol=ncol(nn.idx),nrow=nrow(nn.idx))
  # for (i in 1:nrow(nn.idx)) {
  #   nn.dist[i,]=distMatrix[i,nn.idx[i,]]
  # }

  objNeighbor <- new(
    Class = 'Neighbor',
    nn.idx = knn.mat,
    nn.dist = knd.mat,
    alg.info = list(),
    cell.names = rownames(x = distMatrix)
  )
  return(objNeighbor)
}

# NeighborDist=function(distMatrix,objNeighbor) {
#   n.cells=nrow(objNeighbor@nn.idx)
#   distSummary=structure(rep(NA,n.cells),names=objNeighbor@cell.names)
#   distToNearest=structure(rep(NA,n.cells),names=objNeighbor@cell.names)
#   for (i in 1:n.cells) {
#     #distSummary[i]=mean(distMatrix[i,objNeighbor@nn.idx[i,]])
#     distSummary[i]=mean(distMatrix[i,objNeighbor@nn.idx[i,]][-1]) #need to remove self
#     distToNearest[i]=(distMatrix[i,objNeighbor@nn.idx[i,]][2])  #1 is self
#   }
#   return(data.frame(distSummary,distToNearest))
# }

NeighborDist=function(distMatrix,neighborList,removeFirst=TRUE) {
  if (class(neighborList)=="Neighbor") { #make it a list
    neighborList=lapply(data.frame(t(neighborList@nn.idx)),function(x) x)
  }
  n.cells=length(neighborList)
  k=length(neighborList[[1]])
  distSummary=structure(rep(NA,n.cells),names=row.names(distMatrix))
  distToNearest=structure(rep(NA,n.cells),names=row.names(distMatrix))
  if (removeFirst) {
    allNeighborInd=2:k
    nearestInd=2
  } else {
    allNeighborInd=1:k
    nearestInd=1
  }
  for (i in 1:n.cells) {
    #distSummary[i]=mean(distMatrix[i,objNeighbor@nn.idx[i,]])
    distSummary[i]=mean(distMatrix[i,neighborList[[i]]][allNeighborInd]) #need to remove self
    distToNearest[i]=(distMatrix[i,neighborList[[i]]][nearestInd])  #1 is self
  }
  return(data.frame(distSummary,distToNearest))
}




