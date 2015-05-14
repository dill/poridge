insert.mds.generic <- function(mds.obj, new.points, old.points, dist_fn){
  
  ## NB mds.obj$points is double centred
  
  
  # mudge old and new data into one frame, keep track of the indices
  big.points <- rbind(old.points, new.points)
  ind <- 1:nrow(old.points)
  ind2 <- (nrow(old.points)+1):nrow(big.points)
  
  # calculate the distance between the data used to fit the model
  # and the data we want to predict for
  new.dist <- as.matrix(dist_fn(big.points))
  
  # get the new distances that we want
  new.dist <- new.dist[ind,][,ind2]
  
  # calculate S=XX'
  X <- mds.obj$x
  S <- -1/2*X
  
  # calculate d matrix
  #d <- -(new.dist^2-diag(S))
  # take into account additive constant
  d<- new.dist
  non.diag <- row(d) != col(d)
  d[non.diag] <- (d[non.diag] + mds.obj$ac)
  d<- -(d^2-diag(S))
  
  # lambda^-1, 1/eigenvalues in a diagonal matrix
  ev <- mds.obj$eig[1:ncol(mds.obj$points)]
  lambda.inverse <- diag(1/ev)
  
  # conversion from maths into code of equation 10 in Gower 1968
  mds.points <- t(1/2*(lambda.inverse %*% t(mds.obj$points) %*% d))
  
  return(mds.points)
}