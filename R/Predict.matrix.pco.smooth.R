Predict.matrix.pco.smooth <- function(object, data){

  # somehow data gets passed as a list here :-\
  dat <- matrix(NA, ncol=length(data), nrow=length(data[[1]]))
  for(i in seq_along(data)){
    dat[,i] <- data[[i]]
  }

  mds.data <- as.matrix(insert.mds.generic(object$xt$mds.obj, dat,
                                           object$xt$realdata,
                                           object$xt$dist_fn))

  # make some variable names up
  colnames(mds.data) <- paste0("pco_", 1:ncol(mds.data))

  return(mds.data)
}
