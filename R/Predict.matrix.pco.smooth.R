Predict.matrix.pco.smooth <- function(object, data){

  # somehow data gets passed as a list here :-\
  dat <- matrix(data[[object$term]], ncol=object$dim)

  return(dat)
}
