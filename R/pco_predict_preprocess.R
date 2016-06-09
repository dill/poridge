#' Make predictions using the pco basis
#'
#' Making predictions using the \code{pco} basis requires some pre-processing of the data first. The function \code{pco_predict_preprocess} builds a \code{data.frame} (or augments an existing one) to be used with the usual \code{predict} function.
#'
#' The function takes uses the distances from the data in the model to the prediction points and uses Gower's interpolation to "insert" the prediction points into the existing multidimensional scaling projection. These new coordinates are returned.
#'
#' @param model a fitted model with at least one term of class \code{pco.smooth}
#' @param newdata prediction data
#' @param dist_list a list of distance matrices, each of which has distances from observations to prediction points. Each matrix should have number of rows identical to the number of observations in the model and as many columns as there are prediction points. List entries need the name of the term they refer to as their element name (e.g., if the term is \code{s(x)}, then the corresponding matrix element will be named "\code{x}").
#'
#' @return a \code{data.frame}, either the \code{newdata} augmented with new columns or a \code{data.frame} with just those columns.
#' @author David L Miller
#' @export
#' @references
#' Gower, J. C. (1968). Adding a point to vector diagrams in multivariate analysis. Biometrika, 55(3), 582. http://doi.org/10.2307/2334268
#'
#' Miller, D. L. (2012). On Smooth Models for Complex Domains and Distances. (S. N. Wood, Ed.). Retrieved from http://opus.bath.ac.uk/31800/
pco_predict_preprocess <- function(model, newdata=NULL, dist_list){

  # populate newdata
  destroy_col <- FALSE
  if(is.null(newdata)){
    newdata <- data.frame(newdata_dummy_data = rep(NA,ncol(dist_list[[1]])))
    # reminder to destroy that extra column
    destroy_col <- TRUE
  }

  # which terms in the model are pco terms?
  which_pco <- which(unlist(lapply(model$smooth,
                             function(x) any(class(x)=="pco.smooth"))))

  # die if there are no pco terms
  if(length(which_pco)==0){
    stop("There are no pco smooths in the model")
  }

  # loop over smooths of the right type
  for(i in which_pco){

    # this term
    term_lab <- model$smooth[[i]]$term

    # get this distance matrix
    distmat <- dist_list[[term_lab]]

    # goofy - build the data to put into the equation later on
    mds.obj <- model$smooth[[i]]$xt$mds.obj
    X <- mds.obj$x
    S <- -1/2*X
    d <- distmat
    non.diag <- row(d) != col(d)
    d[non.diag] <- (d[non.diag] + mds.obj$ac)
    d <- -(d^2-diag(S))

    # lambda^-1, 1/eigenvalues in a diagonal matrix
    ev <- mds.obj$eig[1:ncol(mds.obj$points)]
    if(length(ev)>1){
      lambda.inverse <- diag(1/ev)
    }else{
      lambda.inverse <- as.matrix(1/ev)
    }

    # conversion from maths into code of equation 10 in Gower 1968
    mds.points <- t(1/2*(lambda.inverse %*% t(mds.obj$points) %*% d))

    # make the prediction frame for this term
    preddat <- as.data.frame(mds.points)

    # give the columns names
    names(preddat) <- paste0("xdummy", "_" ,1:ncol(mds.points))

    # splodge it into the prediction frame
    newdata[[term_lab]] <- mds.points
  }

  # get rid of the column we needed to setup the data
  if(destroy_col){
    newdata[["newdata_dummy_data"]] <- NULL
  }

  return(newdata)
}

