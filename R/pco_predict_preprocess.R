#' Make predictions using pco basis terms
#'
#' This function performs the necessary preprocessing for making predictions with \code{\link[mgcv]{gam}} models that include \code{\link{pco}} basis terms. The function \code{pco_predict_preprocess} builds a \code{data.frame} (or augments an existing one) to be used with the usual \code{predict} function.
#'
#' Models with \code{\link{pco}} basis terms are fitted by inputting distances among the observations and then regressing (with a ridge penalty) on leading principal coordinates arising from these distances. To perform prediction, we must input the distances from the new data points to the original points, and then "insert" the former into the principal coordinate space by the interpolation method of Gower (1968).  
#'
#' @param model a fitted \code{\link[mgcv]{gam}} model with at least one term of class "\code{pco.smooth}"
#' @param newdata data frame including the new values for any non-\code{\link{pco}} terms in the original fit. If there were none, this can be left as \code{NULL}.
#' @param dist_list a list of \code{n*\times n} matrices, one per \code{\link{pco}} term in the model, giving the distances from the \code{n*} prediction points to the \code{n} design points (original observations). List entry names should correspond to the names of the terms in the model (e.g., if the model includes a \code{s(x)} term, \code{dist_list} must include an element named "\code{x}").
#'
#' @return a \code{\link{data.frame}} with the coordinates for the new data inserted into principal coordinate space, in addition to the supplied \code{newdata} if this was non-\code{NULL}. This can be used as the \code{newdata} argument in a call to \code{\link[mgcv]{predict.gam}}.
#' @author David L Miller
#' @export
#' @references
#' Gower, J. C. (1968). Adding a point to vector diagrams in multivariate analysis. Biometrika, 55(3), 582-585. \url{http://doi.org/10.2307/2334268}
#'
#' Miller, D. L. (2012). On smooth models for complex domains and distances. PhD dissertation, Department of Mathematical Sciences, University of Bath. Available at \url{http://opus.bath.ac.uk/31800/}
#' @seealso \code{\link{smooth.construct.pco.smooth.spec}}
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

