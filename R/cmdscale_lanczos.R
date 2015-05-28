#' Faster multi-dimensional scaling
#'
#' This is a stand-in replacement for \code{\link{cmdscale}} that uses the Lanczos procedure (\code{slanczos}) in \code{\link{mgcv}} instead of \code{eigen}.
#'
#' Much of the code is based on \code{cmdscale}.
#' @param d a distance structure such as that returned by 'dist' or a full symmetric matrix containing the dissimilarities.
#' @param k the maximum dimension of the space which the data are to be represented in; must be in \code{\{1, 2, ..., n-1\}}.
#' @param eig indicates whether eigenvalues should be returned.
#' @param add logical indicating if an additive constant c* should be computed, and added to the non-diagonal dissimilarities such that the modified dissimilarities are Euclidean.
#' @param x.ret indicates whether the doubly centred symmetric distance matrix should be returned.
#' @return as \code{\link{cmdscale}}
#'
#' @author David L Miller
# @importFrom mgcv slanczos
cmdscale_lanczos <- function(d, k = 2, eig = FALSE, add = FALSE, x.ret = FALSE){

  if (anyNA(d))
    stop("NA values not allowed in 'd'")
  if (is.null(n <- attr(d, "Size"))) {
    if(add) d <- as.matrix(d)
    x <- as.matrix(d^2)
    storage.mode(x) <- "double"
    if ((n <- nrow(x)) != ncol(x))
      stop("distances must be result of 'dist' or a square matrix")
    rn <- rownames(x)
  } else {
    rn <- attr(d, "Labels")
    x <- matrix(0, n, n) # must be double
    if (add) d0 <- x
    x[row(x) > col(x)] <- d^2
    x <- x + t(x)
    if (add) {
      d0[row(x) > col(x)] <- d
      d <- d0 + t(d0)
    }
  }
  n <- as.integer(n)

  ## we need to handle nxn internally in dblcen
  if(is.na(n) || n > 46340) stop("invalid value of 'n'")
  if((k <- as.integer(k)) > n - 1 || k < 1)
    stop("'k' must be in {1, 2, .. n - 1}")

  ## NB: this alters argument x, which is OK as it is re-assigned.
  #x <- .Call(stats::C_DoubleCentre, x)
  x <- scale(t(scale(t(x), scale=FALSE)),scale=FALSE)

  if(add) { ## solve the additive constant problem
    ## it is c* = largest eigenvalue of 2 x 2 (n x n) block matrix Z:
    i2 <- n + (i <- 1L:n)
    Z <- matrix(0, 2L*n, 2L*n)
    Z[cbind(i2,i)] <- -1
    Z[ i, i2] <- -x
#    Z[i2, i2] <- .Call(stats::C_DoubleCentre, 2*d)
    Z[i2, i2] <- scale(t(scale(t(2*d), scale=FALSE)),scale=FALSE)

    ###### this is where Dave modified things
    add.c <- max(slanczos(Z, k=1, kl=1)$values)
    #e <- eigen(Z, symmetric = FALSE, only.values = TRUE)$values
    #add.c <- max(Re(e))
    ## and construct a new x[,] matrix:
    x <- matrix(double(n*n), n, n)
    non.diag <- row(d) != col(d)
    x[non.diag] <- (d[non.diag] + add.c)^2
    #x <- .Call(stats::C_DoubleCentre, x)
    x <- scale(t(scale(t(x), scale=FALSE)),scale=FALSE)
  }

  ###### this is where Dave modified things
  e <- slanczos(-x/2, k=k)
  ev <- e$values#[seq_len(k)]
  evec <- e$vectors#[, seq_len(k), drop = FALSE]
  k1 <- sum(ev > 0)

  if(k1 < k) {
    warning(gettextf("only %d of the first %d eigenvalues are > 0", k1, k),
            domain = NA)
    evec <- evec[, ev > 0, drop = FALSE]
    ev <- ev[ev > 0]
  }

  points <- evec * rep(sqrt(ev), each=n)
  dimnames(points) <- list(rn, NULL)
  if (eig || x.ret || add) {
    evalus <- e$values # Cox & Cox have sum up to n-1, though
    list(points = points, eig = if(eig) evalus, x = if(x.ret) x,
         ac = if(add) add.c else 0,
         GOF = sum(ev)/c(sum(abs(evalus)), sum(pmax(evalus, 0))) )
  } else points
}
