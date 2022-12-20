#' @title Generate Random Beta matrix for Effective ij Contact Rates
#' @param N integer; population size
#' @noMd
#' @export

genRandomBetaMat <- function(N) {
  ret <- matrix(runif(N^2), ncol = N, nrow = N)
  diag(ret) <- 0
  return(ret)
}

#' @title Bind SIR trajectories
#' @param x list of vectors matrix
#' @noMd
#' @export
tidy_traj_out <- function(x) {
  if (length(x) == 1) {
    ret <- matrix(data = unlist(x), ncol = 1)
    rownames(ret) <- NULL
    colnames(ret) <- "1"
  } else {
    ret <- t(do.call("cbind", (lapply(x, as.data.frame))))
    rownames(ret) <- NULL
    colnames(ret) <- as.character(1:ncol(ret))
  }
  return(as.data.frame(ret))
}
