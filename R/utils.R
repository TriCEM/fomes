#' @title Update  Network Connections
#' @inheritParams sim_Gillespie_SIR
#' @param adjmat matrix; adjacency matrix from network of connections
#' @noMd
#' @export
updateNetworkConnections <- function(adjmat = NULL, rho = 0.5) {
  #............................................................
  # rewiring probability
  #   NB edge density does not change per node; just arcs between nodes
  #...........................................................
  shuffle_connections <- function(x, rho) {
    nupd <- floor(rho * sum(x))
    old <- which(x == 1)
    new <- which(x == 0)
    # old losses conn but new gains conn
    x[ sample(old, size = nupd) ] <- 0
    x[ sample(new, size = nupd) ] <- 1
    return(x)
  }

}


#' @title Generate Random Network for Connections
#' @inheritParams sim_Gillespie_SIR
#' @noMd
#' @export

genRandomNetworkConnections <- function(N, rho, initNC) {
  net <- igraph::watts.strogatz.game(dim = 1,
                                     size = N,
                                     nei = initNC,
                                     p = rho)

  adjmat <- igraph::as_adjacency_matrix(net,
                                        sparse = FALSE)

  return(adjmat)
}



#' @title Generate Random Beta matrix
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
