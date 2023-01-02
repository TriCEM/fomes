#' @title Update  Network Connections
#' @inheritParams sim_Gillespie_SIR
#' @param adjmat matrix; adjacency matrix from network of connections
#' @noMd
#' @export
updateNetworkConnections <- function(adjmat = NULL, rho = 0.5) {
  prho <- 1 - exp(-rho) # rate to prob
  #............................................................
  # rewiring probability
  #   NB edge density does not change per node; just arcs between nodes
  #...........................................................
  niters <- nrow(adjmat)
  adjmat_new <- adjmat

  # in order to preserve symmetry, will work with just the upper or lower triangle
  # toggle upper or lower for balance of edge densities
  #......................
  # upper triangle
  #......................
  if (rbinom(1, 1, 0.5)) {

    for(i in 1:(niters-1)) {
      for (j in (i+1):niters) {
        r <- runif(1) # generate random prob

        # establishing new connection
        if (r < prho & adjmat[i,j] == 0) {
          adjmat_new[i,j] <- 1
        }

        # dissolving old connection
        if (r < (1 - prho) & adjmat[i,j] == 1) {
          adjmat_new[i,j] <- 0
        }

      }
    } # end ij for upper
    # now make symmetric
    adjmat_new[lower.tri(adjmat_new)] <- t(adjmat_new)[lower.tri(adjmat_new)]

    #......................
    # lower triangle
    #......................
  } else {
    for(i in niters:2) {
      for (j in (i-1):1) {
        r <- runif(1) # generate random prob

        # establishing new connection
        if (r < prho & adjmat[i,j] == 0) {
          adjmat_new[i,j] <- 1
        }

        # dissolving old connection
        if (r < (1 - prho) & adjmat[i,j] == 1) {
          adjmat_new[i,j] <- 0
        }

      }
    } # end ij for lower triangle
    # now make symmetric
    adjmat_new[upper.tri(adjmat_new)] <- t(adjmat_new)[upper.tri(adjmat_new)]
  }

  #......................
  # out
  #......................
  return(adjmat_new)
}






#' @title Generate Random Network for Connections
#' @inheritParams sim_Gillespie_SIR
#' @noMd
#' @export

genRandomNetworkConnections <- function(N, rho, initNC) {

  prho <- 1 - exp(-rho) # rate to prob
  net <- igraph::watts.strogatz.game(dim = 1,
                                     size = N,
                                     nei = initNC,
                                     p = prho)

  adjmat <- igraph::as_adjacency_matrix(net,
                                        sparse = FALSE) # sparse here is type of
  # isSymmetric(adjmat) - is symmetric by default

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
