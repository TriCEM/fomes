#' @title Tidy out the sim_Gillespie_SIR Model
#' @inheritParams sim_Gillespie_SIR
#' @noMd
#' @export
#TODO overload this once determine class type - see github issue #5
tidy_sim_Gillespie_SIR <- function(simout) {
  out <- data.frame(
    time = simout$time,
    numSusc = rowSums(simout$S_traj),
    numInfxn = rowSums(simout$I_traj),
    numRecov = rowSums(simout$R_traj)
  )
  return(out)
}



#' @title Update Network Connections
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

  # in order to preserve symmetry, will work with just the upper (and lower triangle)
  for(i in 1:(niters-1)) {
    for (j in (i+1):niters) {
      r <- runif(1) # generate random prob

      if (r < prho & adjmat[i,j] == 0) {
        # establish this new connection (do this first in case no connections exist, which we will undo in next few lines)
        adjmat_new[i,j] <- 1
        adjmat_new[j,i] <- 1
        # dissolve an old connection
        dissolve <- sample(which(adjmat_new[i,] == 1), 1)
        adjmat_new[dissolve,j] <- 0
        adjmat_new[j,dissolve] <- 0

      } else if (r < prho & adjmat[i,j] == 1) {
        # establish this new connection
        adjmat_new[i,j] <- 0
        adjmat_new[j,i] <- 0
        # dissolve an old connection
        dissolve <- sample(which(adjmat_new[i,] == 0), 1)
        adjmat_new[dissolve,j] <- 1
        adjmat_new[j,dissolve] <- 1
      }
    }
  } # end ij for upper

  # out
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
