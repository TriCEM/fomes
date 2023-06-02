
#' @title Initialize Adjacency Matrix Connections
#' @inheritParams sim_Gillespie_nSIR
#' @param sparseMatrix boolean; Whether the returned adjacency matrix should be formatted
#' as a sparseMatrix
#' @noMd
#' @export
genInitialConnections <- function(initNC, N, sparseMatrix = TRUE) {

  # initial contacts, assume a binomial prob dist
  #initedgedens <- ceiling( N * (1 - exp(-rho)) ) # round to nearest whole number for edge (ceiling so always at least 1)
  #initedgedens <- initNC

  #............................................................
  # using igraph as workhorse here
  #...........................................................
  net <- igraph::degree.sequence.game(rep(initNC, N), method = "vl")
  contactmat <- igraph::as_adjacency_matrix(net,
                                            type = "both",
                                            names = FALSE,
                                            sparse = sparseMatrix)
  return(contactmat)
}




#' @title Update Adjacency Matrix Connections
#' @inheritParams sim_Gillespie_nSIR
#' @param adjmat sparse matrix; adjacency contact matrix that will be rewired using the neighbor exchange model
#' @param sparseMatrix boolean; Whether the returned adjacency matrix should be formatted
#' as a sparseMatrix
#' @details Will sample a single pair of connections
#' @noMd
#' @export
rewireNEnodes <- function(adjmat = NULL, N, sparseMatrix = TRUE) {

  # temporarily convert sparseMatrix to base matrix for which
  if (any(c("dgeMatrix", "dgCMatrix", "dsCMatrix") %in% class(adjmat))) {
    adjmat <- as.matrix(adjmat)
  }

  #............................................................
  # rewiring structure by randomly sampling two nodes with connection (arc)
  #   NB edge density does not change per node; just arcs between nodes
  #...........................................................
  ctch <- TRUE # init catch for proper node selection
  while (ctch) {
    currconn <- which(adjmat == 1)
    # sample two current connections
    currconn_s <- sample(currconn, size = 2, replace = F)
    currconn_i <- currconn_s %% N # get i
    currconn_i[currconn_i == 0] <- N # no remainder means last row
    currconn_j <- ceiling(currconn_s/N) # get j
    ab <- c(currconn_i[1], currconn_j[1])
    cd <- c(currconn_i[2], currconn_j[2])
    #......................
    # catches
    #......................
    # can't be a diagonal (granted diagonals are coded as 0s so this is to protect any future changes)
    diagctch <- length(unique(ab)) == 1 | length(unique(cd)) == 1
    # can't be on opposite sides of the triangle
    trigctch <- paste(sort(ab), collapse = "") == paste(sort(cd), collapse = "")
    # can't share nodes between pairs
    ndctch <- length(unique(c(ab, cd))) != 4
    ctch <- any(c(diagctch, trigctch, ndctch))
  }

  # erase prior connections
  adjmat[ab[1], ab[2]] <- 0
  adjmat[ab[2], ab[1]] <- 0
  adjmat[cd[1], cd[2]] <- 0
  adjmat[cd[2], cd[1]] <- 0

  # identify new connections
  ctch <- TRUE # init catch for ensuring non-redundant connections
  while (ctch) {
    newconn_index <- sample(1:4, size = 4, replace = F)
    abnew <- c(ab,cd)[ c(newconn_index[1:2]) ] # new arc 1
    cdnew <- c(ab,cd)[ c(newconn_index[3:4]) ] # new arc 2

    # catch if there is already a connection there avoiding redundancy and loss of consistent edge density
    # NB this will allow for original connections to reform
    ctch <- adjmat[abnew[1], abnew[2]] == 1 | adjmat[cdnew[1], cdnew[2]] == 1

  }

  # make new connections
  adjmat[abnew[1], abnew[2]] <- 1
  adjmat[abnew[2], abnew[1]] <- 1
  adjmat[cdnew[1], cdnew[2]] <- 1
  adjmat[cdnew[2], cdnew[1]] <- 1

  # out
  if (sparseMatrix) {
    adjmat <- as(adjmat, "sparseMatrix")
  }
  return(adjmat)
}


