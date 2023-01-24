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





#' @title Update Adjacency Matrix Connections
#' @inheritParams sim_Gillespie_SIR
#' @details Will sample a single pair of connections
#' @noMd
#' @export
rewireNEnodes <- function(adjmat = NULL, N) {
  #............................................................
  # rewiring structure by randomly sampling two nodes with connection (arc)
  #   NB edge density does not change per node; just arcs between nodes
  #...........................................................
  ctch <- TRUE # init catch for proper node selection
  while (ctch) {
    currconn <- which(adjmat == 1)
    # sample two current connections
    currconn_s <- sample(currconn, size = 2, replace = F)
    currconn_i <- ceiling(currconn_s/N) # get i
    currconn_j <- currconn_s %% N # get j
    ab <- c(currconn_i[1], currconn_j[1])
    cd <- c(currconn_i[2], currconn_j[2])

    # catch for diagonal or same connection on opposite sides of triangle (lower vs upper)
    ctch <- length(unique(ab)) == 1 | length(unique(cd)) == 1 | paste(sort(ab), collapse = "") == paste(sort(cd), collapse = "")
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
  return(adjmat)
}






#' @title Initialize Adjacency Matrix Connections
#' @inheritParams sim_Gillespie_SIR
#' @noMd
#' @export
genInitialConnections <- function(initNC, N) {

  # initial contacts, assume a binomial prob dist
  #initedgedens <- ceiling( N * (1 - exp(-rho)) ) # round to nearest whole number for edge (ceiling so always at least 1)
  #initedgedens <- initNC

  #......................
  # greedy approach to make initial adjacency matrix
  #......................
  greedy_contactmat_generator <- function(conn) {
    for(i in 1:N) { # for each node
      while (sum(conn[i,]) < initNC) { # until we have init edge density
        newconns <- which(colSums(conn) != initNC)
        newconns <- newconns[newconns != i]
        if (length(newconns) == 0) { # all connections saturated
          break
        } else {
          j <- sample(newconns, 1) # sample a connection
          if (sum(conn[,j]) < initNC) { # look ahead to make sure new node is not saturated (since the check is outside while)
            conn[i,j] <- 1
            conn[j,i] <- 1
          }
        }
      }
    } # end for loop
    return(conn)
  }
  #......................
  # run function
  #......................
  contactmat <- matrix(0, N, N)
  while(all(rowSums(contactmat) != initNC)) { # recursively doing this for transitivity issue leading to early saturation w/out all nodes having equal edge density
    contactmat <- greedy_contactmat_generator(conn = contactmat)
  }

  # isSymmetric(conn) # must be true
  # out
  return(contactmat)
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
