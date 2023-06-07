#' @title SIR NE Gillespie-Algorithm Simulation Model
#' @param Iseed integer; Number of initial infections in the population
#' @param N integer; population size
#' @param beta numeric vector, matrix, or sparseMatrix; Probability of transmission given contact for S_i I_j elements
#' @param dur_I numeric; Duration of infection
#' @param term_time numeric; time termination condition
#' @param init_contact_mat matrix or sparseMatrix; User specified initial contact, or adjacency matrix.
#' The contact matrix must be a symmetric square matrix, with a diagonal of 0 (for self),
#' and have the same edge density, number of connections, for each node, or individual.
#' Default value of NULL results in internal function generating initial contact matrix.
#' @param initNC NC integer; initial connections; NB only used if an initial matrix is not provided.
#' @param rho numeric; network rewiring rate
#' @param return_contact_matrices boolean; Option to return contact matrices throughout dynamic iterations.
#' Warning, this can be memory intensive.
#' @description Gillespie Simulation
#' @details Assuming:
#' \itemize{
#'     closed population
#' }
#' @details If the user does not input an initial contact matrix, the package will generate an initial contat matrix
#' with the same edge density for each node using the \code{igraph::degree.sequence.game} function
#' @returns The return:
#' \itemize{
#'     conn_traj boolean whether or not event included a rewiring
#' }
#' @export

sim_Gillespie_nSIR <- function(Iseed = 1, N = 10,
                              beta = rep(1, 10),
                              dur_I = 1,
                              init_contact_mat = NULL,
                              rho = 0.05,
                              initNC = NA,
                              term_time = 500,
                              return_contact_matrices = FALSE) {
  # runtime
  starttime <- Sys.time()
  #............................................................
  # assertions
  #...........................................................
  goodegg::assert_single_int(Iseed)
  goodegg::assert_single_int(N)
  goodegg::assert_single_numeric(term_time)
  goodegg::assert_logical(return_contact_matrices)

  #.........................
  # initial contact adjacency matrix
  #.........................
  goodegg::assert_single_numeric(rho)
  goodegg::assert_gr(rho, 0)
  if (!is.null(init_contact_mat)) { # User inputted matrix
    # custom assert
    goodegg::assert_eq(sum(any(c("dgeMatrix", "dgCMatrix", "dsCMatrix", "matrix") %in% class(init_contact_mat))),
                       1,
                       message = "If user provides an initial contact matrix, it must be of class
                                  matrix or sparseMatrix (from the Matrix package).")
    if(!is.na(initNC)) {
      warning("The user has provided an initial contact matrix, the initial connection rate will be ignored")
    }

  }

  if (is.matrix(init_contact_mat)) { # catch different type of contact matrix options
    goodegg::assert_square_matrix(init_contact_mat)
    goodegg::assert_symmetric_matrix(init_contact_mat)
    goodegg::assert_dim(init_contact_mat, c(N,N))
    goodegg::assert_eq(sum(diag(init_contact_mat)), 0,
                       message = "Diagonal should be population with 0 to represent self")
    # goodegg::assert_length(unique(rowSums(init_contact_mat)), 1,
    #                         message = "Each node, or individual, must have the same edge density for the initial contact matrix")
    requireNamespace("Matrix")
    conn <- as(init_contact_mat, "sparseMatrix") # rename

  } else if (any(c("dgeMatrix", "dgCMatrix", "dsCMatrix") %in% class(init_contact_mat))) {
    goodegg::assert_square_matrix(as.matrix(init_contact_mat))
    goodegg::assert_symmetric_matrix(as.matrix(init_contact_mat))
    goodegg::assert_dim(as.matrix(init_contact_mat), c(N,N))
    goodegg::assert_eq(sum(diag(as.matrix(init_contact_mat))), 0,
                       message = "Diagonal should be population with 0 to represent self")
    # goodegg::assert_length(unique(rowSums(init_contact_mat)), 1,
    #                         message = "Each node, or individual, must have the same edge density for the initial contact matrix")
    conn <- init_contact_mat # rename

  } else { # user did not input matrix and we will generate one using the degree sequence game
    goodegg::assert_single_int(initNC)
    goodegg::assert_le(initNC, N)
    # simulate initial contact matrix
    conn <- genInitialConnections(initNC, N, sparseMatrix = TRUE)
  }
  #.........................
  # SIR params
  #.........................
  # beta
  if (is.vector(beta)) { # if beta is a vector make matrix but if it is a matrix, leave it alone
    goodegg::assert_vector(beta)
    goodegg::assert_length(beta, N, message = "Beta vector must be of same length as population size, N")
  } else if (is.matrix(beta)) {
    goodegg::assert_square_matrix(beta)
    goodegg::assert_length(nrow(beta), N, message = "Beta matrix must have square dimensions corresponding to population size, NxN")

  } else if (any(c("dgeMatrix", "dgCMatrix", "dsCMatrix") %in% class(beta))) {
    goodegg::assert_square_matrix(as.matrix(beta))
    goodegg::assert_length(nrow(as.matrix(beta)), N, message = "Beta matrix must have square dimensions corresponding to population size, NxN")
  } else {
    stop("Beta must either be a vector, matrix, or sparseMatrix")
  }
  # duration
  goodegg::assert_numeric(dur_I)


  #............................................................
  # initialize and storage
  #...........................................................
  # vector inits
  I_now <- rep(0, N)
  I_now[1:Iseed] <- 1
  S_now <- 1 - I_now
  R_now <- rep(0, N)
  Inds <- 1:N

  # time keeping
  Time_traj <- currtime <- 0.0

  # trajectories for score keeping and plotting
  S_traj <- list(S_now)
  I_traj <- list(I_now)
  R_traj <- list(R_now)
  event_traj <- c("init")
  if (return_contact_matrices) {
    contact_store <- list(conn)
  }

  # lift over beta want B1 - BN as a column copied N times (square matrix with B varying by rows, so columns are replicates)
  beta <- outer(beta, rep(1,N))
  diag(beta) <- 0
  #............................................................
  # run simulation based on rates and events
  #...........................................................

  while (currtime < term_time) {

    # catch - if no infxns, exit loop
    if (sum(I_now) == 0) {
      break
    }

    # transmission rates
    # NB: betaSI has elements beta_i,j * S_j * I_i * connections
    betaSI <- as.matrix( beta * outer(I_now, S_now) * conn ) # need to reduce class to matrix for downstream indexing
    # frequency dependent
    rate_t <- sum(betaSI) # transmission rate depending on overall kinetics

    # recovery rates
    now_dur_I <- (1/dur_I) * I_now
    rate_r <- sum(now_dur_I)

    # Calculate time until each of the events occurs
    event <- c("rewire" = Inf,
               "transmission" = Inf,
               "recovery" = Inf)
    if (rate_t > 0) {
      event[["rewire"]] <- stats::rexp(1, rho)
    }
    if (rate_t > 0) {
      event[["transmission"]] <- stats::rexp(1, rate_t)
    }
    if (rate_r > 0) {
      event[["recovery"]] <- stats::rexp(1, rate_r)
    }

    # Get the event that will occur first, and the time that will take
    next_event <- names(event)[which(event == min(event))]

    # Jump to that time
    currtime <- currtime + event[[next_event]]

    # make change for earliest event
    if (next_event == "rewire") {
      # rewire a pair of nodes
      conn <- rewireNEnodes(adjmat = conn, N = N, sparseMatrix = T)
    } else if (next_event == "transmission") {
      # Choose infector node based on transmission rates
      ind_probs <- rowSums(betaSI) / rate_t # sum is over the columns (i.e. the S_pops)
      parent_pop <- sample(Inds, size = 1, prob = ind_probs)

      # Choose susceptible based on probability of becoming infxn
      ind_probs <- betaSI[parent_pop,] / sum(betaSI[parent_pop,])
      child_pop <- sample(Inds, size = 1, prob = ind_probs)

      # population level updates
      S_now[child_pop] <- S_now[child_pop] - 1
      I_now[child_pop] <- I_now[child_pop] + 1

    } else if (next_event == "recovery") { # If recovery occurs first
      # Choose recovery pop based on recovery rates
      ind_probs <- now_dur_I / rate_r
      recoveree_pop <- sample(Inds, 1, prob = ind_probs)

      # population level updates
      I_now[recoveree_pop] <- I_now[recoveree_pop] - 1
      R_now[recoveree_pop] <- R_now[recoveree_pop] + 1
    }

    # Update lists and time vec
    Time_traj <- c(Time_traj, currtime)
    S_traj <- append(S_traj, list(S_now))
    I_traj <- append(I_traj, list(I_now))
    R_traj <- append(R_traj, list(R_now))
    event_traj <- c(event_traj, next_event)
    if (return_contact_matrices) {
      contact_store <- append(contact_store, list(conn))
    }
  }

  # out
  out <- list(
    S_traj = tidy_traj_out(S_traj),
    I_traj = tidy_traj_out(I_traj),
    R_traj = tidy_traj_out(R_traj),
    Event_traj = event_traj,
    Time_traj = Time_traj,
    runTime = round(Sys.time() - starttime, 2)
  )

  if (return_contact_matrices) {
    out <- append(out, list("contact_store" = contact_store))
  }

  # add S3 class structure
  attr(out, "class") <- "GillespieSIRne"
  return(out)
}

