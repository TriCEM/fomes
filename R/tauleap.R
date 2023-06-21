#' @title SIR NE Gillespie-Tau-Leaping-Algorithm Simulation Model
#' @inheritParams sim_Gillespie_nSIR
#' @param tau numeric; Tau step for efficient (but approximate) stochastic simulation\#' @description Gillespie Simulation
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



sim_tauGillespie_nSIR <- function(Iseed = 1, N = 10,
                                  beta = rep(1, 10),
                                  dur_I = 1,
                                  init_contact_mat = NULL,
                                  rho = 0.05,
                                  tau = 1,
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
  goodegg::assert_single_numeric(tau)
  #.........................
  # initial contact adjacency matrix
  #.........................
  goodegg::assert_single_numeric(rho)
  goodegg::assert_gr(rho, 0)
  goodegg::assert_non_null(init_contact_mat, message = "User must input an initial contact matrix")
  if (!is.null(init_contact_mat)) { # User inputted matrix
    # custom assert
    goodegg::assert_eq(sum(any(c("dgeMatrix", "dgCMatrix", "dsCMatrix", "matrix") %in% class(init_contact_mat))),
                       1,
                       message = "Tbe initial contact matrix must be of class
                                  matrix or sparseMatrix (from the Matrix package).")
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
  event_traj <- list("init" = list("Ic" = I_now, "Rc" = R_now, "Rewirec" = 0))
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

    # Determine number of events in period of time: t, t+tau
    num_rw <- rpois(1, rho * tau)
    num_SI <- rpois(1, rate_t * tau)
    num_IR <- rpois(1, rate_r * tau)

    #..............................
    # implement tau jumped changes
    #..............................
    #......................
    # rewiring changes
    #......................
    if (num_rw > 0) {
      for (i in 1:num_rw) {
        conn <- rewireNEnodes(adjmat = conn, N = N, sparseMatrix = T)
      }
    }
    #......................
    # transmission changes
    #......................
    if (num_SI > 0) {
      # catch if we have exceed our susceptible pool - but b/c we sample w/ replacement based on limited susc by contact matrix, this can be excluded
      # num_SI <- ifelse(num_SI > sum(S_now), sum(S_now), num_SI)
      # Choose infector node based on transmission rates
      ind_probs <- rowSums(betaSI) / rate_t # sum is over the columns (i.e. the S_pops)
      # now draw parent pop where infxns are coming from
      parent_pop <- sample(Inds, size = num_SI, prob = ind_probs, replace = T) # same individual can infect in our time interval
      # Choose susceptible based on probability of becoming infxn
      ind_probs <- betaSI[parent_pop,] / rowSums(betaSI[parent_pop,])
      # but this must be new infxns, so must iteratively update our ind_probs
      child_pop <- rep(NA, num_SI)
      for (i in 1:num_SI) {
        child_pop[i] <- sample(Inds, size = 1, prob = ind_probs[i,], replace = F)
        # ind_probs[, child_pop[i] ] <- 0 would like to say that you can no longer pick this individual but because of contact matrix, we are limited in our susceptible pool even further
      }

      # population level updates
      S_now[child_pop] <- S_now[child_pop] - 1
      I_now[child_pop] <- I_now[child_pop] + 1
    }
    #......................
    # recovery changes
    #......................
    if (num_IR > 0) {
    # Choose recovery pop based on recovery rates
    ind_probs <- now_dur_I / rate_r
    # but catch if we have exceeded our infected pool
    num_IR <- ifelse(num_IR > sum(I_now), sum(I_now), num_IR)
    recoveree_pop <- sample(Inds, size = num_IR, prob = ind_probs)

    # population level updates
    I_now[recoveree_pop] <- I_now[recoveree_pop] - 1
    R_now[recoveree_pop] <- R_now[recoveree_pop] + 1
    }

    #......................
    # Add time
    #......................
    currtime <- currtime + tau


    # Update lists and time vec
    Time_traj <- c(Time_traj, currtime)
    S_traj <- append(S_traj, list(S_now))
    I_traj <- append(I_traj, list(I_now))
    R_traj <- append(R_traj, list(R_now))
    new_event <- list(list("Ic" = num_SI, "Rc" = num_IR, "Rewirec" = num_rw))
    names(new_event) <- paste0("currtime", currtime)
    event_traj <- c(event_traj, new_event)
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
  attr(out, "class") <- "TauGillespieSIRne"
  return(out)
}
