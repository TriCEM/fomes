#' @title SIR-Network Discrete Case Discrete Time Simulation Model
#' @param Iseed integer; Number of initial infections in the population
#' @param N integer; population size
#' @param beta numeric vector, matrix, or sparseMatrix; Probability of transmission given contact for S_i I_j elements
#' @param dur_I numeric; Duration of infection
#' @param term_time numeric; number of steps before time termination condition is met
#' @param init_contact_mat matrix or sparseMatrix; User specified initial contact, or adjacency matrix.
#' The contact matrix must be a symmetric square matrix, with a diagonal of 0 (for self),
#' and have the same edge density, number of connections, for each node, or individual.
#' Default value of NULL results in internal function generating initial contact matrix.
#' @description Gillespie simulation with discrete time and discrete infections
#' @details Assuming:
#' \itemize{
#'    static network
#'     closed population
#' }
#' @noMd
#' @noRd
# make invisible as use-case is for tests not general-purpose

sim_DTDC_nSIR <- function(Iseed = 1, N = 10,
                          beta = rep(1, 10),
                          dur_I = 1,
                          init_contact_mat = NULL,
                          time_steps = 500) {
  #............................................................
  # checks
  #............................................................
  goodegg::assert_single_int(Iseed)
  goodegg::assert_single_int(N)
  goodegg::assert_single_int(time_steps)

  # initial contact adjacency matrix
  goodegg::assert_non_null(init_contact_mat)
  if (!is.null(init_contact_mat)) { # check user inputted matrix
    # custom assert
    goodegg::assert_eq(sum(any(c("dgeMatrix", "dgCMatrix", "dsCMatrix", "matrix") %in% class(init_contact_mat))),
                       1,
                       message = "User must provide an initial contact matrix, and it must be of class
                                  matrix or sparseMatrix (from the Matrix package).")
  }

  if (is.matrix(init_contact_mat)) { # catch different type of contact matrix options
    goodegg::assert_square_matrix(init_contact_mat)
    goodegg::assert_symmetric_matrix(init_contact_mat)
    goodegg::assert_dim(init_contact_mat, c(N,N))
    goodegg::assert_eq(sum(diag(init_contact_mat)), 0,
                       message = "Diagonal should be population with 0 to represent self in your contact matrix")
    # goodegg::assert_length(unique(rowSums(init_contact_mat)), 1,
    #                         message = "Each node, or individual, must have the same edge density for the initial contact matrix")
    conn <- init_contact_mat # rename

  } else if (any(c("dgeMatrix", "dgCMatrix", "dsCMatrix") %in% class(init_contact_mat))) {
    goodegg::assert_square_matrix(as.matrix(init_contact_mat))
    goodegg::assert_symmetric_matrix(as.matrix(init_contact_mat))
    goodegg::assert_dim(as.matrix(init_contact_mat), c(N,N))
    goodegg::assert_eq(sum(diag(as.matrix(init_contact_mat))), 0,
                       message = "Diagonal should be population with 0 to represent self in your contact matrix")
    # goodegg::assert_length(unique(rowSums(init_contact_mat)), 1,
    #                         message = "Each node, or individual, must have the same edge density for the initial contact matrix")
    conn <- as.matrix(init_contact_mat) # remake for indexing below
  } else {
    stop("Inappropriate contact matrix - user must provide an initial contact matrix, and it must be of class
                                  matrix or sparseMatrix (from the Matrix package).")
  }


  #............................................................
  # setup (const, storage, etc)
  #............................................................
  # init state
  currstates <- rep("S", N)
  currstates[1:Iseed] <- "I"

  # time keeping - using 1-based indexing
  currtime <- 1

  # trajectories for score keeping and plotting
  S_traj <- sum(currstates == "S")
  I_traj <- sum(currstates == "I")
  R_traj <- sum(currstates == "R")

  #............................................................
  # core
  #............................................................
  # new infections and recoveries are happening in a binomial process
  while (currtime <= time_steps & sum(currstates == "I") != 0) {  #termination conditions
    #......................
    # each time step has chance to infect and recover given node & state
    #......................
    for (i in 1:length(currstates)) {
      if (currstates[i] == "S") {
        betai <- beta[i] # store Suscind beta
        ptedges <- which(conn[i,] == 1) # potential transmission edges
        nt <- sum(currstates[ptedges] == "I")
        if (nt > 0) { # number of infectors that have a trial to infect S
          if( any( rbinom(1, size = nt, prob = 1 - exp(-betai)) ) ) { # if an infection has occured
            currstates[i] <- "I"
          }
        }
      } # end S condition
      else if (currstates[i] == "I") {
        if (rbinom(1,1, 1 - exp(-1/dur_I))){
          currstates[i] <- "R"
        }
      } # end recover
    } # end state for loop
    #......................
    # update counters before next for loop
    #......................
    S_traj <- c(S_traj, sum(currstates == "S"))
    I_traj <- c(I_traj, sum(currstates == "I"))
    R_traj <- c(R_traj, sum(currstates == "R"))
    currtime <- currtime + 1
  }
  #............................................................
  # tidy & out
  #............................................................
  ret <- tibble::tibble(
    time = 1:currtime,
    Susc = S_traj,
    Infxns = I_traj,
    Recov = R_traj)
  return(ret)
}



#' @title Traditional Gillespie Simulation
#' @param Iseed integer; Number of initial infections in the population
#' @param N integer; population size
#' @param beta numeric vector; Probability of transmission given contact for S_i I_j elements
#' @param dur_I numeric; Duration of infection
#' @param term_time numeric; time termination condition
#' @description
#' Internal function for checking mass action of dynamic gillespie
#' notably drawing uncoupled times to events and events
#' @noMd
#' @noRd
# make invisible as use-case is for tests not general-purpose

sim_Gillespie_SIR <- function(Iseed = 1, N = 10,
                              beta = 1,
                              dur_I = 1,
                              term_time = 500) {

  #............................................................
  # assertions
  #...........................................................
  goodegg::assert_single_int(Iseed)
  goodegg::assert_single_int(N)
  goodegg::assert_single_numeric(dur_I)
  goodegg::assert_single_numeric(beta)

  #............................................................
  # initialize and storage
  #...........................................................
  # vector inits
  I_traj <- I_now <- Iseed
  S_traj <- S_now <- N - Iseed
  R_traj <- R_now <- 0
  Inds <- 1:N
  # time keeping
  t <- time <- 0.0

  #............................................................
  # run simulation based on rates and events
  #...........................................................

  while (time < term_time) {

    # catch - if no infxns, exit loop
    if (I_now == 0) {
      break
    }

    # transmission rates
    rate_t <- beta * 1/N * I_now * S_now    # density dependent
    # recovery rates
    rate_r <- (1/dur_I) * I_now

    # Calculate time until each of the events occurs
    event <- c("transmission" = Inf,
               "recovery" = Inf)
    if (rate_t > 0) {
      event[["transmission"]] <- stats::rexp(1, rate_t)
    }
    if (rate_r > 0) {
      event[["recovery"]] <- stats::rexp(1, rate_r)
    }

    # Get the event that will occur first, and the time that will take
    next_event <- names(event)[which(event == min(event))]

    # Jump to that time
    time <- time + event[[next_event]]

    if (next_event == "transmission") {
      # population level updates
      S_now <- S_now - 1
      I_now <- I_now + 1

    } else if (next_event == "recovery") { # If recovery occurs first
      I_now <- I_now - 1
      R_now <- R_now + 1
    }

    # Update lists and time vec
    t <- c(t, time)
    S_traj <- c(S_traj, S_now)
    I_traj <- c(I_traj, I_now)
    R_traj <- c(R_traj, R_now)
  }

  # out
  out <- data.frame(
    Time = t,
    Susc = S_traj,
    Infxn = I_traj,
    Recov = R_traj
  )
  return(out)
}


