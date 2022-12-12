


#' @title Gillespie Simulation
#' @param Iseed integer; Number of initial infections in the population
#' @param N integer; population size
#' @param beta double vector; Effective contact rate for S_i I_j elements
#' @description
#' Assuming:
#' \itemize{
#'     closed population
#' }

sim_Gillespie_SIS <- function(Iseed = 1, N = 1e3, beta) {

  #............................................................
  # assertions
  #...........................................................
  assert_single_int(Iseed)
  assert_single_int(N)
  assert_square_matrix(beta)
  assert_dim(beta, c(N,N))

  #............................................................
  # initialize and storage
  #...........................................................
  # vector manips
  I_now <- rep(0, N)
  I_now[1:Iseed] <- 1
  S_now <- abs(1-I_now)
  R_now <- rep(0, N)

  # trajectories for score keeping and plotting
  S_traj <- c()
  I_traj <- c()
  R_traj <- c()

  #............................................................
  # transmission rates and events
  #...........................................................
  betaSI <- beta * outer(S_now, I_now)
  R_t <- sum(betaSI) # transmission rate depending on overall kinetics
  # events
  events <- list()
}
