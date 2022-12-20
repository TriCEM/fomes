#' @title Gillespie Simulation
#' @param Iseed integer; Number of initial infections in the population
#' @param N integer; population size
#' @param beta numeric vector; Effective contact rate for S_i I_j elements
#' @param dur_I numeric; Duration of infection
#' @param term_time numeric; time termination condition
#' @description
#' Assuming:
#' \itemize{
#'     closed population
#' }
#' @export

sim_Gillespie_SIR <- function(Iseed = 1, N = 1e3,
                              beta = genRandomBetaMat(N),
                              dur_I = 1,
                              term_time = 500) {

  #............................................................
  # assertions
  #...........................................................
  goodegg::assert_single_int(Iseed)
  goodegg::assert_single_int(N)
  goodegg::assert_numeric(dur_I)
  goodegg::assert_square_matrix(beta)
  goodegg::assert_dim(beta, c(N,N))

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
  i <- 2 # i of 1 is our initial conditions and time of 0
  t <- time <- 0.0

  # trajectories for score keeping and plotting
  S_traj <- list(S_now)
  I_traj <- list(I_now)
  R_traj <- list(R_now)

  #............................................................
  # run simulation based on rates and events
  #...........................................................

  while (time < term_time) {

    # catch - if no infxns, exit loop
    if (sum(I_now) == 0) {
      break
    }

    # transmission rates
    betaSI <- beta * outer(S_now, I_now) # betaSI has elements beta_i,j * S_j * I_i
    rate_t <- sum(betaSI) # transmission rate depending on overall kinetics

    # recovery rates
    now_dur_I = dur_I * I_now
    rate_r <- sum(now_dur_I)

    # Calculate time until each of the events occurs
    event <- c("transmission" = Inf,
                "recovery" = Inf)
    if (rate_t > 0) {
      event[["transmission"]] <- rexp(1, 1 / rate_t)
    }
    if (rate_r > 0) {
      event[["recovery"]] <- rexp(1, 1 / rate_r)
    }

    # Get the event that will occur first, and the time that will take
    next_event <- names(event)[which(event == min(event))]

    # Jump to that time
    time <- time + event[[next_event]]

    if (next_event == "transmission") {
      #Choose infector pop based on transmission rates
      ind_probs <- rowSums(betaSI) / rate_t # sum is over the columns (i.e. the S_pops)
      parent_pop <- sample(Inds, size = 1, prob = ind_probs)

      ind_probs <- betaSI[parent_pop,] / sum(betaSI[parent_pop,]) # going across columns
      child_pop <- sample(Inds, size = 1, prob = ind_probs)

      # population level updates
      S_now[child_pop] <- S_now[child_pop] - 1
      I_now[child_pop] <- I_now[child_pop] + 1 # had parent_pop (incorrectly) here for some reason?

    } else if (next_event == "recovery") { # If recovery occurs first
      # Choose recovery pop based on recovery rates
      ind_probs <- now_dur_I / rate_r
      recoveree_pop <- sample(Inds, 1, p = ind_probs)

      # population level updates
      I_now[recoveree_pop] <- I_now[recoveree_pop] - 1
      R_now[recoveree_pop] <- R_now[recoveree_pop] + 1
    }

    # Update lists and time vec
    t <- c(t, time)
    S_traj <- append(S_traj, list(S_now))
    I_traj <- append(I_traj, list(I_now))
    R_traj <- append(R_traj, list(R_now))
    # update counter
    i <- i+ 1

  }

  # out
  out <- list(
    S_traj = tidy_traj_out(S_traj),
    I_traj = tidy_traj_out(I_traj),
    R_traj = tidy_traj_out(R_traj),
    time = t
  )
  return(out)
}

