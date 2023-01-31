#' @title Traditional Gillespie Simulation
#' @param Iseed integer; Number of initial infections in the population
#' @param N integer; population size
#' @param beta numeric vector; Probability of transmission given contact for S_i I_j elements
#' @param dur_I numeric; Duration of infection
#' @param term_time numeric; time termination condition
#' @description
#' Internal function for checking mass action of dynamic gillespie
#' notably drawing uncoupled times to events and events
#' @export
#' @noMd
#' @noRd

tradsim_Gillespie_SIR <- function(Iseed = 1, N = 10,
                                  beta = 1,
                                  dur_I = 1,
                                  term_time = 500) {

  #............................................................
  # assertions
  #...........................................................
  goodegg::assert_single_int(Iseed)
  goodegg::assert_single_int(N)
  goodegg::assert_numeric(dur_I)
  goodegg::assert_numeric(beta)

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
    rate_t <- beta/N * I_now * S_now    # density dependent
    # recovery rates
    rate_r <- dur_I * I_now

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
