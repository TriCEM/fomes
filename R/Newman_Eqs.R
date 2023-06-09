
#' @title Probability Generating Function for Degree Distribution
#' @param x numeric vector; discerete random variable (ie vertices)
#' @param degdist numeric vector; degree distribution
#' @description Newman 2002
#' @details
#' @returns
#' @noRd

get_newman_G0x <- function(x, degdist) {
  # checks
  goodegg::assert_vector_numeric(x)
  goodegg::assert_vector_numeric(degdist)

  # core
  p <- 0 # init
  for (k in (seq_along(degdist) - 1)) { # \sum_0^inf (but actually k) p_k * x^k
    p <- p + degdist[k + 1] * x^k
  }
  return(p)
}

#' @title First Moment of G0: Mean Vertex Degrees
#' @inheritParams get_G0
#' @description Newman 2002
#' @details
#' @returns
#' @noRd

get_newman_G0x_firstmoment <- function(x, degdist) {
  # checks
  goodegg::assert_vector_numeric(x)
  goodegg::assert_vector_numeric(degdist)
  # core
  p <- 0 # init
  for (k in (seq_along(degdist) - 1)) { # this G_0'
    p <- p + k * degdist[k + 1] * x^(k-1)
  }
  return(p)
}


#' @title Second Moment of G0: Variance of Vertex Degrees
#' @inheritParams get_G0
#' @description Newman 2002
#' @details
#' @returns
#' @noRd

get_newman_G0x_secondmoment <- function(x, degdist) {
  # checks
  goodegg::assert_vector_numeric(x)
  goodegg::assert_vector_numeric(degdist)
  # core
  # core
  p <- 0 # init
  for (k in (seq_along(degdist) - 1)) { # this is G_0''
    p <- p + (k-1) * k * degdist[k + 1] * x^(k - 2)
  }
  return(p)
}




#' @title Probability Generating Function for Occupied Edges Leaving a Vertex Arrived at by following a Randomly Chosen Edge
#' @inheritParams get_G0
#' @description Newman 2002
#' @details
#' @returns
#' @noRd

get_newman_G1x <- function(x, degdist) {
  # checks
  goodegg::assert_vector_numeric(x)
  goodegg::assert_vector_numeric(degdist)
  # core
  p <- 0 # init
  for (k in (seq_along(degdist) - 1)) { # this is G_0'
    p <- p + k * degdist[k + 1] * x^(k - 1)
  }
  # z average vertex degree
  z <- get_newman_G0x_firstmoment(1, degdist)
  p <- p / z
  # out
  return(p)
}

#' @title Newman Epidemic Transition
#' @inheritParams get_G0
#' @description Newman 2002
#' @details
#' @returns
#' @noRd
get_newman_Tc <- function(degdist) {
  # checks
  goodegg::assert_vector_numeric(degdist)
  # core
  Tc <- get_newman_G0x_firstmoment(1, degdist) / get_newman_G0x_secondmoment(1, degdist)
  # out
  return(Tc)
}

#' @title Newman Mean Final Size below Tc
#' @inheritParams get_G0
#' @param Transm
#' @description Newman 2002
#' @details
#' @returns
#' @noRd
get_newman_mean_epidemic_size_below_Tc <- function(Transm, degdist) {
  # checks
  goodegg::assert_single_numeric(Transm)
  goodegg::assert_vector_numeric(degdist)
  # catch
  if (Transm > get_newman_Tc(degdist)) {
    stop("Transmission exceeds epidemic transition point and the equation is no longer valid")
  }

  # core
  num <- Transm * get_newman_G0x_firstmoment(1, degdist)
  denom <- 1 - Transm * (1/get_newman_Tc(degdist))
  s <- 1 + num/denom
  # out
  return(s)
}





#' @title Solving Self Consistency Relationship
#' @param transm numeric; Transmissibility of process under consideration
#' @inheritParams get_G0
#' @param initu numeric; initial value of u bounded by [0,1]
#' @param iters integer; number of steps to consider for evaluating u
#' @param tol numeric; convergence tolerance/escape for u
#' @description
#' @details Newman 2002: u = H1(1; T), and is the solution of the self-consistency relation u = G1(u; T)
#' @returns
#' @noRd

solve_newman_u <- function(transm, degdist, initu = 0.5, iters = 1e6, tol = 1e-5) {
  # checks
  goodegg::assert_single_numeric(transm)
  goodegg::assert_bounded(transm, left = 0, right = 1, inclusive_left = T, inclusive_right = T,
                          message = "Transmission must be between 0 and 1")
  goodegg::assert_vector_numeric(degdist)
  goodegg::assert_single_numeric(tol)
  goodegg::assert_single_int(iters)
  goodegg::assert_single_numeric(initu)
  goodegg::assert_bounded(initu, left = 0, right = 1, inclusive_left = T, inclusive_right = T)
  # core
  u <- initu  # starting value
  for (i in 1:iters) {
    u_new <- get_newman_G1x(1 - transm + u *transm, degdist)
    if (abs(u_new - u) < tol) break  # break if change is very small, reached self-consistency
    u <- u_new
  }
  # out
  return(u)
}


#' @title Transmissibility from Disease Causing Contact Rate, Duration of Illness from Newman
#' @param r numeric vector; rate of disease causing contacts between individuals i and j
#' @param tau numeric vector; duration of illness per by individiaul i
#' @description Newman 2002
#' @details
#' @returns
#' @noRd

get_newman_transmissiblity <- function(taui, rij) {
  #......................
  # checks
  #......................
  goodegg::assert_numeric(taui)
  goodegg::assert_numeric(rij)

  # rij is average rate of disease-causing contacts between infected i and susceptible j
  # tau is the duration of illness for individual i
  # like Newman, will assume simple uniform dist for 0 <= x <= rmax
  # although Newman summation goes to Inf, we only have support from a-b on our Unif Dist
  #......................
  # integration
  #......................
  P_r <- function(r){
    dunif(r, min = 0, max = max_r)
  }
  # tau is the duration of illness for individual i
  # assuming discrete time for simplicity and uniform dist bouned by max time
  # therefore although Newman summation goes to Inf, we only have support from a-b on our Unif Dist
  P_t <- function(t){
    dunif(t, min = 0, max = max_t)
  }
  # set maxes from empiric data
  max_r <- max(rij)
  mat_t <- max(taui)

  # Define the integrand
  integrand <- function(rt) {
    r <- rt[1]
    t <- rt[2]
    return(P_r(r) * P_t(t) * exp(-r*t))
  }

  # Perform the double integral
  retinteg <- cubature::adaptIntegrate(integrand,
                                       lowerLimit = c(0, 0),
                                       upperLimit = c(Inf, Inf))

  #......................
  # out
  #......................
  return(1 - retinteg$integral)
}



#' @title Transmissibility Approximation
#' @param r numeric vector; rate of disease causing contacts between individuals i and j
#' @param tau numeric vector; duration of illness per by individiaul i
#' @description
#' @details
#' @returns
#' @noRd

get_approx_transmissiblity <- function(taui, rij) {
  #......................
  # checks
  #......................
  goodegg::assert_numeric(taui)
  goodegg::assert_numeric(rij)

  # discretize both of these for simiplicity
  timeperiods <- as.numeric(names(table(taui)))
  timewi <- as.numeric(table(taui))/length(taui)
  contactrates <- as.numeric(names(table(rij)))
  connwi <- as.numeric(table(rij))/length(rij)
  # init
  Tout <- 0
  for (i in 1:length(contactrates)) {
    for (j in 1:length(timeperiods)) {
      #Tout <- Tout + exp(- (contactrates[i] * connwi[i]) * (timeperiods[j] * timewi[j]))
      Tout <- Tout + ((1-contactrates[i])^timeperiods[j])*connwi[i]*timewi[j]
    }
  }
  return(1 - Tout)
}



#' @title Mean Final Outbreak Size by Degree Distribution from Newman
#' @inheritParams get_G0
#' @inheritParams get_newman_transmissiblity
#' @param transmApprox whether to use an approximation for transmission or the exact integral solution
#' @description
#' @details
#' @returns
#' @noRd

get_newman_mean_final_epidemic_size <- function(graph, taui, rij,
                                                transmApprox = TRUE,
                                                initu = 0.5, iters = 1e6, tol = 1e-5) {
  #......................
  # checks
  #......................
  goodegg::assert_numeric(taui)
  goodegg::assert_numeric(rij)
  goodegg::assert_eq(class(graph), "igraph",
                     message = "Graph must be of class igraph")
  goodegg::assert_single_numeric(tol)
  goodegg::assert_single_int(iters)
  goodegg::assert_single_numeric(initu)
  goodegg::assert_bounded(initu, left = 0, right = 1, inclusive_left = T, inclusive_right = T)


  #......................
  # core
  #......................
  degdist <- igraph::degree_distribution(graph = graph)

  if (transmApprox) {
    transmnow <- get_approx_transmissiblity(taui, rij)
  } else {
    transmnow <- get_newman_transmissiblity(taui, rij)
  }
  unow <- solve_newman_u(transm = transmnow, degdist = degdist,
                         initu = initu, iters = iters, tol = tol)
  ret <- 1 - get_newman_G0x(x = unow, degdist = degdist)

  #......................
  # out
  #......................
  return(ret)
}

