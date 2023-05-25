test_that("Network Mass Action vs Traditional Gillespie are Essentially Same", {
  #............................................................
  # Magic Numbers
  #...........................................................
  popsize <- 100
  duration_of_I <- 5
  initial_infxns <- 1
  betaind <- 1
  rhoconn <- 1e-100
  conmat <- matrix(1, popsize, popsize)
  diag(conmat) <- 0

  # storage
  niters <- 1e3
  combouts <- as.data.frame(matrix(NA, nrow = niters, ncol = 5))
  colnames(combouts) <- c("iter",
                          "NEfinalsize", "NEfinaltime",
                          "TDfinalsize", "TDfinaltime")
  combouts$iter <- 1:niters

  #............................................................
  # run through
  #...........................................................
  for (i in 1:niters) {
    #......................
    # dynamic network
    #......................
    NEdynSIR <- suppressWarnings( fomes::sim_Gillespie_SIR(Iseed = initial_infxns,
                                         N = popsize,
                                         beta = rep(betaind, popsize),
                                         dur_I = duration_of_I,
                                         rho = rhoconn,
                                         init_contact_mat = conmat,
                                         term_time = Inf)
    ) # warnings from igraph that we saturated network

    tidyNEdynSIR <- summary(NEdynSIR)
    # storage
    combouts[i, "NEfinalsize"] <- tidyNEdynSIR$FinalEpidemicSize
    combouts[i, "NEfinaltime"] <- max(tidyNEdynSIR$CumEvents$Time)


    #......................
    # traditional model
    #......................
    tradSIR <- fomes:::tradsim_Gillespie_SIR(Iseed = initial_infxns,
                                             N = popsize,
                                             beta = betaind,
                                             dur_I = duration_of_I,
                                             term_time = Inf)

    # storage
    combouts[i, "TDfinalsize"] <- sum(tradSIR[1,c("Susc", "Infxn")]) - tradSIR[nrow(tradSIR), "Susc"]
    combouts[i, "TDfinaltime"] <- tradSIR[nrow(tradSIR), "Time"]

  }

  #............................................................
  # analyze results
  # know that this a parameter space that is not in a phase transition region
  # so primarily should be all or nothing, and we can mostly treat it as a discrete
  # count space, and do the a cheap KL
  #...........................................................
  tNE <- as.data.frame(table(combouts$NEfinalsize), stringsAsFactors = F)
  tMA <- as.data.frame(table(combouts$TDfinalsize), stringsAsFactors = F)
  difftab <- dplyr::full_join(tNE, tMA, by = "Var1") %>%
    dplyr::mutate(Freq.x = ifelse(is.na(Freq.x), 0, Freq.x),
                  Freq.y = ifelse(is.na(Freq.y), 0, Freq.y)) %>%
    dplyr::mutate(Var1 = as.numeric(Var1)) %>%
    dplyr::arrange(Var1)

  # do cheap KL
  p <- difftab$Freq.x
  p <- ifelse(p == 0, 1e-10, p)
  q <- difftab$Freq.y
  q <- ifelse(q == 0, 1e-10, q)
  kl_div <- sum(p * log(p / q, base = exp(1)))
  testthat::expect_lt(kl_div, 500) # higher tolerance given 0s





})
