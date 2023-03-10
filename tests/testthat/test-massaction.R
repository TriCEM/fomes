test_that("Network Mass Action vs Traditional Gillespie are Essentially Same", {
  #............................................................
  # Magic Numbers
  #...........................................................
  popsize <- 100
  duration_of_I <- 5
  initial_infxns <- 1
  initNCval <- 99
  betaind <- 1
  rhoconn <- 1e-100

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
                                         initNC = initNCval,
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
  #...........................................................
  testthat::expect_gt(wilcox.test(combouts$NEfinalsize, combouts$TDfinalsize)$p.value, 0.05)
  testthat::expect_gt(wilcox.test(combouts$NEfinaltime, combouts$TDfinaltime)$p.value, 0.05)
  # some made up tolerance for divergence
  mytolsize <- 1
  mytoltime <- 25
  # calculate cheap, not stat robust KL
  fsKL <- KLdivergence_UnifDist(p = combouts$NEfinalsize, q = combouts$TDfinalsize)
  testthat::expect_lt(fsKL, mytolsize)


  ftKL <- KLdivergence_UnifDist(p = round(combouts$NEfinaltime), q = round(combouts$TDfinaltime))
  testthat::expect_lt(ftKL, mytoltime)


})
