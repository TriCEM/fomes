test_that("Network continuous time and network discrete time are similar", {
  #............................................................
  # Magic Numbers
  #...........................................................
  N <- 100
  duration_of_I <- 5
  initial_infxns <- 1
  betaind <- 0.15
  rhoconn <- 1e-100
  N <- 1e2
  conmat <- igraph::as_adjacency_matrix(
    igraph::degree.sequence.game(
      out.deg = rep(floor(0.25*N), N), method = "vl"
    ), sparse = F)
  diag(conmat) <- 0

  # storage
  niters <- 1e3
  combouts <- as.data.frame(matrix(NA, nrow = niters, ncol = 3))
  colnames(combouts) <- c("iter",
                          "ContTimefinalsize",
                          "DiscTimefinalsize")
  combouts$iter <- 1:niters

  #............................................................
  # run through
  #...........................................................
  for (i in 1:niters) {
    #......................
    # cont network
    #......................
    CTnSIR <- suppressWarnings( fomes::sim_Gillespie_nSIR(Iseed = initial_infxns,
                                                          N = N,
                                                          beta = rep(betaind, N),
                                                          dur_I = duration_of_I,
                                                          rho = rhoconn,
                                                          init_contact_mat = conmat,
                                                          term_time = Inf)
    ) # warnings from igraph that we saturated network

    tidyCTnSIR <- summary(CTnSIR)
    # storage
    combouts[i, "ContTimefinalsize"] <- tidyCTnSIR$FinalEpidemicSize

    #......................
    # discrete network
    #......................
    tradSIR <- fomes:::sim_DTDC_nSIR(Iseed = initial_infxns,
                                     N = N,
                                     beta = rep(betaind, N),
                                     dur_I = duration_of_I,
                                     init_contact_mat = conmat,
                                     time_steps = Inf)

    # storage

    combouts[i, "DiscTimefinalsize"] <- N - tradSIR$Susc[length(tradSIR$Susc)]
  }

  #............................................................
  # analyze results
  # know that this a parameter space that is not in a phase transition region
  # so primarily should be all or nothing, and we can mostly treat it as a discrete
  # count space, and do the a cheap KL
  #...........................................................
  tCT <- as.data.frame(table(combouts$ContTimefinalsize), stringsAsFactors = F)
  tDT <- as.data.frame(table(combouts$DiscTimefinalsize), stringsAsFactors = F)
  difftab <- dplyr::full_join(tCT, tDT, by = "Var1") %>%
    dplyr::mutate(Var1 = as.numeric(Var1)) %>%
    dplyr::arrange(Var1)

  # do cheap KL
  p <- difftab$Freq.x
  q <- difftab$Freq.y
  kl_div <- sum(p * log(p / q, base = exp(1)), na.rm = T)
  testthat::expect_lt(kl_div, 100) # higher tolerance given 0s

   #......................
   # man viz
   #......................
   p1d <- tibble::tibble(fs = combouts$ContTimefinalsize, mod = "Cont")
   p2d <- tibble::tibble(fs = combouts$DiscTimefinalsize, mod = "Disc")
  dplyr::bind_rows(p1d, p2d)%>%
     ggplot( aes(x = fs, fill = mod)) +
     geom_histogram(alpha = 0.4, position = 'identity') +
     scale_fill_manual(values = c("#fee391", "#a6bddb")) +
     theme_linedraw()




})
