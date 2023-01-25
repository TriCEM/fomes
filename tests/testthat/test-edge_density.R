test_that("Contact Mat Rewiring is less than or equal to rewiring events", {
    out <- sim_Gillespie_SIR(Iseed = 5, N = 10,
                             beta = rep(0.8, 10),
                             dur_I = 5,
                             rho = 100, # high rewiring rate
                             initNC = 3,
                             term_time = 50,
                             return_contact_matrices = T)
    rewire_count <- sum(out$Event_traj == "rewire")
    unique_contact_mat <- length(unique(out$contact_store))
    testthat::expect_lte(unique_contact_mat, rewire_count)
})


test_that("Consistent Edge Density", {
  initNCit <- 3
  out <- sim_Gillespie_SIR(Iseed = 5, N = 10,
                           beta = rep(0.8, 10),
                           dur_I = 5,
                           rho = 100, # high rewiring rate
                           initNC = initNCit,
                           term_time = 50,
                           return_contact_matrices = T)

  edden <- unique(unlist(lapply(out$contact_store, rowSums)))
  testthat::expect_length(edden, 1)
  testthat::expect_equal(edden, initNCit)
})


test_that("Consistent Switches", {
  # tracker
  count_switch_nodes <- list()
  # do a series of iters
  for (i in 1:100){
    contact_mat_length <- 0
    while (contact_mat_length < 2) {

      out <- sim_Gillespie_SIR(Iseed = 5, N = 10,
                               beta = rep(0.8, 10),
                               dur_I = 5,
                               rho = 100, # high rewiring rate
                               initNC = 3,
                               term_time = 50,
                               return_contact_matrices = T)
      contact_mat_length <- length(unique(out$contact_store))
    } # end while

    swtch <- unique(out$contact_store)
    reldiff <- swtch[[1]][upper.tri(swtch[[1]])] - swtch[[2]][upper.tri(swtch[[2]])]
    reldiff <- sum(reldiff == -1)
    count_switch_nodes <- append(count_switch_nodes, reldiff)
  } # end iters

  #......................
  # make sure we have only 2 nodes (or no switch) per iteration
  #......................
  chck_swtch <- all(unlist(count_switch_nodes) %in% c(0,2))
  testthat::expect_true(chck_swtch)
})








