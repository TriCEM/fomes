
test_that("Event types behaving: aka Rewiring matrices are less than or equal to rewiring events", {
  out <- sim_Gillespie_SIR(Iseed = 1, N = 10,
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


test_that("Initial Network has Consistent Edge Density", {
  initNCit <- 3 # edge density that should be cannon
  firstmat_edge_den <- list()
  # iter it out
  for (i in 1:10) {
    out <- sim_Gillespie_SIR(Iseed = 1, N = 10,
                             beta = rep(0.8, 10),
                             dur_I = 5,
                             rho = 100, # high rewiring rate
                             initNC = initNCit,
                             term_time = 50,
                             return_contact_matrices = T)
    firstmat_edge_den <- append(firstmat_edge_den, out$contact_store[1])
  }

  firstedden <- unique(unlist(lapply(firstmat_edge_den, rowSums)))
  testthat::expect_length(firstedden, 1)
  testthat::expect_equal(firstedden, initNCit)
})


test_that("Rewiring Networks have Consistent Edge Density", {
  initNCit <- 3 # edge density that should be cannon
  edge_den <- list()
  # iter it out
  for (i in 1:10) {
    out <- sim_Gillespie_SIR(Iseed = 1, N = 10,
                             beta = rep(0.8, 10),
                             dur_I = 5,
                             rho = 100, # high rewiring rate
                             initNC = initNCit,
                             term_time = 50,
                             return_contact_matrices = T)
    edge_den <- append(edge_den, unique(out$contact_store))
  }

  edden <- unique(unlist(lapply(edge_den, rowSums)))
  testthat::expect_length(edden, 1)
  testthat::expect_equal(edden, initNCit)
})


test_that("Rewiring Produces Consistent Switches", {
  # tracker
  count_switch_nodes <- list()
  # do a series of iters
  for (i in 1:10){
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










