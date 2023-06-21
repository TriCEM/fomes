
test_that("Event types behaving: aka Rewiring matrices are less than or equal to rewiring events", {
  # network
  N <- 1e2
  init_contact_mat <- igraph::as_adjacency_matrix(
    igraph::degree.sequence.game(
      out.deg = rep(floor(0.25*N), N), method = "vl"
    ), sparse = F)
  # run of model
   out <- sim_Gillespie_nSIR(Iseed = 1, N = N,
                           beta = rep(0.8, N),
                           dur_I = 5,
                           rho = 100, # high rewiring rate
                           init_contact_mat = init_contact_mat,
                           term_time = 50,
                           return_contact_matrices = F)
  rewire_count <- sum(out$Event_traj == "rewire")
  unique_contact_mat <- length(unique(out$contact_store))
  testthat::expect_lte(unique_contact_mat, rewire_count)
})


test_that("Initial Network has Consistent MODE/MEDIAN Edge Density", {
  # network
  N <- 1e2
  mydef <- 3
  init_contact_mat <- igraph::as_adjacency_matrix(
    igraph::degree.sequence.game(
      out.deg = rep(mydef, N), method = "vl" # edge density that should be cannon
    ), sparse = F)

  # check
  firstmat_edge_den <- list()
  # iter it out
  for (i in 1:10) {
    out <- sim_Gillespie_nSIR(Iseed = 1, N = N,
                             beta = rep(0.8, N),
                             dur_I = 5,
                             rho = 1,
                             init_contact_mat = init_contact_mat,
                             term_time = 50,
                             return_contact_matrices = T)
    firstmat_edge_den <- append(firstmat_edge_den, out$contact_store[1])
  }
  getmode <- function(x) {
    ret <- names(table(x))[which(table(x) == max(table(x)))]
    ret <- as.numeric(ret)
    return(ret)
  }
  # liftover
  firstmat_edge_den <- lapply(firstmat_edge_den, as.matrix)
  firstmode <- getmode(unlist(lapply(firstmat_edge_den, rowSums)))
  firstmed <- median(unlist(lapply(firstmat_edge_den, rowSums)))
  testthat::expect_equal(firstmode, mydef)
  testthat::expect_equal(firstmed, mydef)
})



test_that("Rewiring Produces Consistent Switches", {
  # network
  N <- 1e2
  init_contact_mat <- igraph::as_adjacency_matrix(
    igraph::degree.sequence.game(
      out.deg = rep(3, N), method = "vl" # edge density that should be cannon
    ), sparse = F)



  # tracker
  count_switch_nodes <- list()
  # do a series of iters
  for (i in 1:10){
    contact_mat_length <- 0
    while (contact_mat_length < 2) {

      out <- sim_Gillespie_nSIR(Iseed = 5, N = N,
                               beta = rep(0.8, N),
                               dur_I = 5,
                               rho = 10, # medium rewiring rate
                               init_contact_mat = init_contact_mat,
                               term_time = 50,
                               return_contact_matrices = T)
      contact_mat_length <- length(unique(out$contact_store))
    } # end while

    swtch <- unique(out$contact_store)
    swtch <- lapply(swtch, as.matrix)
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










