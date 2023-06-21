test_that("Gillespie Network CT model runs", {
  N <- 1e2
  init_contact_mat <- igraph::as_adjacency_matrix(
    igraph::degree.sequence.game(
      out.deg = rep(floor(0.25*N), N), method = "vl"
    ), sparse = F)

  out <- sim_Gillespie_nSIR(Iseed = 1, N = N,
                            beta = rep(0.5, N),
                            dur_I = 5,
                            rho = 5,
                            init_contact_mat = init_contact_mat,
                            term_time = 500)
  # tests
  testthat::expect_equal(length(out), 6)
  testthat::expect_type(out, "list")
})




test_that("tradiational Gillespie MA model runs", {

  out <- sim_Gillespie_SIR(Iseed = 1, N = 1e2,
                           beta = 0.5,
                           dur_I = 5,
                           term_time = Inf)
  # tests
  testthat::expect_equal(length(out), 4)
  testthat::expect_equal(class(out), "data.frame")

})


test_that("Gillespie Network DTDC model runs", {
  N <- 1e2
  init_contact_mat <- igraph::as_adjacency_matrix(
    igraph::degree.sequence.game(
      out.deg = rep(floor(0.25*N), N), method = "vl"
    ), sparse = F)


  out <- fomes:::sim_DTDC_nSIR(Iseed = 1, N = N,
                               beta = rep(1, N),
                               dur_I = 5,
                               init_contact_mat = init_contact_mat,
                               time_steps = 500)

  # tests
  testthat::expect_equal(length(out), 4)
  testthat::expect_type(out, "list")
})




test_that("Gillespie Network Tau Leaping model runs", {
  N <- 1e2
  init_contact_mat <- igraph::as_adjacency_matrix(
    igraph::degree.sequence.game(
      out.deg = rep(floor(0.25*N), N), method = "vl"
    ), sparse = F)

  out <- sim_tauGillespie_nSIR(Iseed = 1, N = N,
                               beta = rep(0.5, N),
                               dur_I = 5,
                               rho = 1,
                               tau = 1,
                               init_contact_mat = init_contact_mat,
                               term_time = 500)
  # tests
  testthat::expect_equal(length(out), 6)
  testthat::expect_type(out, "list")
})
