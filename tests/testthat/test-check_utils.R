test_that("S3 simple overload behaves", {
  # network
  N <- 1e2
  init_contact_mat <- igraph::as_adjacency_matrix(
    igraph::degree.sequence.game(
      out.deg = rep(floor(0.25*N), N), method = "vl"
    ), sparse = F)
  # run of model
  out <- sim_Gillespie_nSIR(Iseed = 1, N = N,
                            beta = rep(1, N),
                            dur_I = 5,
                            rho = 5,
                            init_contact_mat = init_contact_mat,
                            term_time = 500)
  testthat::expect_equal(is(out), "GillespieSIRne")
  testthat::expect_match(capture.output(out)[[1]], "Final Epidemic Size")
  testthat::expect_match(capture.output(out)[[2]], "Run Time")
  testthat::expect_s3_class(out, "GillespieSIRne")
  testthat::expect_length(summary(out), 2)

})
