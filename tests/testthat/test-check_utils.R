test_that("S3 simple overload behaves", {
  out <- sim_Gillespie_SIR(Iseed = 1, N = 1e3,
                           beta = rep(1, 1e3),
                           dur_I = 5,
                           rho = 5,
                           initNC = 3,
                           term_time = 500)
  testthat::expect_equal(is(out), "GillespieSIRne")
  testthat::expect_match(capture.output(out)[[1]], "Final Epidemic Size")
  testthat::expect_match(capture.output(out)[[2]], "Run Time")
  testthat::expect_s3_class(out, "GillespieSIRne")
  testthat::expect_length(summary(out), 2)

})
