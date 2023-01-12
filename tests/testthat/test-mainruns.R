test_that("main model runs", {

  out <- sim_Gillespie_SIR(Iseed = 1, N = 1e3,
                    beta = genRandomBetaMat(1e3),
                    dur_I = 5,
                    rho = 1e-3,
                    term_time = 50)

  testthat::expect_equal(length(out), 4)

})
