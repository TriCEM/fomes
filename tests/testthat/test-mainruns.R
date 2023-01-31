test_that("main model runs", {

  out <- sim_Gillespie_SIR(Iseed = 1, N = 1e3,
                    beta = rep(1, 1e3),
                    dur_I = 5,
                    rho = 5,
                    initNC = 3,
                    term_time = 500)

  testthat::expect_equal(length(out), 6)


  # traditional model for internal checks and comparisons
 out <- fomes:::tradsim_Gillespie_SIR(Iseed = 1, N = 10,
                        beta = 1,
                        dur_I = 1,
                        term_time = 500)
 testthat::expect_type(out, "list")
 testthat::expect_equal(class(out), "data.frame")
})
