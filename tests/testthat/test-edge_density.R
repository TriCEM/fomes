test_that("Edge Density Stays Consistent", {

    out <- sim_Gillespie_SIR(Iseed = 5, N = 10,
                             beta = rep(0.8, 10),
                             dur_I = 5,
                             rho = 0.5,
                             initNC = 3,
                             term_time = 50,
                             return_contact_matrices = T)
    out$Event_traj



})
