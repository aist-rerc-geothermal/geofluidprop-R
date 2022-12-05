context("IAPWS95")

test_that("IAPWS95_rho_pT", {
  rho <- iapws95_rho_pT(p=1e6, T=293.15)
  expect_equal(998.6184, rho, tolerance = 1e2*1e-3)
})

