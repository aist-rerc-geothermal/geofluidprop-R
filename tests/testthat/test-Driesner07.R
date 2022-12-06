context("DRIESNER07")

test_that("DRIESNER07_rho_pT", {
  rho <- driesner07_H2O_NaCl_rho_pTx(p=10e6, TK=273.15+200, x=H2ONaCl_massfrac_to_x(0.1e-2))
  print(rho)
  expect_equal(871.8296, rho, tolerance = 1e2*1e-3)
})

