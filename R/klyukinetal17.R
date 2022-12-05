
#' Dynamic Viscosity of Salinewater, Function of Density, Temperature, and Salinity
#'
#' @references Klyukin, Y.I., Lowell R.P., Bodnar, R.J. (2017) A revised empirical
#' model to calculate the dynamic viscosity of H2O-NaCl fluids at elevated
#' temperatures and pressures (<=1000 C, <=500 MPa, 0-100 wt % NaCl). Fluid Phase
#' Equilibria 433, 193-205.
#'
#' @param rho Fluid density [kg/m3]
#' @param TK  Temperature [K]
#' @param x   Salinity in mole fraction [mol/mol]
#' @return Dynamic viscosity [Pa s]
#'
#' @export
klyukinetal2017_H2O_NaCl_viscosity_rhoTx <- function(rho, TK, x)
{
  ret <- .C("R_klyukinetal2017_H2O_NaCl_viscosity_rhoTx", rho=as.double(rho), TK=as.double(TK), x=as.double(x), out=as.double(1.0))
  return(ret$out)
}
