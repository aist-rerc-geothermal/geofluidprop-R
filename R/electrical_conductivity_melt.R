
#' Electrical conductivity model for hydrous rhyolitic–granitic melts in Gaillard (2004, EPSL)
#'
#' @param TK Temperature [K]
#' @param pMPa Pressure [MPa]
#' @param wtp_H2O Weight percent of H2O [\%]
#'
#' @return electrical conductivity [S/m]
#' 
#' @export
melt_ec_rhyolite_Gaillard2004_Tpw <- function(TK, pMPa, wtp_H2O)
{
  # unit of wtp_H2O is [%]
  s0 = -78.9*log(wtp_H2O) + 754
  Ea = -2925*log(wtp_H2O) + 64132
  R = 8.31446261815324
  s = s0 * exp((-Ea+2*pMPa)/R/TK)
  s
}

#' Electrical conductivity model for rhyolitic melt in Guo et al (2016, EPSL)
#'
#' @param TK Temperature [K]
#' @param pGPa Pressure [GPa]
#' @param wtp_H2O Weight percent of H2O [\%]
#'
#' @return electrical conductivity [S/m]
#' @export
melt_ec_rhyolite_GuoEtAl2016_Tpw <- function(TK, pGPa, wtp_H2O)
{
  # unit of wtp_H2O is [%]
  logs = 2.983 - 0.0732*wtp_H2O 
  logs = logs - (3528 - 233.8*wtp_H2O + (763-7.5*wtp_H2O^2)* pGPa)  / TK
  10^logs
}
