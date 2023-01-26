
#' Bulk electrical conductivity of fluid-solid mixtures from Hashin-Shtrikman upper bound
#'
#' @param phi porosity [-]
#' @param sigma_f fluid electrical conductivity [S/m]
#' @param sigma_s Solid electrical conductivity [S/m]
#'
#' @return bulk electrical conductivity [S/m]
#' @export
bulk_ec_HSupper <- function(phi, sigma_f, sigma_s)
{
  sigma_f + (1-phi)/(1/(sigma_s - sigma_f)+phi/(3*sigma_f))
}

#' Bulk electrical conductivity of fluid-solid mixtures based on Archie's law
#'
#' @param phi porosity [-]
#' @param sigma_f fluid electrical conductivity [S/m]
#' @param m cementation exponent [-]
#' @param a tortuosity factor [-]
#'
#' @return bulk electrical conductivity [S/m]
#' @export
bulk_ec_Archie <- function(phi, sigma_f, m, a=1)
{
  sigma_bulk = 1/a*sigma_f*phi^m
  sigma_bulk
}

#' Bulk electrical conductivity of fluid-solid mixtures based on modified Archie's law
#'
#' @param phi porosity [-]
#' @param sigma_f fluid electrical conductivity [S/m]
#' @param sigma_s Solid electrical conductivity [S/m]
#' @param m cementation exponent [-]
#'
#' @return bulk electrical conductivity [S/m]
#' @export
bulk_ec_ModifiedArchie <- function(phi, sigma_f, sigma_s, m)
{
  p = log(1-phi^m)/log(1-phi^m)
  sigma_bulk = sigma_f*phi^m + sigma_s*(1-phi)^p
  sigma_bulk
}

#' Calculate porosity from bulk electrical conductivity using Hashin-Shtrikman upper bound
#'
#' @param sb bulk electrical conductivity [S/m]
#' @param sf fluid electrical conductivity [S/m]
#' @param ss Solid electrical conductivity [S/m]
#'
#' @return porosity [-]
#' @export
phi_HSupper_ec <- function(sb, sf, ss)
{
  phi = (1 - (sb - sf)*(1/(ss - sf)))/(1 + (sb - sf)/(3*sf))
  phi
}

#' Calculate porosity from bulk electrical conductivity using Archie's law
#'
#' @param sigma_bulk bulk electrical conductivity [S/m]
#' @param sigma_f fluid electrical conductivity [S/m]
#' @param m cementation exponent [-]
#'
#' @return porosity [-]
#' @export
phi_Archie_ec <- function(sigma_bulk, sigma_f, m)
{
  phim = sigma_bulk/sigma_f
  phi = phim^(1/m)
  phi
}
