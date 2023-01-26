
#' Electrical conductivity model for dry granite based on Olhoeft (1981, JGR) figure
#'
#' @param TK Temperature [K]
#'
#' @return electrical conductivity of dry granite [S/m]
#' @export
granite_ec_Olhoeft1981_T <- function(TK)
{
  TC = TK -273.15
  log_r = 7.34406E-06*TC^2 - 1.95002E-02*TC + 1.35479E+01
  1/(10^log_r)
}

