
#' Water density using IAPWS-95
#'
#' @param pPa  pressure [Pa]
#' @param TK  Temperature [K]
#' @return water density [kg/m^3]
#'
#' @export
iapws95_rho_pT <- function(pPa, TK)
{
  ret <- .C("R_iapws95_rho_pT", pPa=as.double(pPa), TK=as.double(TK), out=as.double(1.0))
  return(ret$out)
}
