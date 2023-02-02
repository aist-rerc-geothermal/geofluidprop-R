
#' Critical pressuer of H2O-NaCl fluid as a function of temperature
#'
#' Returns critical pressure of H2O-NaCl fluid at the given temperature
#'
#' @param TK Temperature [K]
#'
#' @return Critical pressure [Pa]
#'
#' @export
driesner07_H2O_NaCl_pc_T <- function(TK)
{
  ret <- .C("R_driesner07_H2O_NaCl_pc_T", TK=as.double(TK), pc=as.double(1.0))
  return(ret$pc)
}

#' Critical pressuer of H2O-NaCl fluid as a function of temperature
#' using saturation pressure of water below critical temperature
#'
#' Returns critical pressure of H2O-NaCl fluid at the given temperature
#'
#' @param TK Temperature [K]
#'
#' @return Critical pressure [Pa]
#'
#' @export
driesner07_H2O_NaCl_pc_T2 <- function(TK)
{
  ret <- .C("R_driesner07_H2O_NaCl_pc_T2", TK=as.double(TK), pc=as.double(1.0))
  return(ret$pc)
}

#' Critical NaCl composition of H2O-NaCl fluid as a function of temperature
#'
#' Returns critical composition of H2O-NaCl fluid at the given temperature
#'
#' @param TK Temperature [K]
#'
#' @return Critical composition (mole fraction) of NaCl [mol/mol]
#'
#' @export
driesner07_H2O_NaCl_xc_T <- function(TK)
{
  ret <- .C("R_driesner07_H2O_NaCl_xc_T", TK=as.double(TK), xc=as.double(1.0))
  return(ret$xc)
}

#' Critical temperature of H2O-NaCl fluid as a function of NaCl composition
#'
#' Returns critical temperature of H2O-NaCl fluid at the given composition
#'
#' @param x Composition (mole fraction) of NaCl [mol/mol]
#' @return Temperature [K]
#'
#' @export
driesner07_H2O_NaCl_Tc_x <- function(x)
{
  ret <- .C("R_driesner07_H2O_NaCl_Tc_x", x=as.double(x), Tc=as.double(1.0))
  return(ret$Tc)
}

#' Phase id of H2O-NaCl fluid as a function of temperature, pressure, and NaCl composition
#'
#' Returns phase id of H2O-NaCl fluid at the given temperature, pressure, and bulk salinity
#'
#' @param TK Temeprature [K]
#' @param p Pressure [Pa]
#' @param x Bulk composition (mole fraction) of NaCl [mol/mol]
#' @return phase id
#'
#' @export
driesner07_H2O_NaCl_phase_Tpx <- function(TK,p,x)
{
  ret <- .C("R_driesner07_H2O_NaCl_phase_pTx", p=as.double(p), TK=as.double(TK), x=as.double(x), phaseid = as.integer(0))
  if (ret$phaseid == 0) return("V")
  if (ret$phaseid == 1) return("VH")
  if (ret$phaseid == 10) return("L")
  if (ret$phaseid == 11) return("LH")
  if (ret$phaseid == 20) return("VL")
  if (ret$phaseid == 21) return("VLH")
  if (ret$phaseid == 30) return("F")
  if (ret$phaseid == 31) return("FH")
  return(ret$phaseid)
}

#' Saturated liquid composition at V+L surface, Function of Temperature and Pressure
#'
#' The function returns the saturated liquid composition at V+L surface
#' for given temperature and pressure.
#'
#' @param TK Temperature [K]
#' @param p Pressure [Pa]
#'
#' @return Mole fraction [-]
#'
#' @export
driesner07_H2O_NaCl_VL_xl_Tp <- function(TK, p)
{
  ret <- .C("R_driesner07_H2O_NaCl_VL_xl_pT", p=as.double(p), TK=as.double(TK), x=as.double(1.0))
  return(ret$x)
}

#' Saturated vapor composition at V+L surface, Function of Temperature and Pressure
#'
#' The function returns the saturated vapor composition at V+L surface
#' for given temperature and pressure.
#'
#' @param TK Temperature [K]
#' @param p Pressure [Pa]
#'
#' @return Mole fraction [-]
#'
#' @export
driesner07_H2O_NaCl_VL_xv_Tp <- function(TK, p)
{
  ret <- .C("R_driesner07_H2O_NaCl_VL_xv_pT", p=as.double(p), TK=as.double(TK), x=as.double(1.0))
  return(ret$x)
}

#' Saturated liquid density at V+L surface, Function of Temperature and Pressure
#'
#' The function returns the saturated liquid density at V+L surface
#' for given temperature and pressure.
#'
#' @param TK Temperature [K]
#' @param p Pressure [Pa]
#'
#' @return Density [kg/m3]
#'
#' @export
driesner07_H2O_NaCl_VL_rhol_Tp <- function(TK, p)
{
  ret <- .C("R_driesner07_H2O_NaCl_VL_rhol_pT", p=as.double(p), TK=as.double(TK), out=as.double(1.0))
  return(ret$out)
}

#' Saturated vapor density at V+L surface, Function of Temperature and Pressure
#'
#' The function returns the saturated vapor density at V+L surface
#' for given temperature and pressure.
#'
#' @param TK Temperature [K]
#' @param p Pressure [Pa]
#'
#' @return Density [kg/m3]
#'
#' @export
driesner07_H2O_NaCl_VL_rhov_Tp <- function(TK, p)
{
  ret <- .C("R_driesner07_H2O_NaCl_VL_rhov_pT", p=as.double(p), TK=as.double(TK), out=as.double(1.0))
  return(ret$out)
}

#' Liquid mass fraction in V+L phase
#'
#' The function returns the liquid mass fraction in V+L phase
#' for given temperature, pressure, and composition.
#'
#' @param TK Temperature [K]
#' @param p Pressure [Pa]
#' @param x NaCl composition [mol/mol]
#'
#' @return Mass fraction [kg/kg]
#'
#' @export
driesner07_H2O_NaCl_VL_mass_frac_l <- function(TK, p, x)
{
  xl = driesner07_H2O_NaCl_VL_xl_Tp(TK, p)
  xv = driesner07_H2O_NaCl_VL_xv_Tp(TK, p)
  X = H2ONaCl_x_to_massfrac(x)
  Xl = H2ONaCl_x_to_massfrac(xl)
  Xv = H2ONaCl_x_to_massfrac(xv)
  wl = (X-Xv)/(Xl-Xv)
  return(wl)
}

#' Vapor mass fraction in V+L phase
#'
#' The function returns the vapor mass fraction in V+L phase
#' for given temperature, pressure, and composition.
#'
#' @param TK Temperature [K]
#' @param p Pressure [Pa]
#' @param x NaCl composition [mol/mol]
#'
#' @return Mass fraction [kg/kg]
#'
#' @export
driesner07_H2O_NaCl_VL_mass_frac_v <- function(TK, p, x)
{
  wv = 1.-driesner07_H2O_NaCl_VL_mass_frac_l(TK, p, x)
  return(wv)
}

#' Liquid volume fraction in V+L phase
#'
#' The function returns the liquid volume fraction in V+L phase
#' for given temperature, pressure, and composition.
#'
#' @param TK Temperature [K]
#' @param p Pressure [Pa]
#' @param x NaCl composition [mol/mol]
#'
#' @return Volume fraction [m3/m3]
#'
#' @export
driesner07_H2O_NaCl_VL_vol_frac_l <- function(TK, p, x)
{
  xl = driesner07_H2O_NaCl_VL_xl_Tp(TK, p)
  xv = driesner07_H2O_NaCl_VL_xv_Tp(TK, p)
  X = H2ONaCl_x_to_massfrac(x)
  Xl = H2ONaCl_x_to_massfrac(xl)
  Xv = H2ONaCl_x_to_massfrac(xv)
  rhol = driesner07_H2O_NaCl_VL_rhol_Tp(TK, p)
  rhov = driesner07_H2O_NaCl_VL_rhov_Tp(TK, p)
  satl = rhov*(Xv-X)/(rhov*(Xv-X)+rhol*(X-Xl))
  return(satl)
}

#' Vapor volume fraction in V+L phase
#'
#' The function returns the vapor volume fraction in V+L phase
#' for given temperature, pressure, and composition.
#'
#' @param TK Temperature [K]
#' @param p Pressure [Pa]
#' @param x NaCl composition [mol/mol]
#'
#' @return Volume fraction [m3/m3]
#'
#' @export
driesner07_H2O_NaCl_VL_vol_frac_v <- function(TK, p, x)
{
  satv = 1.-driesner07_H2O_NaCl_VL_vol_frac_l(TK, p, x)
  return(satv)
}


#' Saturation pressure at V+L+H surface, Function of Temperature
#'
#' The function returns the saturation pressure at V+L+Halite coexisting surface
#' for given temperature.
#'
#' @param TK Temperature [K]
#' @return Pressure [Pa]
#'
#' @export
driesner07_H2O_NaCl_VLH_p_T <- function(TK)
{
  ret <- .C("R_driesner07_H2O_NaCl_VLH_p_T", TK=as.double(TK), out=as.double(1.0))
  return(ret$out)
}


#' Saturated liquid composition at V+L+H surface, Function of Temperature
#'
#' The function returns the saturated liquid composition at V+L+Halite coexisting surface
#' for given temperature.
#'
#' @param TK Temperature [K]
#' @return Mole fraction [-]
#'
#' @export
driesner07_H2O_NaCl_VLH_xl_T <- function(TK)
{
  ret <- .C("R_driesner07_H2O_NaCl_VLH_xl_T", TK=as.double(TK), out=as.double(1.0))
  return(ret$out)
}

#' Saturated vapor composition at V+L+H surface, Function of Temperature
#'
#' The function returns the saturated vapor composition at V+L+Halite coexisting surface
#' for given temperature.
#'
#' @param TK Temperature [K]
#' @return Mole fraction [-]
#'
#' @export
driesner07_H2O_NaCl_VLH_xv_T <- function(TK)
{
  ret <- .C("R_driesner07_H2O_NaCl_VLH_xv_T", TK=as.double(TK), out=as.double(1.0))
  return(ret$out)
}


#' Saturated vapor composition at V+H coexisting, Function of Temperature and Pressure
#'
#' The function returns the saturated vapor composition at V+Halite coexisting surface
#' for given temperature and pressure.
#'
#' @param TK Temperature [K]
#' @param p Pressure [Pa]
#' @return Mole fraction [-]
#'
#' @export
driesner07_H2O_NaCl_VH_xv_Tp <- function(TK, p)
{
  ret <- .C("R_driesner07_H2O_NaCl_VH_xv_pT", p=as.double(p), TK=as.double(TK), out=as.double(1.0))
  return(ret$out)
}

#' Halite volume fraction in V+H phase
#'
#' The function returns the solid volume fraction in V+H phase
#' for given temperature, pressure, and composition.
#'
#' @param TK Temperature [K]
#' @param p Pressure [Pa]
#' @param x NaCl composition [mol/mol]
#'
#' @return Volume fraction [m3/m3]
#'
#' @export
driesner07_H2O_NaCl_VH_vol_frac_h = function(TK,p,x)
{
  xv = driesner07_H2O_NaCl_VH_xv_Tp(TK,p)
  rhov = driesner07_H2O_NaCl_rho_singlephase_pTx(p, TK,xv)
  rhoh = 2163
  X = H2ONaCl_x_to_massfrac(x)
  Xv = H2ONaCl_x_to_massfrac(xv)
  Xh = 1
  (rhov* (X - Xv) )/(rhoh*(Xh - X)+rhov*(X-Xv))
}

#' Saturated liquid composition at L+H coexisting, Function of Temperature and Pressure
#'
#' The function returns the saturated liquid composition at L+Halite coexisting surface
#' for given temperature and pressure.
#'
#' @param TK Temperature [K]
#' @param p Pressure [Pa]
#' @return Mole fraction [-]
#'
#' @export
driesner07_H2O_NaCl_LH_xl_Tp <- function(TK, p)
{
  ret <- .C("R_driesner07_H2O_NaCl_LH_xl_pT", p=as.double(p), TK=as.double(TK), out=as.double(1.0))
  return(ret$out)
}

#' Salinewater density, Function of Temperature, Pressure, and composition
#'
#' The function returns the saline water density for given temperature, pressure,
#' and composition.
#'
#' @param p Pressure [Pa]
#' @param TK Temperature [K]
#' @param x Mole fraction [-]
#' @return Density [kg/m3]
#'
#' @export
driesner07_H2O_NaCl_rho_pTx <- function(p, TK, x)
{
  ret <- .C("R_driesner07_H2O_NaCl_rho_pTx", p=as.double(p), TK=as.double(TK), x=as.double(x), out=as.double(1.0))
  return(ret$out)
}

#' Salinewater density, Function of Temperature, Pressure, and composition
#'
#' The function returns the saline water density for given temperature, pressure,
#' and composition.
#'
#' @param p Pressure [Pa]
#' @param TK Temperature [K]
#' @param x Mole fraction [-]
#' @return Density [kg/m3]
#'
#' @export
driesner07_H2O_NaCl_rho_singlephase_pTx <- function(p, TK, x)
{
  ret <- .C("R_driesner07_H2O_NaCl_rho_singlephase_pTx", p=as.double(p), TK=as.double(TK), x=as.double(x), out=as.double(1.0))
  return(ret$out)
}

#' Salinewater Specific Enthalphy, Function of Temperature, Pressure, and composition
#'
#' The function returns the saline water specific enthalpy for given temperature, pressure,
#' and composition.
#'
#' @param TK Temperature [K]
#' @param p Pressure [Pa]
#' @param x Mole fraction [-]
#' @return Specific Enthalpy [J/kg]
#'
#' @export
driesner07_H2O_NaCl_singlephase_h_Tpx <- function(TK, p, x)
{
  ret <- .C("R_driesner07_H2O_NaCl_singlephase_h_pTx", p=as.double(p), TK=as.double(TK), x=as.double(x), out=as.double(1.0))
  return(ret$out)
}


#' calculate pressure from salinewater density and salinity
#'
#' @param rho Density [kg/m^3]
#' @param TK Temperature [K]
#' @param x Mole fraction [-]
#' @param pVL V-L pressure used as lower bound [Pa]
#' @param pmax upper bound of pressure [Pa]
#' @return pressure [Pa]
#'
#' @importFrom stats uniroot
#' @export
driesner07_H2O_NaCl_singlephase_p_rhoTx = function(rho, TK, x, pVL, pmax=500e6)
{
  ff = function(p) driesner07_H2O_NaCl_rho_singlephase_pTx(p, TK, x) - rho  
  ret = uniroot(ff, lower=pVL, upper=pmax)
  ret$root
}


#' calculate halite density
#'
#' @param pPa Pressure [Pa]
#' @param TK Temperature [K]
#' @return Density [kg/m^3]
#'
#' @export
driesner07_NaCl_H_rho_pT = function(pPa, TK)
{
  ret <- .C("R_driesner07_NaCl_H_rho_pT", p=as.double(pPa), TK=as.double(TK), out=as.double(1.0))
  return(ret$out)
}

#' calculate halite specific enthalpy
#'
#' @param pPa Pressure [Pa]
#' @param TK Temperature [K]
#' @return specific enthalpy [J/kg]
#'
#' @export
driesner07_NaCl_H_h_pT = function(pPa, TK)
{
  ret <- .C("R_driesner07_NaCl_H_h_pT", p=as.double(pPa), TK=as.double(TK), out=as.double(1.0))
  return(ret$out)
}

#' calculate halite specific heat 
#'
#' @param pPa Pressure [Pa]
#' @param TK Temperature [K]
#' @return specific heat [J/kg/K]
#'
#' @export
driesner07_NaCl_H_cp_pT = function(pPa, TK)
{
  ret <- .C("R_driesner07_NaCl_H_cp_pT", p=as.double(pPa), TK=as.double(TK), out=as.double(1.0))
  return(ret$out)
}
