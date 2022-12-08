
geofluidprop_dyn_name = paste("geofluidprop", .Platform$dynlib.ext, sep="")
if (file.exists(geofluidprop_dyn_name))
{
  dyn.load(geofluidprop_dyn_name)
  source("geofluidprop.R")
} else {
  cat(sprintf("ERROR: cannot find %s. Failed to load the lib.\n", geofluidprop_dyn_name))
}


#' C-like printf
#' 
#' @param ... arguments
#'
#' @export
printf <- function(...)
  cat(sprintf(...))

#' convert temperature in degree C to Kelvin
#'
#' @param TC temperature in deg C
#' @return temperature in K
#' @export
toK = function(TC)
  TC + 273.15

#' convert temperature in Kelvin to degree C
#'
#' @param TK temperature in K
#' @return temperature in deg C
#' @export
toC = function(TK)
  TK - 273.15


#' Convert H2O-NaCl concentration unit from molar concentration to mass fraction
#'
#' Convert concentration unit of saline water from molar concentration to mass fraction
#'
#' @param M Molar concentraion [mol/L] = [M]
#' @param rho Solution density [kg/m3]
#'
#' @return Mass fraction [kg/kg]
#'
#' @export
H2ONaCl_M_to_massfrac <- function(M, rho = 998.2)
{
  # g/mol
  MW_H2O = 18
  MW_NaCl = 58.44
  
  # mol/L = n_NaCl/(V_solution)
  # wt = wt_NaCl / wt_total
  
  n_NaCl = M
  X_total = 1e3*rho*1e-3 #g/L
  X = n_NaCl * MW_NaCl / X_total
  X
}

#' Convert H2O-NaCl concentration unit from molar concentration to mass fraction
#'
#' Convert concentration unit of saline water from molar concentration to mass fraction
#'
#' @param M Molar concentraion [mol/L] = [M]
#' @param rho Solution density [kg/m3]
#'
#' @return molality [mol/kg-H2O]
#'
#' @export
H2ONaCl_M_to_b <- function(M, rho = 998.2)
{
  # kg/mol
  MW_NaCl = 58.44e-3
  
  # assum 1L of solution
  n = M
  vol = 1e-3 #m3
  mass_solution = rho*vol
  mass_solute = n * MW_NaCl
  
  m = n/(mass_solution - mass_solute)
  m
}

#' Convert H2O-NaCl concentration unit from molar concentration to molar fraction
#'
#' @param M Molarity [mol/L]
#' @param rho Solution density [kg/m3]
#'
#' @return Mole fraction of NaCl [mol/mol]
#'
#' @export
H2ONaCl_M_to_x <- function(M, rho = 998.2)
{
  # kg/mol
  MW_NaCl = 58.44e-3
  MW_H2O = 18e-3
  
  # assum 1L of solution
  n_NaCl = M
  vol = 1e-3 #m3
  mass_solution = rho*vol
  mass_solute = n_NaCl * MW_NaCl
  mass_H2O = mass_solution - mass_solute
  n_H2O = mass_H2O / MW_H2O
  
  x = n_NaCl/(n_NaCl + n_H2O)
  x
}

#' Convert H2O-NaCl concentration unit from mass fraction to molar concentration
#'
#' @param X Mass fraction [kg/kg]
#' @param rho Solution density [kg/m3]
#'
#' @return Molar concentraion [mol/L]
#'
#' @export
H2ONaCl_massfrac_to_M <- function(X, rho = 998.2)
{
  # g
  wt_NaCl = X
  wt_total = 1.0 # g
  
  # cm3
  v = wt_total/(rho*1e-3)
  
  #mol
  MW_NaCl = 58.44
  m_NaCl = wt_NaCl / MW_NaCl
  
  # mol/L
  molL = m_NaCl / v * 1e3
  molL
}



#' Convert H2O-NaCl concentration unit from mass fraction to mole fraction
#'
#' Convert concentration unit of saline water
#'
#' @param mass_frac Mass fraction [-]
#'
#' @return Mole fraction [-]
#'
#' @export
H2ONaCl_massfrac_to_x <- function(mass_frac)
{
  mNaCl = mass_frac / 58.443e-3;
  mH2O = (1 - mass_frac) / 18.015e-3;
  x = mNaCl / (mNaCl + mH2O);
  x;
}


#' Convert H2O-NaCl concentration unit from mass fraction to molality
#'
#' @param X Mass fraction [kg/kg]
#'
#' @return molality [mol/kg-H2O]
#'
#' @export
H2ONaCl_massfrac_to_b <- function(X)
{
  MW_NaCl = 58.44e-3
  
  # mass for 1kg solution
  mass_NaCl = X
  mass_H2O = 1.-mass_NaCl
  
  # mol
  n_NaCl = mass_NaCl / MW_NaCl
  
  # mol/kg-solvent
  m = n_NaCl / mass_H2O
  m
}

#' Convert H2O-NaCl concentration unit from molality to mole fraction
#'
#' Convert concentration unit of saline water
#'
#' @param molality molality [mol/kg-H2O]
#'
#' @return Mole fraction [-]
#'
#' @export
H2ONaCl_b_to_x <- function(molality)
{
  mNaCl = molality*1.0; # 1kg water
  mH2O = 1./18.015e-3; # mole of water
  x = mNaCl / (mNaCl + mH2O);
  x;
}

#' Convert H2O-NaCl concentration unit from molality to mole fraction
#'
#' Convert concentration unit of saline water
#'
#' @param m molality [mol/kg-H2O]
#' @param rho density [kg/m3]
#'
#' @return Molar concentration [mol/L]
#'
#' @export
H2ONaCl_b_to_M <- function(m, rho=998.2)
{
  # kg/mol
  MW_NaCl = 58.44e-3
  
  # M = 1e3*m*rho / (m*MW_NaCl+1000)
  
  # assum 1kg of solvent
  mass_solvent = 1
  n = m * mass_solvent
  mass_solute = n * MW_NaCl
  vol_solution = (mass_solvent + mass_solute) / rho
  
  M = n / (vol_solution*1e3)
  M
}

#' Convert H2O-NaCl concentration unit from molality to mass fraction
#'
#' @param b molality [mol/kg-H2O]
#'
#' @return mass fraction of NaCl [kg/kg]
#'
#' @export
H2ONaCl_b_to_massfrac <- function(b)
{
  MW_NaCl = 58.44e-3
  
  # assume 1kg solvent
  n_NaCl = b
  m_NaCl = n_NaCl * MW_NaCl
  
  # kg/kg
  X = m_NaCl / (m_NaCl + 1)
  X
}

#' Convert H2O-NaCl concentration unit from mole fraction to mass fraction
#'
#' Convert concentration unit of saline water
#'
#' @param mole_frac Mole fraction [-]
#'
#' @return Mass fraction [-]
#'
#' @export
H2ONaCl_x_to_massfrac <- function(mole_frac)
{
  W_NaCl = mole_frac * 58.443e-3
  W_H2O = (1 - mole_frac) * 18.015e-3
  X = W_NaCl / (W_NaCl + W_H2O)
  X
}




#' Phase name of H2O-NaCl fluid as a function of temperature, pressure, and NaCl composition
#'
#' @param TK Temeprature [K]
#' @param p Pressure [Pa]
#' @param x Bulk composition (mole fraction) of NaCl [mol/mol]
#' @return phase name
#'
#' @export
driesner07_H2O_NaCl_phase_name_Tpx <- function(TK,p,x)
{
  phaseid = driesner07_H2O_NaCl_phase_type(p, TK, x)
  if (phaseid == 0) return("V")
  if (phaseid == 1) return("VH")
  if (phaseid == 10) return("L")
  if (phaseid == 11) return("LH")
  if (phaseid == 20) return("VL")
  if (phaseid == 21) return("VLH")
  if (phaseid == 30) return("F")
  if (phaseid == 31) return("FH")
  return(phaseid)
}

#' Estimate density from concentration
#'
#' @param M0 Molar concentration at room temperature [mol/L]
#'
#' @return density [kg/m3]
#'
#' @export
driesner07_H2O_NaCl_estimate_density_from_M0 <- function(M0)
{
  TK = toK(20)
  p = 0.101325e6
  MW_H2O = 18e-3 #kg/mol
  MW_NaCl = 58.44e-3
  vol = 1e-3 #m3
  nNaCl = M0 # in 1L of solution
  
  nH2O = 1./MW_H2O
  converged = FALSE
  for (i in 1:50)
  {
    x = nNaCl / (nH2O + nNaCl)
    
    rho = driesner07_H2O_NaCl_rho_singlephase_pTx(TK, p, x)
    mass_solution = rho * vol #kg
    #total mass = nH2O * M_H2O  + nNaCl * M_NaCl
    nH2O_new = (mass_solution - nNaCl * MW_NaCl) / MW_H2O
    error = nH2O_new - nH2O
    #printf("%d: nH2O = %g, error = %g\n", i, nH2O_new, error/nH2O_new)
    if (abs(error) < 1e-3*nH2O_new) {
      converged = TRUE
      break
    }
    nH2O = nH2O_new
  }
  
  if (!converged) {
    print("*** error: diverged in H2ONaCl_estimate_density_from_M0()")
    return(NA)
  }
  
  rho
}


