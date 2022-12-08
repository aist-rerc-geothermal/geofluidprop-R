

#' H2O-NaCl fluid electrical conductivity model by Sen & Goode (1992)
#'
#' The model is valid for 20-200 C.
#'
#' @param m NaCl molality [mol/kg]
#' @param TC Temperature [deg C]
#'
#' @return Electrical conductivity [S/m]
#' @export
H2ONaCl_ec_SenGoode1992 <- function(m, TC)
{
  d1 = 5.6; d2 = 0.27; d3=-1.51e-4; d4=2.36; d5=0.099; d6=0.214
  sf = (d1+d2*TC+d3*TC^2)*m - (d4+d5*TC)/(1+d6*m^0.5)*m^1.5
  #  sf = (d1+d2*TC+d3*TC^2)*M - (d4+d5*TC)/(1+d6*M^0.5)*M^1.5
  sf
}


#' H2O-NaCl fluid electrical conductivity model in Sakuma and Ichiki (2016) for high pressures
#'
#' The model is valid for 0.2-2 GPa, 673-2000K, 0.6-9.6wt\%.
#'
#' @references Sakuma, H., Ichiki, M. (2016) Electrical conductivity of NaCl-H2O fluid
#' in the crust. JGR Solid Earth.
#'
#' @param p Pressure [MPa]
#' @param T Temperature [K]
#' @param c Salinity [wt\%]
#'
#' @return Electrical conductivity [ohm^-1 m^-1]
#' @export
H2ONaCl_ec_SakumaIchiki2016_highP <- function(p, T, c)
{
  # MPa, K, wt%
  
  gamma = c(
    -2.76823e-12,
    2.86668e-11,
    -1.01120e-11,
    6.32515e-9,
    -6.35950e-8,
    2.14326e-8,
    -2.92588e-6,
    2.69121e-5,
    -9.20740e-6,
    6.52051e-9,
    -7.43514e-8,
    2.23618e-8,
    -1.47966e-5,
    1.67038e-4,
    -4.54299e-5,
    6.88977e-3,
    -7.25629e-2,
    1.89836e-2,
    -2.60077e-6,
    3.64027e-5,
    -7.50611e-6,
    6.12874e-3,
    -9.01143e-2,
    1.51621e-2,
    -3.17282,
    50.2186,
    -6.22277
  )
  beta = 1:9
  ig = 0
  for (i in 1:9)
  {
    beta[ig+1] = gamma[ig*3+1]*c^2 + gamma[ig*3+2]*c + gamma[ig*3+3]
    ig = ig + 1
  }
  
  alpha1 = beta[1]*T^2 + beta[2]*T + beta[3]
  alpha2 = beta[4]*T^2 + beta[5]*T + beta[6]
  alpha3 = beta[7]*T^2 + beta[8]*T + beta[9]
  s = alpha1*p^2 + alpha2*p + alpha3
  s
}

#' H2O-NaCl fluid electrical conductivity model in Sakuma and Ichiki (2016) for low pressures
#'
#' The model is valid for <0.2 GPa, <600K, 0.6-9.6wt\% NaCl
#'
#' @references Sakuma, H., Ichiki, M. (2016) Electrical conductivity of NaCl-H2O fluid
#' in the crust. JGR Solid Earth.
#'
#' @param p Pressure [MPa]
#' @param TK Temperature [K]
#' @param c Salinity [wt\%]
#'
#' @return Electrical conductivity [ohm^-1 m^-1]
#' @export
H2ONaCl_ec_SakumaIchiki2016_lowP <- function(p, TK, c)
{
  # MPa, K, wt%
  
  phi = c(
    -1.61994e-12,
    4.32808e-11,
    1.15235e-10,
    2.52257e-10,
    1.88235e-9,
    -5.82409e-8,
    -3.37538e-7,
    -4.53779e-7,
    -5.65158e-7,
    2.70538e-5,
    2.40270e-4,
    2.97574e-4,
    4.64690e-5,
    -6.70560e-3,
    -2.69091e-2,
    -8.37212e-2,
    2.58834e-3,
    6.92510e-1,
    -3.22923,
    8.48091
  )
  delta = 1:5
  for (i in 0:(length(delta)-1))
    delta[i+1] = phi[i*4+1]*c^3 + phi[i*4+2]*c^2 + phi[i*4+3]*c + phi[i*4+4]
  
  s = delta[1]*TK^4 + delta[2]*TK^3 + delta[3]*TK^2 + delta[4]*TK + delta[5]
  s
}


#' H2O-NaCl fluid electrical conductivity model by Sinmyo and Keppler (2017)
#'
#' @references Sinmyo, R., Keppler, H. (2017) Electrical conductivity of NaCl-bearing aqueous
#' fluids to 600C and 1GPa. Contrib Mineral Petrol 172:4
#'
#' @param pMPa Pressure [MPa]
#' @param TK Temperature [K]
#' @param c_wtp Salinity [wt\%]
#'
#' @return Electrical conductivity [ohm^-1 m^-1]
#' @export
H2ONaCl_ec_SinmyoKeppler2017 <- function(pMPa, TK, c_wtp)
{
  rho = iapws95_rho_pT(pMPa*1e6, TK)*1e-3 # g/cm3
  kappa0 = 1573. -1212.*rho + 537062./TK - 208122721./TK^2
  log_s = -1.706 -93.78/TK + 0.8075*log10(c_wtp) + 3.0781*log10(rho) + log10(kappa0)
  s = 10^log_s # S/m
  #  printf("rhow=%g, kappa = %g, logs = %g, s=%g\n", rho, kappa0, log_s, s)
  s
}


#' Molar electrical conductivity model of H2O-NaCl fluids in Watanabe et al (2021, FPE)
#'
#' @param vis viscosity [Pa s]
#' @param m NaCl molality [mol/kg-H2O]
#'
#' @return Molar conductivity [Sm^2/mol]
#' @export
H2ONaCl_molar_ec_WatanabeEtAl2021_vism <- function(vis, m)
{
  coeff = rep(0, 10)
  coeff[1] = 4.169750911984e-03
  coeff[2] = -5.082058147867e-03
  coeff[3] = 5.755884911077e-01
  coeff[4] = 1.004220435951e+00
  coeff[5] = 2.550084629367e+01
  coeff[6] = 6.049110839797e-02
  coeff[7] = 2.518605467877e+06
  coeff[8] = 4.309522921111e-01
  coeff[9] = -4.892452234519e-10
  coeff[10] = -1.753387741372e-11
  
  invvis = 1./vis
  i = 0
  A = coeff[i+1] + (coeff[i+2] - coeff[i+1])/(1+(m/coeff[i+3])^coeff[i+4])
  i = i + 4
  B = coeff[i+1] + (coeff[i+2] - coeff[i+1])/(1+(m^0.5/coeff[i+3])^coeff[i+4])
  B = 1/B*1e-6
  i = i + 4
  C = coeff[i+1] + coeff[i+2]*m
  ms = A + B*invvis + C*invvis^2
  ms
}

#' H2O-NaCl fluid electrical conductivity model in Watanabe et al (2021, FPE)
#'
#' @param TK Temperature [K]
#' @param pMPa Pressure [MPa]
#' @param m NaCl molality [mol/kg-H2O]
#'
#' @return Electrical conductivity [S/m]
#' @export
H2ONaCl_ec_WatanabeEtAl2021_Tpm <- function(TK, pMPa, m)
{
  x = H2ONaCl_b_to_x(m)
  phase = driesner07_H2O_NaCl_phase_name_Tpx(TK, pMPa*1e6, x)
  # if (phase%in%c("VL","VH")) {
  #   printf("invalid TPx condition: phase=%s at TK=%g, p=%g, m=%g \n", phase, TK, pMPa, m)
  #   return(NA)
  # }
  rho = driesner07_H2O_NaCl_rho_singlephase_pTx(pMPa*1e6, TK, x)
  vis = klyukinetal2017_viscosity(rho, TK, x)
  ms = H2ONaCl_molar_ec_WatanabeEtAl2021_vism(vis, m)
  if (ms<=0) return(1e-6)
  M = H2ONaCl_b_to_M(m, rho)
  s = ms * M *1e3
  #s = max(s, 1e-6)
  if (s<1e-6) s = 1e-6
  s
}
