
 
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
 
 
