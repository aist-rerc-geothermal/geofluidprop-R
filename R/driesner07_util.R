

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
  phaseid = driesner07_H2O_NaCl_phase_Tpx(TK, p, x)
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

#' search possibly multiple intervals of an input argument where
#' thee returned value of a given function has a different sign
#'
search_intervals = function(f, p0, p1, n = 100)
{
  dat = data.frame(p = seq(p0, p1, length.out = n))
  dat$f = mapply(f, dat$p)
  dat = dplyr::filter(dat, !is.na(f))

  sign0 = dat$f[1]
  if (is.na(sign0)) {
    #printf("is.na(sign0) at p0=%g, p1=%g\n", p0*1e-6, p1*1e-6)
    return(c())
  }
  p_intervals = c(p0)
  for (i in 2:nrow(dat))
  {
    if (is.na(dat$f[i])) {
      printf("is.na(dat$f[i]) at p0=%g, p1=%g\n", p0*1e-6, p1*1e-6)
    }
    if (sign0 * dat$f[i] < 0) {
      sign0 = dat$f[i]
      p_intervals = c(p_intervals, dat$p[i])
    }
  }
  p_intervals
}

#' get vapor phase pressure(s) in V+L region at the given T and x
#'
#' @param TK Temperature [K]
#' @param x Bulk composition (mole fraction) of NaCl [mol/mol]
#'
#' @return a list of vapor pressures [Pa]
#'
#' @export
driesner07_H2O_NaCl_VL_pv_Tx = function(TK, x)
{
  xc = driesner07_H2O_NaCl_xc_T(TK)
  pc = driesner07_H2O_NaCl_pc_T2(TK)
  p_VLH = driesner07_H2O_NaCl_VLH_p_T(TK)
  if (abs(x - xc)<1e-8*xc) return(c(pc, NA))

  f = function(p) (driesner07_H2O_NaCl_VL_xv_Tp(TK, p) - x)
  p_intervals = search_intervals(f, p_VLH, pc*0.999999)
  if (length(p_intervals) <= 1)
  {
    #printf("no interval found\n")
    return(c(NA,NA))
  } else if (length(p_intervals) == 2) {
    ret = uniroot(f, c(p_intervals[1], p_intervals[2]))
    return( c(ret$root, NA) )
  } else if (length(p_intervals) == 3) {
    ret = uniroot(f, c(p_intervals[1], p_intervals[2]))
    pv1 = ret$root
    ret = uniroot(f, c(p_intervals[2], p_intervals[3]))
    pv2 = ret$root
    return( c(pv1, pv2) )
  }
}

#' get liquid phase pressure in V+L region at the given T and x
#'
#' @param TK Temperature [K]
#' @param x Bulk composition (mole fraction) of NaCl [mol/mol]
#'
#' @return a list of liquid pressures [Pa]
#'
#' @export
driesner07_H2O_NaCl_VL_pl_Tx = function(TK, x)
{
  xc = driesner07_H2O_NaCl_xc_T(TK)
  pc = driesner07_H2O_NaCl_pc_T2(TK)
  p_VLH = driesner07_H2O_NaCl_VLH_p_T(TK)
  if (x == xc) return(pc)

  # printf("TC=%g, x=%g, xc=%g\n", toC(TK), x, xc)

  f = function(p) (driesner07_H2O_NaCl_VL_xl_Tp(TK, p) - x)
  p_intervals = search_intervals(f, p_VLH, pc)
  if (length(p_intervals) <= 1)
  {
    #printf("no interval found\n")
    return(NA)
  } else if (length(p_intervals) == 2) {
    ret = uniroot(f, c(p_intervals[1], p_intervals[2]))
    return( c(ret$root) )
  } else if (length(p_intervals) == 3) {
    ret = uniroot(f, c(p_intervals[1], p_intervals[2]))
    pv1 = ret$root
    ret = uniroot(f, c(p_intervals[2], p_intervals[3]))
    pv2 = ret$root
    return( c(pv1, pv2) )
  }
}

#' get temperature of VLH curve at the given vapor x
#'
#' @param xv NaCl composition (mole fraction) of vapor [mol/mol]
#'
#' @return Temperature [K]
#'
#' @export
driesner07_H2O_NaCl_VLH_T_xv = function(xv)
{
  # if (xv>driesner07_H2O_NaCl_wt_to_x(0.1e-2))
  #   return(NA)

  f = function(TK) (driesner07_H2O_NaCl_VLH_xv_T(TK) - xv)

  intervals = search_intervals(f, toK(10), toK(750))
  if (length(intervals)<=1)
    return(c(NA))

  if (length(intervals)==2)
  {
    TK0 = intervals[1]
    TK1 = intervals[2]
    ret = try(uniroot(f, c(TK0, TK1)), silent = TRUE)
    if (class(ret) != "try-error")
      return( c(ret$root) )
  }
  else if (length(intervals)==3)
  {
    TK0 = intervals[1]
    TK1 = intervals[2]
    ret = try(uniroot(f, c(TK0, TK1)), silent = TRUE)
    if (class(ret) != "try-error")
      out = c(ret$root)
    TK0 = intervals[2]
    TK1 = intervals[3]
    ret = try(uniroot(f, c(TK0, TK1)), silent = TRUE)
    if (class(ret) != "try-error")
      out = c(out, ret$root)
    return(out)
  } else {
    print("nr. intervals > 3 is not supported in driesner07_H2O_NaCl_VLH_T_xv()")
  }

  NA
}

#' get an array of T and p values along a phase boundary between F and VL regions
#'
#' @param massfrac_NaCl NaCl mass fraction (bulk) [kg/kg]
#' @param TCmax Maximum temperature to be calculated [K]
#'
#' @return dataframe object including calculated TK and p values
#'
#' @importFrom dplyr select filter
#' @importFrom rlang .data
#' @export
driesner07_H2O_NaCl_get_Tp_curve_on_F_VL_boundary = function(massfrac_NaCl, TCmax)
{
  # phase
  x = H2ONaCl_massfrac_to_x(massfrac_NaCl)
  Tc = driesner07_H2O_NaCl_Tc_x(x)
  pc = driesner07_H2O_NaCl_pc_T2(Tc)
  TK_VLH = driesner07_H2O_NaCl_VLH_T_xv(x)[1]

  # generate a boiling curve
  #printf("-> boiling curve\n")
  dat = data.frame(TK=seq(toK(10), Tc-0.01, length.out = 20))
  dat$p = mapply(driesner07_H2O_NaCl_VL_pl_Tx, dat$TK, x)
  dall = dat

  # add critical point
  dat = data.frame(TK=c(Tc), p=c(pc))
  dall = rbind(dall, dat)

  # generate a condensation curve
  #printf("-> condensation curve\n")
  dat = data.frame(TK=seq(Tc+0.01, toK(TCmax), length.out = 100))
  dat$xv_VLH = mapply(driesner07_H2O_NaCl_VLH_xv_T, dat$TK)
  dat = dplyr::filter(dat, .data$xv_VLH <= x)
  if (nrow(dat)>0)
  {
    if (!is.na(TK_VLH))
      TT = seq(350, toC(TK_VLH), 0.2)
    else
      TT = seq(350, TCmax, 1)
    #    TT = seq(350, max(dat$TC), length.out = 80)
    #    TT = c(seq(374, 374.9, 0.1), TT)
    dat = data.frame(TC=TT)
    #dat = data.frame(TC=seq(min(dat$TC), max(dat$TC), length.out = 100))
    dat$TK = toK(dat$TC)

    # draw lower curve
    dat$p = mapply(function(TK,x)driesner07_H2O_NaCl_VL_pv_Tx(TK,x)[1], dat$TK, x)
    dat2 = dplyr::filter(dat, !is.na(.data$p))
    if (nrow(dat2)>0) {
      dall = rbind(dall, dplyr::select(dat2, .data$TK, .data$p))
    }

    # draw upper curve
    if (nrow(dat)>0)
    {
      dat2 = dat
      dat2$p = mapply(function(TK,x)driesner07_H2O_NaCl_VL_pv_Tx(TK,x)[2], dat$TK, x)
      dat2 = dplyr::filter(dat2, !is.na(.data$p))
      if (nrow(dat2)>0)
      {
        dall = rbind(dall, dplyr::select(dat2, .data$TK, .data$p))

        # fill gap
        dat_min = dplyr::filter(dat, .data$TK==min(dat$TK))
        dat2_min = dplyr::filter(dat2, .data$TK==min(dat2$TK))
        dat3 = rbind(dat_min, dat2_min)
      }
    }
  }

  dall
}

#' get an array of T and p values along a phase boundary between VL and VH regions
#'
#' @param TCmax Maximum temperature to be calculated [K]
#' @param n the number of sampling temperatures
#'
#' @return dataframe object including calculated TK and p values
#'
#' @importFrom dplyr select
#' @importFrom rlang .data
#' @export
driesner07_H2O_NaCl_get_Tp_curve_on_VL_VH_boundary = function(TCmax, n=50)
{
  dat = data.frame(TC=seq(10, TCmax, length.out = n))
  dat$TK = toK(dat$TC)
  dat$p = mapply(driesner07_H2O_NaCl_VLH_p_T, dat$TK)
  dat$xv = mapply(driesner07_H2O_NaCl_VLH_xv_T, dat$TK)
  dat$wtp = mapply(H2ONaCl_x_to_massfrac, dat$xv) * 1e2
  dplyr::select(dat, .data$TK, .data$p)
}

#' get information about H2O-NaCl phase relation at the given NaCl mass fraction, based on Driesner & Heinrich (2007)
#'
#' @param massfrac_NaCl NaCl mass fraction (bulk) [kg/kg]
#' @param TC_max Maximum temperature to be calculated [K]
#'
#' @return a list object including all the phase relation information
#'
#' @export
driesner07_H2O_NaCl_get_phase_relation_on_Tp_space = function(massfrac_NaCl, TC_max = 800)
{
  x = H2ONaCl_massfrac_to_x(massfrac_NaCl)

  dat = list()
  dat$massfrac_NaCl = massfrac_NaCl

  # critical point
  Tc = driesner07_H2O_NaCl_Tc_x(x)
  pc = driesner07_H2O_NaCl_pc_T2(Tc)
  dat$cp = data.frame(TK = Tc, p = pc)

  # Boundary curve between LV and VH regions
  dat$boundary_LV_VH = driesner07_H2O_NaCl_get_Tp_curve_on_VL_VH_boundary(TC_max)

  # Boundary curve between F and LV regions
  dat$boundary_F_LV = driesner07_H2O_NaCl_get_Tp_curve_on_F_VL_boundary(massfrac_NaCl, TC_max)

  dat
}

#' plot H2O-NaCl phase diagram on T-p space
#'
#' @param dat H2O-NaCl phase relation information
#' @param TC_range_max Maximum temperature to be plotted [K]
#'
#' @importFrom ggplot2 ggplot aes geom_line geom_point coord_cartesian xlab ylab ggtitle
#' @export
driesner07_H2O_NaCl_plot_phase_relation_on_Tp_space = function(dat, TC_range_max = 800)
{
  g = ggplot2::ggplot()
  g = g + ggplot2::geom_line(data=dat$boundary_LV_VH, ggplot2::aes(x=toC(TK), y=p*1e-6))
  #g = g + geom_point(data=dat$boundary_LV_VH, aes(x=toC(TK), y=p*1e-6), size=1)
  g = g + ggplot2::geom_line(data=dat$boundary_F_LV, ggplot2::aes(x=toC(TK), y=p*1e-6), color="blue")
  #g = g + geom_point(data=dat$boundary_F_LV, aes(x=toC(TK), y=p*1e-6), size=1, color="blue")
  g = g + ggplot2::geom_point(data=dat$cp, ggplot2::aes(x=toC(TK), y=p*1e-6), size=3, shape=21, fill="white")
  g = g + ggplot2::coord_cartesian(xlim=c(0, TC_range_max), ylim=c(0,150))
  g = g + ggplot2::xlab("T [deg. C]") + ggplot2::ylab("p [MPa]")
  g = g + ggplot2::ggtitle(sprintf("%g wt%% NaCl", dat$massfrac_NaCl*1e2))

  print(g)
}

#' @noRd
write_lines = function(filename, lines)
{
  fileConn <- file(filename, open = "wt")
  writeLines(lines, fileConn)
  close(fileConn)
}

#' @noRd
append_lines = function(filename, lines)
{
  fileConn <- file(filename, open = "at")
  writeLines(lines, fileConn)
  close(fileConn)
}

#' @noRd
to_TC_pMPa = function(dat)
{
  dat$TC = toC(dat$TK)
  dat$pMPa = dat$p * 1e-6
  dplyr::select(dat, .data$TC, .data$pMPa)
}

#' save H2O-NaCl phase relations on T-p space to a text file
#'
#' @param dat H2O-NaCl phase relation information
#' @param filename Output file name
#'
#' @importFrom dplyr select
#' @importFrom rlang .data
#' @export
driesner07_H2O_NaCl_write_phase_relation_on_Tp_space = function(dat, filename)
{
  
  write_lines(filename, sprintf("H2O-NaCl phase relations at %g wt%%NaCl", dat$massfrac_NaCl*1e2))

  append_lines(filename, "# LV-VH boundary")
  write.table(to_TC_pMPa(dat$boundary_LV_VH), filename,
            append = TRUE, row.names = FALSE, quote = FALSE)

  append_lines(filename, "# F-LV boundary")
  write.table(to_TC_pMPa(dat$boundary_F_LV), filename,
            append = TRUE, row.names = FALSE, quote = FALSE)

  append_lines(filename, "# C.P.")
  write.table(to_TC_pMPa(dat$cp), filename,
            append = TRUE, row.names = FALSE, quote = FALSE)

}

