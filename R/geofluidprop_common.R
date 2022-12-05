#' geofluidprop
#' 
#' R package to geofluidprop
#' 
#' @docType package
#' @author Norihiro Watanabe
#' @useDynLib geofluidprop
#' @exportPattern "^[[:alpha:]]+"
#' @name geofluidprop
NULL

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
 
 