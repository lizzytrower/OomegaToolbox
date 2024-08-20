#' Kellerhals et al. Thin Section Grain Diameter Correction
#' 
#' @description
#' This function applies a correction to a measurement of median intermediate axis grain diameter made from thin sections.
#' 
#' @details
#' This is an application of an equation from Kellerhals et al. (1975):
#' 
#' Kellerhals, R., Shaw, J., & Arora, V. K. (1975). On Grain Size from Thin Sections. The Journal of Geology, 83(1), 79–96.
#' 
#' Due to the well-rounded and highly spherical nature of ooids, we assume that the ooids are rotationally symmetrical around the major axis (A), meaning that the intermediate (B) and minor (C) axis dimensions are equal. Therefore, the parameter k<sub>2</sub> = 1. A correction factor based on Fig. 5b in Kellerhals et al. (1975) is used in the calculation. For other types of sand where the assumption B=C is not true, this is probably not a good correction.
#' 
#' @param b50 median intermediate axis dimension as measured from thin sections, any units
#' @returns approximation of true median sieve diameter, in the same length units as the input measurement
#' @examples
#' D50_1 <- kellerhals(500)
#' D50_2 <- kellerhals(c(500,550,600))
#' 
#' @export

kellerhals <- function(b50) {
  
  k2 <- 1 #implies B = C, reasonable assumption for ooids
  Cb_corr <- -21.29 #from Fig. 5b in Kellerhals et al. for k2 = 1
  D50 <- (1 - Cb_corr/100)*b50/(2*k2)*(2*(1 + k2^2))^0.5
  
  return(D50)
}
