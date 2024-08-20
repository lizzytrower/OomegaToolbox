#' Seawater Dynamic Viscosity Calculator
#' 
#' @description
#' Calculates seawater dynamic viscosity as a function of temperature in Celsius and salinity in g/kg
#' 
#' @details
#' Function from J. D. Isdale, C. M. Spence, and J. S. Tudhope, 
#' Physical properties of sea water solutions: viscosity, 
#' Desalination, 10(4), 319 - 328, 1972.
#' 
#' @param tempC temperature in degrees Celsius; valid for 10 to 180 degrees C
#' @param sal practical salinity in g/kg; valid for 0 to 150 g/kg
#' @returns dynamic viscosity in Pa*s
#' @examples
#' mu_sw <- seawaterdynamicviscosity(25,35)
#' 
#' @export

seawaterdynamicviscosity <- function(tempC,sal) {
  
  A <- 1.474*10^-3 + 1.5*10^-5*tempC - 3.927*10^-8*tempC^2
  B <- 1.073*10^-5 - 8.5*10^-8*tempC + 2.23*10^-10*tempC^2
  mu_w <- exp(-10.7019 + 604.129/(139.18 + tempC))
  mu_sw <- mu_w*(1 + A*sal + B*sal^2)
  
  return(mu_sw)
}