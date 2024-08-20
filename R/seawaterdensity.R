#' Seawater Density Calculator
#' 
#' @description
#' Calculates seawater density as a function of temperature in Celsius and salinity in g/kg
#' 
#' @details
#' This is a function from F. J. Millero, and A. Poisson, 
#' International one-atmosphere equation of state of seawater. 
#' Deep-Sea Research, 28A (6), 625 – 629, 1981.
#' 
#' @param tempC temperature in degrees Celsius; valid for -2 to 40 degrees C
#' @param sal practical salinity in g/kg; valid for 0 to 42 g/kg
#' @returns density in kg/m^3
#' @examples
#' rho_sw1 <- seawaterdensity(25,35)
#' rho_sw2 <- seawaterdensity(c(25,30,35),35)
#' 
#' @export

seawaterdensity <- function(tempC,sal) {
  
  A <- 0.824496 - 4.0899*10^-3*tempC + 7.6438*10^-5*tempC^2 - 8.2467*10^-7*tempC^2 + 5.3875*10^-9*tempC^4
  B <- -5.72466*10^-3 + 1.0227*10^-4*tempC - 1.6546*10^-6*tempC^2
  C <- 4.8314*10^-4
  rho_w <- 999.842594 + 6.793952*10^-2*tempC - 9.09529*10^-3*tempC^2 + 1.001685*10^-4*tempC^3 - 1.120083*10^-6*tempC^4 + 6.536336*10^-9*tempC^5
  rho_sw <- rho_w + A*sal + B*sal^1.5 + C*sal^2
  
  return(rho_sw)
}