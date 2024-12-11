#' Critical Intermittency Calculator
#' 
#' @description
#' Calculates the critical intermittency, below which ooids will cement together rather than continue to grow
#' 
#' @details
#' This function uses a bedload transport equation and an estimate of the cement thickness necessary to cement grains together during transport to estimate the minimum intermittency, below which ooids will be cemented together and can no longer continue to grow as individual coated grains. Field observations suggest that the critical cement thickness ratio is 2-4% of the grain diameter (cement_thickness_ratio = 0.02 to 0.04) for sand-sized marine ooids in Ambergris shoal.
#' 
#' @param D ooid diameter in m
#' @param Omega \eqn{\Omega}, the CaCO<sub>3</sub> mineral saturation state, typically with respect to aragonite or calcite
#' @param cement_thickness_ratio parameter that scales the thickness of cement needed to form the hardground; this factor is multiplied by the grain diameter, so a value of 1 means that a thickness = grain diameter is required, while a value of 0.01 means that a thickness of 1% of the grain diameter is required
#' @param k CaCO<sub>3</sub> precipitation rate constant in umol/m<sup>2</sup>/hr, this is sensitive to temperature and mineralogy and can be set with an optional helper function
#' @param n CaCO<sub>3</sub> precipitation reaction order (unitless), this is also sensitive to temperature and mineralogy and can be set with an optional helper function
#' @param rho_s density of sediment in kg/m<sup>3</sup>, this is sensitive to mineralogy
#' @param M_min molar mass of mineral (g/mol), default value is for CaCO<sub>3</sub>
#' @param rho_f density of fluid (kg/m<sup>3</sup>), this is sensitive to temperature and salinity, helper function for calculating seawater density is available in Oomega Toolbox package
#' @param nu kinematic viscosity of fluid in m<sup>2</sup>/s, this is sensitive to temperature and salinity, helper function for calculating dynamic viscosity for seawater is available in Oomega Toolbox package; kinematic viscosity can be calculated from dynamic viscosity and density
#' @param H water depth in m
#' @param L_dune characteristic dune wavelength (m)
#' @param H_dune characteristic dune height (m)
#' @param H_swave significant wave height (m)
#' @param T_pwave peak wave period (s)
#' @param L_wave wave length (m)
#' @returns critical intermittency (unitless)
#' @examples
#' f_crit1 <- f_crit_solver(500*10^-6,5)
#' f_crit2 <- mapply(f_crit_solver,D=500*10^-6,Omega = seq(from=2,to=10,by=0.5))
#' f_crit3 <- mapply(f_crit_solver,D=seq(from=400,to=500,by=10)*10^-6,Omega=5)
#' 
#' @export

f_crit_solver <- function(D,
                          Omega,
                          cement_thickness_ratio = 0.04, #parameter that scales the thickness of cement needed to form the hardground; this factor is multiplied by the grain diameter, so a value of 1 means that a thickness = grain diameter is required, while a value of 0.01 means that a thickness of 1% of the grain diameter is required
                          k = 10^1.11, #rate constant (umol/m^2/hr)
                          n = 2.26, #reaction order
                          rho_s = 2800, #density of sediment (kg/m^3)
                          M_min = 100.0869, #molar mass of mineral (g/mol)
                          rho_f = 1025, #density of fluid (kg/m^3)
                          H = 1, #water depth (m)
                          nu = 9.37*10^-7, #kinematic viscosity of fluid (m^2/s)
                          L_dune = 28, #dune wavelength (m)
                          H_dune = 0.4, #dune height (m)
                          H_swave = 0.2, #significant wave height (m)
                          T_pwave = 10, #peak wave period (s)
                          L_wave = 10 #wave length (m)
) {
  g <- 9.81 #acceleration due to gravity (m/s)
  Rouse <- 2.5 #(dimensionless)
  CSF <- 1 #1 is for spheres, 0.8 is for natural
  PS <- 6 #6 is for spheres, 3.5 is for natural
  
  R_precip <- k*(Omega-1)^n #(umol/m^2/hr)
  R_precip_vol <- R_precip*M_min*10^(-9)/rho_s #(m^3/m^2/hr)
  R_precip_lin <- R_precip_vol/3600 #(m/s)
  
  alphab <- 0.015
  R <- (rho_s - rho_f)/rho_f #submerged specific density
  
  #calculate settling velocity
  Dstar <- (R*g*D^3)/(nu^2)
  X <- log10(Dstar)
  R1 <- -3.76715+1.92944*X - 0.09815*(X^2) - 0.00575*(X^3) + 0.00056*(X^4)
  R2 <- log10(1-((1-CSF)/0.85))-(((1-CSF)^2.3)*tanh(X-4.6)) + 0.3*(0.5-CSF)*((1-CSF)^2)*(X-4.6)
  R3 <- (0.65-((CSF/2.83)*tanh(X-4.6)))^(1+((3.5-PS)/2.5))
  Wstar <- R3*10^(R2+R1)
  ws <- (R*g*nu*Wstar)^(1/3) #settling velocity (m/s)
  
  ustar <- ws/(0.41*Rouse) #(m/s)
  
  #compute flow velocity
  z0 <- 3*D/30 #roughness coefficient (m)
  dz <- (H-z0)/1000 #(m)
  z <- seq(from=z0,to=H,by=dz)
  Uf <- sum((ustar/0.41)*log(z/z0)*dz/H) #depth-averaged flow velocity (m/s)
  
  #calculate mobility parameter
  k_wave <- 2*pi/L_wave
  Uw <- pi*H_swave/(T_pwave*sinh(k_wave*H)) #peak orbital wave velocity (m/s)
  Ue <- Uf + 0.4*Uw #effective velocity (m/s)
  Beta <- Uf/(Uf + Uw)
  if (D <= 500*10^-6) {
    Ucr_c <- 0.19*D^0.1*log(12*H/(3*D*1.2)) #(m/s)
    Ucr_w <- 0.24*((R)*g)^0.66*D^0.33*T_pwave^0.33 #(m/s)
  }
  if (D > 500*10^-6) {
    Ucr_c <- 8.5*D^0.6*log(12*H/(3*D*1.2)) #(m/s)
    Ucr_w <- 0.95*((R)*g)^0.57*D^0.43*T_pwave^0.14 #(m/s)
  }
  Ucr <- Beta*Ucr_c + (1-Beta)*Ucr_w #(m/s)
  Me <- (Ue-Ucr)/((R)*g*D)^0.5 #(dimensionless)
  
  #calculate bedload transport (kg/m/s)
  qb <- alphab*rho_s*Uf*H*(D/H)^1.2*Me^1.5
  qb <- qb/rho_s #(m^3/m/s)
  
  f_crit <- (L_dune*H_dune*R_precip_lin)/(qb*D*cement_thickness_ratio)
  if (is.na(f_crit)) {
    f_crit <- 1.
  }
  if (f_crit >1) {
    f_crit <- 1.
  }
  return(f_crit)
}