#' Oomega Solver
#' 
#' @description
#' Calculates the CaCO<sub>3</sub> saturation state (\eqn{\Omega}) value at which precipitation and abrasion rates are balanced for the given ooid diameter.
#' 
#' @details
#' This function is an R implementation of the OomegaSolver function described in the following references:
#' 
#' Ingalls, M., Fetrow, A. C., Snell, K. E., Frantz, C. M., & Trower, E. J. (2022). Lake level controls the recurrence of giant stromatolite facies. Sedimentology, 69(4), 1649–1674.
#' 
#' Trower, E. J., Smith, B. P., Koeshidayatullah, A. I., & Payne, J. L. (2022). Marine ooid sizes record Phanerozoic seawater carbonate chemistry. Geophysical Research Letters, 49(22), e2022GL100800.
#' 
#' The basic premise is to use measurements of ooid diameter (which may need to be corrected if measured from thin sections) to reconstruct the carbonate chemistry of the waters they formed in.
#' 
#' @param D ooid diameter in m
#' @param k CaCO<sub>3</sub> precipitation rate constant in umol/m<sup>2</sup>/hr, this is sensitive to temperature and mineralogy and can be set with an optional helper function
#' @param n CaCO<sub>3</sub> precipitation reaction order (unitless), this is also sensitive to temperature and mineralogy and can be set with an optional helper function
#' @param rho_s density of sediment in kg/m<sup>3</sup>, this is sensitive to mineralogy
#' @param M_min molar mass of mineral (g/mol), default value is for CaCO<sub>3</sub>
#' @param rho_f density of fluid (kg/m<sup>3</sup>), this is sensitive to temperature and salinity, helper function for calculating seawater density is available in Oomega Toolbox package
#' @param nu kinematic viscosity of fluid in m<sup>2</sup>/s, this is sensitive to temperature and salinity, helper function for calculating dynamic viscosity for seawater is available in Oomega Toolbox package; kinematic viscosity can be calculated from dynamic viscosity and density
#' @param H water depth in m
#' @param intermittency intermittency of movement, must be in the range 0-1 where 0 implies no movement and 1 implies constant movement
#' @returns calculated \eqn{\Omega} value
#' @examples
#' Omega <- OomegaSolver(D = 500*10^-6,
#'                       H = 2,
#'                       intermittency = 0.1)
#'                       
#' Omega <- mapply(OomegaSolver, 
#'                 D = c(400,500,600)*10^-6, 
#'                 H = 1.5)
#' 
#' @export

OomegaSolver <- function(
    D, #grain diameter (m)
    k = 10^1.11, #rate constant (umol/m^2/hr)
    n = 2.26, #reaction order
    rho_s = 2800, #density of sediment (kg/m^3)
    M_min = 100.0869, #molar mass of mineral (g/mol)
    rho_f = 1025, #density of fluid (kg/m^3)
    nu = 9.37*10^-7, #kinematic viscosity of fluid (m^2/s)
    H = 1, #water depth (m)
    intermittency = 0.15) {
  
  R <- (rho_s - rho_f)/rho_f #submerged specific density
  A1 <- 0.36 #(dimensionless)
  kv <- 9*10^5 #(dimensionless)
  young <- 20*10^9 #Young's modulus (Pa)
  strength <- 1*10^6 #tensile strength (Pa)
  Rouse <- 2.5 #(dimensionless)
  Stc <- 9 #Stokes threshold (dimensionless)
  CSF <- 1 #Corey Shape Factor, 1 is for spheres, 0.8 is for natural
  PS <- 6 #Powers roundness, 6 is for spheres, 3.5 is for natural
  g <- 9.81 #acceleration due to gravity (m/s^2)
  
  #set critical Shields number (following Li et al., 2021):
  if (D <= 2*10^-3){
    tau_c <- 0.03
  }
  if (D > 2*10^3 && D <= 4*10^-3) {
    tau_c <- 0.04
  }
  if (D > 4*10^-3 && D <= 8*10^-3) {
    tau_c <- 0.043
  }
  if (D > 8*10^-3) {
    tau_c <- 0.045
  }
  
  #calculate settling velocity
  Dstar <- (R*g*D^3)/(nu^2)
  X <- log10(Dstar)
  R1 <- -3.76715+1.92944*X - 0.09815*(X^2) - 0.00575*(X^3) + 0.00056*(X^4)
  R2 <- log10(1-((1-CSF)/0.85))-(((1-CSF)^2.3)*tanh(X-4.6)) + 0.3*(0.5-CSF)*((1-CSF)^2)*(X-4.6)
  R3 <- (0.65-((CSF/2.83)*tanh(X-4.6)))^(1+((3.5-PS)/2.5))
  Wstar <- R3*10^(R2+R1)
  ws <- (R*g*nu*Wstar)^(1/3) #settling velocity (m/s)
  
  cdrag <- (4/3)*(R*g*D)/(ws^2) #coefficient of drag (units)
  ustar <- ws/(0.41*Rouse) #(m/s)
  tau <- ustar^2/(R*g*D)
  tstage <- tau/tau_c
  
  #compute flow velocity
  z0 <- 3*D/30 #roughness coefficient (m)
  dz <- (H-z0)/1000 #(m)
  z <- seq(from=z0,to=H,by=dz)
  Uf <- sum((ustar/0.41)*log(z/z0)*dz/H) #depth-averaged flow velocity (m/s)
  
  #compute bed load height and velocity
  hb <- D*1.44*(tstage-1)^0.5 #height of bed load layer (m)
  Us <- (R*g*D)^0.5*1.56*(tstage-1)^0.56 #bed load velocity (m/s)
  
  if (Us > Uf) {
    Us <- Uf
  }
  
  #compute suspended load
  if (hb < H) {
    hb[hb < D] <- D
    b <- hb
    res <- 1000
    di5 <- (log(H)-log(b))/res
    i5 <- seq(from=log(b),to=log(H),by=di5)
    z <- exp(i5)
    z[length(z)] <- H
    dz <- diff(z)
    dz <- c(dz[1],dz)
    a1 <- sum((((1-(z[z>z0]/H))/(1-(b/H)))*(b/z[z>z0]))^Rouse*log(z[z>z0]/z0))/(Uf*H)*(ustar/0.41)
    cb <- 1/(Us*hb+Uf*H*a1)
    
    #find concentration profile
    c <- 0
    c[1] <- cb
    c[2:(length(z)+1)] <- cb*(((1-(z/H))/(1-(b/H)))*(b/z))^Rouse
    z <- c(0,z)
    c[z==H] <- 0
    
    #calculate fall distance
    gradc <- rep(0,length(c))
    gradc[2:length(c)] <- -diff(c)
    Hfall <- (1/cb)*sum(z*gradc)
  }
  
  else {
    hb <- H
    cb <- 1/(Us*hb)
    Hfall <- hb
    a1 <- 0
  }
  
  if (cb == 0) {
    Hfall <- 0
  }
  
  sig <- ustar
  dx <- sig/100
  X <- seq(from=-6*sig,to=6*sig,by=dx)
  f <- dnorm(X,mean=0,sd=sig)
  X <- X/ws
  
  Scos <- 1
  
  wfall <- Scos*((2*(2/3)*D*g/cdrag*R)*(1-exp(-cdrag*rho_f/rho_s*(Hfall/Scos)/(2/3*D))))^0.5
  
  psifall <- wfall/ws
  settlematrix <- psifall + X
  settlematrix1 <- settlematrix
  settlematrix[settlematrix<0] <- 0
  psifall_turb <- sum((settlematrix)*f)*dx
  psi_fall3 <- sum((settlematrix^3)*f)*dx
  E1 <- psi_fall3
  
  wi_st <- settlematrix
  wi_st[(D*wi_st*ws*rho_s/(9*nu*rho_f))<Stc] <- 0
  psi_fall3_st <- sum((wi_st^3)*f)*dx
  E1_st <- psi_fall3_st
  
  ti <- D/Hfall
  ti[Hfall<= 0.5*D] <- 0
  
  En_suspt_st <- 1/kv*(ws/(g*D)^0.5)^3*E1_st*ti/6
  En_suspt_st[En_suspt_st<0] <- 0
  
  Efactor <- 60*60*rho_s*young*(g*D)^(3/2)/(strength)^2
  
  if (tstage<=1) {
    En_suspt_st <- 0
  }
  
  if (H<D) {
    En_suspt_st <- 0
  }
  
  R_abrasion_val <- En_suspt_st*Efactor*A1*4*pi*(D/2)^2
  
  SSA_ooid <- pi*23*D^2 #(m^2)
  
  Oomega <- (R_abrasion_val*intermittency*rho_s/(k*M_min*10^-9*SSA_ooid))^(1/n)+1
  
  return(Oomega)
}