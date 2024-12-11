#' Maximum Intermittency Calculator
#' 
#' @description
#' Calculates the maximum intermittency, above which precipitation cannot keep up with rapid abrasion
#' 
#' @details
#' This function uses a similar framework as OomegaSolver and DeqSolver to solve for the intermittency value at a theoretical maximum \eqn{\Omega} value for a system. The idea is that this maximum Omega represents the maximum precipitation rate in that system; therefore the intermittency at such a precipitation rate represents the maximum possible intermittency because beyond this value, precipitation will always be outpaced by abrasion. Note that for many conditions, this function may return a value of 1 because the scaling of abrasion rate with size means that this is not a useful constraint for smaller ooids. It can be useful for larger ooids.
#' 
#' @param D ooid diameter in m
#' @param Omega_max \eqn{\Omega}, the CaCO<sub>3</sub> mineral saturation state, typically with respect to aragonite or calcite; here you should choose a maximum plausible value. For example, in normal seawater for aragonite precipitation, 20 is a reasonable maximum because it is the threshold of homogeneous nucleation.
#' @param k CaCO<sub>3</sub> precipitation rate constant in umol/m<sup>2</sup>/hr, this is sensitive to temperature and mineralogy and can be set with an optional helper function
#' @param n CaCO<sub>3</sub> precipitation reaction order (unitless), this is also sensitive to temperature and mineralogy and can be set with an optional helper function
#' @param rho_s density of sediment in kg/m<sup>3</sup>, this is sensitive to mineralogy
#' @param M_min molar mass of mineral (g/mol), default value is for CaCO<sub>3</sub>
#' @param rho_f density of fluid (kg/m<sup>3</sup>), this is sensitive to temperature and salinity, helper function for calculating seawater density is available in Oomega Toolbox package
#' @param nu kinematic viscosity of fluid in m<sup>2</sup>/s, this is sensitive to temperature and salinity, helper function for calculating dynamic viscosity for seawater is available in Oomega Toolbox package; kinematic viscosity can be calculated from dynamic viscosity and density
#' @param H water depth in m
#' @returns maximum intermittency (unitless)
#' @examples
#' f_max1 <- f_max_solver(2000*10^-6,20)
#' f_max2 <- mapply(f_max_solver,D=2000*10^-6,Omega_max = seq(from=20,to=25,by=0.5))
#' f_max3 <- mapply(f_max_solver,D=seq(from=1000,to=2000,by=100)*10^-6,Omega_max=20)
#' 
#' @export

f_max_solver <- function(
    D, #grain diameter (m)
    Omega_max,
    k = 10^1.11, #rate constant (umol/m^2/hr)
    n = 2.26, #reaction order
    rho_s = 2800, #density of sediment (kg/m^3)
    M_min = 100.0869, #molar mass of mineral (g/mol)
    rho_f = 1025, #density of fluid (kg/m^3)
    nu = 9.37*10^-7, #kinematic viscosity of fluid (m^2/s)
    H = 1 #water depth (m)
) {
  
  R <- (rho_s - rho_f)/rho_f #submerged specific density
  A1 <- 0.36 #(dimensionless)
  kv <- 9*10^5 #(dimensionless)
  young <- 20*10^9 #Youngs' modulus (Pa)
  strength <- 1*10^6 #tensile strength (Pa)
  g <- 9.81 #acceleration due to gravity (m/s^2)
  Rouse <- 2.5 #(dimensionless)
  Stc <- 9 #Stokes threshold (dimensionless)
  CSF <- 1 #1 is for spheres, 0.8 is for natural
  PS <- 6 #6 is for spheres, 3.5 is for natural
  
  #set critical Shields number (following Li et al., 2021):
  if (D <= 2*10^-3){
    tau_c <- 0.03
  }
  if (D > 2*10^-3 && D <= 4*10^-3) {
    tau_c <- 0.04
  }
  if (D > 4*10^-3 && D <= 8*10^-3) {
    tau_c <- 0.043
  }
  if (D > 8*10^-3) {
    tau_c <- 0.045
  }
  
  #standard carbonate precipitation rate function
  R_precip <- k*(Omega_max-1)^n #(umol/m^2/hr)
  #calculate specific surface area
  SSA_ooid <- pi*23*D^2 #(m^2)
  R_growth <- R_precip*M_min*10^(-9)*SSA_ooid/rho_s
  
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
  
  f <- R_growth/R_abrasion_val
  if (f>1) {
    f <- 1
  }
  
  return(f)
}