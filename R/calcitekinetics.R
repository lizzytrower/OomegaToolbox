#' Calcite Kinetics Interpolation Function
#' 
#' @description
#' Calculates the precipitation rate parameters (k and n) for calcite as a function of temperature.
#' 
#' @details
#' This function interpolates between the experimental constrains of the rate constant (k) and reaction order (n) for calcite precipitation, with the choice of two references to use:
#' 
#' Burton, E. A., & Walter, L. M. (1987). Relative precipitation rates of aragonite and Mg calcite from seawater: Temperature or carbonate ion control? Geology, 15, 111–114.
#' 
#' Lopez, O., Zuddas, P., & Faivre, D. (2009). The influence of temperature and seawater composition on calcite crystal growth mechanisms and kinetics: Implications for Mg incorporation in calcite lattice. Geochimica et Cosmochimica Acta, 73(2), 337–347.
#' 
#' @param tempC temperature in degrees Celsius
#' @param datasource string with the reference used for the calculation; the two options are Burton and Walter, 1989 ("BurtonWalter"), which is valid for 5 to 37 C, or Lopez et al. (2009) ("Lopezetal"), which is valid for 5 to 55 C. The default ("BurtonWalter") is recommended if you want to compare outputs for aragonite vs. calcite mineralogy, because Lopez et al. (2009) did not conduct experiments for aragonite.
#' @returns a list with two items, k (umol/m<sup>2</sup>/hr) and n (unitless)
#' @examples
#' calc1 <- calcitekinetics(25)
#' calc2 <- mapply(calcitekinetics,c(25,26,27))
#' calc3 <- calcitekinetics(25,datasource = "Lopezetal")
#' 
#' @export

calcitekinetics <- function(tempC,
                            datasource = "BurtonWalter") {
  
  if (datasource == "BurtonWalter") {
    caln <- c(0.6, 1.9, 2.3)
    calk <- c(14, 3.9, 3.7)
    calT <- c(5, 25, 37)
    
    if (tempC > 37) {
      k <- 3.7
      n <- 2.3
    }
    
    if (tempC < 5) {
      k <- 14
      n <- 0.6
    }
    
    if (tempC >= 5 && tempC <= 37) {
      k <- pracma::interp1(calT, calk, tempC)
      n <- pracma::interp1(calT, caln, tempC)
    }
  }
  
  if (datasource == "Lopezetal") {
    caln <- c(1.55, 1.84 ,2.3, 2.55)
    callogk <- c(0.2, 0.33, 0.4, 0.51)
    calT <- c(5, 25, 40, 55)
    
    if (tempC < 5) {
      k <- 10^0.2
      n <- 1.55
    }
    
    if (tempC > 55) {
      k <- 10^0.51
      n <- 2.55
    }
    
    if (tempC >= 5 && tempC <= 55) {
      k <- 10^pracma::interp1(calT, callogk, tempC)
      n <- pracma::interp1(calT, caln, tempC)
    }
  }
  
  return(list(k = k, n = n))
}
