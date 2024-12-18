% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Kellerhalscorrection.R
\name{kellerhals}
\alias{kellerhals}
\title{Kellerhals et al. Thin Section Grain Diameter Correction}
\usage{
kellerhals(b50)
}
\arguments{
\item{b50}{median intermediate axis dimension as measured from thin sections, any units}
}
\value{
approximation of true median sieve diameter, in the same length units as the input measurement
}
\description{
This function applies a correction to a measurement of median intermediate axis grain diameter made from thin sections.
}
\details{
This is an application of an equation from Kellerhals et al. (1975):

Kellerhals, R., Shaw, J., & Arora, V. K. (1975). On Grain Size from Thin Sections. The Journal of Geology, 83(1), 79–96.

Due to the well-rounded and highly spherical nature of ooids, we assume that the ooids are rotationally symmetrical around the major axis (A), meaning that the intermediate (B) and minor (C) axis dimensions are equal. Therefore, the parameter k\if{html}{\out{<sub>}}2\if{html}{\out{</sub>}} = 1. A correction factor based on Fig. 5b in Kellerhals et al. (1975) is used in the calculation. For other types of sand where the assumption B=C is not true, this is probably not a good correction.
}
\examples{
D50_1 <- kellerhals(500)
D50_2 <- kellerhals(c(500,550,600))

}
