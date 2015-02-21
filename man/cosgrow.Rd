\name{cosgrow}
\alias{cosgrow}
\alias{cosgrowz}
\alias{cosgrowa}
\alias{cosgrowOmegaM}
\alias{cosgrowOmegaL}
\alias{cosgrowOmegaK}
\alias{cosgrowH}
\alias{cosgrowFactor}
\alias{cosgrowRate}
\alias{cosgrowFactorApprox}
\alias{cosgrowRateApprox}
\alias{cosgrowRhoCrit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Cosmological Growth and Evolution Calculations
}
\description{
  These functions allows comoving, angular size and luminosity distances to be calculated for a given redshift, it can also return look back time. They use curvature correctly, calculated internally using the relation OmegaM+OmegaL+OmegaK=1, but by default they assume a flat Universe where only OmegaM needs to be specified.
}
\usage{
  cosgrow(z=1, H0=100, OmegaM=0.3, OmegaL=1-OmegaM)
  cosgrowz(z = 1)
  cosgrowa(z = 1)
  cosgrowH(z = 1, H0 = 100, OmegaM = 0.3, OmegaL = 1 - OmegaM)
  cosgrowOmegaM(z = 1, OmegaM = 0.3, OmegaL = 1-OmegaM)
  cosgrowOmegaL(z = 1, OmegaM = 0.3, OmegaL = 1-OmegaM)
  cosgrowOmegaK(z = 1, OmegaM = 0.3, OmegaL = 1-OmegaM)
  cosgrowFactor(z = 1, OmegaM = 0.3, OmegaL = 1-OmegaM)
  cosgrowRate(z = 1, OmegaM = 0.3, OmegaL = 1-OmegaM)
  cosgrowFactorApprox(z = 1, OmegaM = 0.3, OmegaL = 1-OmegaM)
  cosgrowRateApprox(z = 1, OmegaM = 0.3, OmegaL = 1-OmegaM)
  cosgrowRhoCrit(z = 1, H0 = 100, OmegaM = 0.3, OmegaL = 1 - OmegaM)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{z}{
    Redshift, where z must be > -1 (can be a vector)
  }
  \item{H0}{
    Hubble constant as defined at z=0 (default is H0=100 (km/s)/Mpc)
  }
  \item{OmegaM}{
    Omega matter (default is 0.3)
  }
  \item{OmegaL}{
    Omega Lambda (default is for a flat Universe with OmegaL=1-OmegaM)
  }
}
\details{
  The above functions are heavily based on the equations in Hamilton A.J.S., 2001, MNRAS 322 419 and Lahav O., et al., 1991, MNRAS, 251, 136.
  
  Negative value of z> -1 are allowed, which produces future predictions based on present day cosmology.
  
  The approximation routines are generally accurate to sub 1 percent, and since they do not involve numericla integration they are substantially faster when computing large grids of numbers, i.e. they are recommended for plots, since the accuracy is sub the line width.
}
\value{
cosgrow function returns a data.frame (even if only 1 redshift if requested) with the following parameters evaluated at the respective redshift/s:
\tabular{ll}{
z       \tab Requested redshift\cr
a       \tab Universe expansion factor, as given by a=1/(1+z)\cr
H       \tab Hubble value in units of (km/s)/Mpc\cr
OmegaM  \tab Omega matter\cr
OmegaL  \tab Omega lambda\cr
OmegaK  \tab Omega curvature\cr
Factor  \tab Exact growth factor (g, see cosgrowFactor below for details)\cr
Rate    \tab Exact growth rate (f, see cosgrowRate below for details)\cr
RhoCrit \tab critical mass density of the Universe at z, in units of Msol/Mpc^3\cr
}

  cosgrowz function returns the input redshift (only included for clarity).
  
  cosgrowa function returns the Universe expansion factor, as given by a=1/(1+z).
  
  cosgrowH function returns the value of the Hubble constant at z, in units of km/s/Mpc.
  
  cosgrowOmegaM function returns the value of Omega matter at z.
  
  cosgrowOmegaL function returns the value of Omega lambda at z.
  
  cosgrowOmegaK function returns the value of Omega curvature at z.
  
  cosgrowFactor function returns the exact value of the growth factor (typically referred to as 'g' in the astronomy literature), at z.
  
  cosgrowFactorApprox function returns the approximate value of the growth factor (typically referred to as 'g' in the astronomy literature), at z.
  
  cosgrowRate function returns the exact value of the growth rate (typically referred to as 'f' in the astronomy literature), at z.
  
  cosgrowRateApprox function returns the approximate value of the growth rate (typically referred to as 'f' in the astronomy literature), at z.
  
  cosgrowRhoCrit function returns the critical mass density of the Universe at z, in units of Msol/Mpc^3.
  
}
\references{
  Based on the equations in Hamilton A.J.S., 2001, MNRAS 322 419 and Lahav O., et al., 1991, MNRAS, 251, 136
}
\author{
  Aaron Robotham
}

\seealso{
  \code{\link{cosvol}}, \code{\link{cosmap}},  \code{\link{cosdist}}
}

\examples{
  cosgrow(0.3)
  cosgrowz(0.3)
  cosgrowa(0.3)
  cosgrowH(0.3)
  cosgrowOmegaM(0.3)
  cosgrowOmegaL(0.3)
  cosgrowOmegaK(0.3)
  sum(cosgrowOmegaM(0.3)+cosgrowOmegaL(0.3)+cosgrowOmegaK(0.3)) #Still 1.
  cosgrowFactor(0.3)
  cosgrowFactorApprox(0.3) #The approximation better than 1% for reasonable cosmologies.
  cosgrowRate(0.3)
  cosgrowRateApprox(0.3) #The approximation better than 1% for reasonable cosmologies.
  cosgrowRhoCrit(0.3)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{cosmo}
\keyword{cosmology}
\keyword{growth}% __ONLY ONE__ keyword per line

