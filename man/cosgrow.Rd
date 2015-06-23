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
\alias{cosgrowSigma8}
\alias{cosgrowFactorApprox}
\alias{cosgrowRateApprox}
\alias{cosgrowSigma8Approx}
\alias{cosgrowRhoCrit}
\alias{cosgrowRhoMean}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Cosmological Growth and Evolution Calculations
}
\description{
  These functions allows comoving, angular size and luminosity distances to be calculated for a given redshift, it can also return look back time. They use curvature correctly, calculated internally using the relation OmegaM+OmegaL+OmegaK=1, but by default they assume a flat Universe where only OmegaM needs to be specified.
}
\usage{
  cosgrow(z=1, H0=100, OmegaM=0.3, OmegaL=1-OmegaM, Sigma8=0.8, fSigma8=FALSE, ref)
  cosgrowz(z = 1)
  cosgrowa(z = 1)
  cosgrowH(z=1, H0=100, OmegaM=0.3, OmegaL=1-OmegaM, ref)
  cosgrowOmegaM(z=1, OmegaM=0.3, OmegaL=1-OmegaM, ref)
  cosgrowOmegaL(z=1, OmegaM=0.3, OmegaL=1-OmegaM, ref)
  cosgrowOmegaK(z=1, OmegaM=0.3, OmegaL=1-OmegaM, ref)
  cosgrowFactor(z=1, OmegaM=0.3, OmegaL=1-OmegaM, ref)
  cosgrowRate(z=1, OmegaM=0.3, OmegaL=1-OmegaM, Sigma8=0.8, fSigma8=FALSE, ref)
  cosgrowSigma8(z=1, OmegaM=0.3, OmegaL=1-OmegaM, Sigma8=0.8, ref)
  cosgrowFactorApprox(z=1, OmegaM=0.3, OmegaL=1-OmegaM, ref)
  cosgrowRateApprox(z=1, OmegaM=0.3, OmegaL=1-OmegaM, Sigma8=0.8, fSigma8=FALSE, ref)
  cosgrowSigma8Approx(z=1, OmegaM=0.3, OmegaL=1-OmegaM, Sigma8=0.8, ref)
  cosgrowRhoCrit(z=1, H0=100, OmegaM=0.3, OmegaL=1-OmegaM, ref)
  cosgrowRhoMean(z=1, H0=100, OmegaM=0.3, OmegaL=1-OmegaM, ref)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{z}{
    Redshift, where z must be > -1 (can be a vector).
  }
  \item{H0}{
    Hubble constant as defined at z=0 (default is H0=100 (km/s)/Mpc).
  }
  \item{OmegaM}{
    Omega matter (default is 0.3).
  }
  \item{OmegaL}{
    Omega Lambda (default is for a flat Universe with OmegaL=1-OmegaM).
  }
  \item{Sigma8}{
    The value of Sigma8 to use if fsigma8=TRUE (by default this is a reasonable 0.8 for simplicity).
  }
  \item{fSigma8}{
    Logical to express whether the growth rate of structure calculated by cosgrow, cosgrowRate or cosgrowRateApprox is given as f*Sigma8 (TRUE) or simply f (FALSE). This is useful for redshift space distortion comparisons (RSD), since RSD strictly measures f*Sigma8.
  }
  \item{ref}{
    The name of a reference cosmology to use, one of 137 / 737 / Planck / WMAP9 / WMAP7 / WMAP5 / WMAP3 / WMAP1 / Millennium / GiggleZ. See \code{\link{cosref}} for details. This overrides any other settings for H0/ OmegaM and OmegaL. If ref=137 or ref=737 no specific Sigma8 is assumed, instead Sigma8 is set to whatever the input value is set to (by default this is 0.8).
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
Rate    \tab Exact growth rate (f or f*Sigma8, see cosgrowRate below for details)\cr
Sigma8  \tab Power spectrum fluctuation amplitude on the scale 8 Mpc/z\cr
RhoCrit \tab Critical energy density of the Universe at z [3H(z)^2/8.pi.G], in units of Msol/Mpc^3\cr
RhoMean \tab Mean mass density of the Universe at z [RhoCrit.OmegaM(z)], in units of Msol/Mpc^3\cr
}

  cosgrowz function returns the input redshift (only included for clarity).
  
  cosgrowa function returns the Universe expansion factor, as given by a=1/(1+z).
  
  cosgrowH function returns the value of the Hubble constant at z, in units of km/s/Mpc.
  
  cosgrowOmegaM function returns the value of Omega matter at z.
  
  cosgrowOmegaL function returns the value of Omega lambda at z.
  
  cosgrowOmegaK function returns the value of Omega curvature at z.
  
  cosgrowFactor function returns the exact value of the growth factor (typically referred to as 'g' in the astronomy literature), at z.
  
  cosgrowRate function returns either the true (typically referred to as 'f' in the astronomy literature) or RSD type (f*Sigma8) value of the growth rate of structure, at z.
  
  cosgrowSigma8 function returns the power spectrum fluctuation amplitude on the scale 8 Mpc/z at z, and is unitless.
  
  cosgrowFactorApprox function returns the approximate value of the growth factor (typically referred to as 'g' in the astronomy literature), at z.
  
  cosgrowRateApprox function returns either the approximate true (typically referred to as 'f' in the astronomy literature) or approximate RSD type (f*Sigma8) value of the growth rate of structure, at z.
  
  cosgrowSigma8Approx function returns the approximate power spectrum fluctuation amplitude on the scale 8 Mpc/z at z, and is unitless.
  
  cosgrowRhoCrit function returns the critical energy density of the Universe at z [3H(z)^2/8.pi.G], in units of Msol/Mpc^3.
  
  cosgrowRhoMean function returns the mean mass density of the Universe at z [RhoCrit.OmegaM(z)], in units of Msol/Mpc^3.
}
\note{
  RhoCrit and RhoMean return Msol/Mpc^3, where Mpc is in *physical* units (hence both get denser at higher redshift). To get their respective values in comoving units you shoudl simply set z=0, since by definition comoving densities equal physcial densities at z=0. The difference between RhoCrit and RhoMean at z=0 is simply RhoMean=RhoCrit*OmegaM. Corrected for 1/(1+z)^3 RhoMean stays constant with redshift (as it should- we do not gain or lose mass within a cmovoing volume).
}
\references{
  Based on the equations in Hamilton A.J.S., 2001, MNRAS 322 419 and Lahav O., et al., 1991, MNRAS, 251, 136
}
\author{
  Aaron Robotham
}

\seealso{
  \code{\link{cosvol}}, \code{\link{cosmap}},  \code{\link{cosdist}}, \code{\link{cosref}}
}

\examples{
  cosgrow(0.3)
  cosgrow(0.3,ref='Planck')
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
  cosgrowRhoMean(0.3)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{cosmo}
\keyword{cosmology}
\keyword{growth}% __ONLY ONE__ keyword per line

