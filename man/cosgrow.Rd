\name{cosgrow}
\alias{cosgrow}
\alias{cosgrowz}
\alias{cosgrowa}
\alias{cosgrowOmegaM}
\alias{cosgrowOmegaL}
\alias{cosgrowOmegaR}
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
  Cosmological growth and evolution calculations
}
\description{
  These functions allows comoving, angular size and luminosity distances to be calculated for a given redshift, it can also return look back time. They use curvature correctly, calculated internally using the relation OmegaM+OmegaL+OmegaR+OmegaK=1, but by default they assume a flat Universe where only OmegaM needs to be specified and OmegaR=0 (so no radiation pressure at any epoch).
}
\usage{
  cosgrow(z=1, H0=100, OmegaM=0.3, OmegaL=1-OmegaM-OmegaR, OmegaR=0, Sigma8=0.8,
  fSigma8=FALSE, Dist='Co',
  Mass='Msun', ref)
  cosgrowz(z = 1)
  cosgrowa(z = 1)
  cosgrowH(z=1, H0=100, OmegaM=0.3, OmegaL=1-OmegaM-OmegaR, OmegaR=0, ref)
  cosgrowOmegaM(z=1, OmegaM=0.3, OmegaL=1-OmegaM-OmegaR, OmegaR=0, ref)
  cosgrowOmegaL(z=1, OmegaM=0.3, OmegaL=1-OmegaM-OmegaR, OmegaR=0, ref)
  cosgrowOmegaR(z=1, OmegaM=0.3, OmegaL=1-OmegaM-OmegaR, OmegaR=0, ref)
  cosgrowOmegaK(z=1, OmegaM=0.3, OmegaL=1-OmegaM-OmegaR, OmegaR=0, ref)
  cosgrowFactor(z=1, OmegaM=0.3, OmegaL=1-OmegaM, ref)
  cosgrowRate(z=1, OmegaM=0.3, OmegaL=1-OmegaM, Sigma8=0.8,
  fSigma8=FALSE, ref)
  cosgrowSigma8(z=1, OmegaM=0.3, OmegaL=1-OmegaM, Sigma8=0.8, ref)
  cosgrowFactorApprox(z=1, OmegaM=0.3, OmegaL=1-OmegaM, ref)
  cosgrowRateApprox(z=1, OmegaM=0.3, OmegaL=1-OmegaM, Sigma8=0.8,
  fSigma8=FALSE, ref)
  cosgrowSigma8Approx(z=1, OmegaM=0.3, OmegaL=1-OmegaM, Sigma8=0.8, ref)
  cosgrowRhoCrit(z=1, H0=100, OmegaM=0.3, OmegaL=1-OmegaM-OmegaR, OmegaR=0, Dist='Co',
  Mass='Msun', ref)
  cosgrowRhoMean(z=1, H0=100, OmegaM=0.3, OmegaL=1-OmegaM-OmegaR, OmegaR=0, Dist='Co',
  Mass='Msun', ref)
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
    Omega Lambda (default is for a flat Universe with OmegaL = 1-OmegaM-OmegaR = 0.7).
  }
  \item{OmegaR}{
    Omega Radiation (default is 0, but ~OmegaM/3400 is typical).
  }
  \item{Sigma8}{
    The value of Sigma8 to use if fsigma8=TRUE (by default this is a reasonable 0.8 for simplicity).
  }
  \item{fSigma8}{
    Logical to express whether the growth rate of structure calculated by cosgrow, cosgrowRate or cosgrowRateApprox is given as f*Sigma8 (TRUE) or simply f (FALSE). This is useful for redshift space distortion comparisons (RSD), since RSD strictly measures f*Sigma8.
  }
  \item{Dist}{
    Determines the distance type, i.e. whether the Rho critical energy or mean mass densities are calculated with respect to angular / physical distances (Ang) or with respect to comoving distances (Co).
  }
  \item{Mass}{
    Determines the mass type, i.e. whether Rho ritical energy or mean mass densities are calculated with respect to solar masses (Msun) or with respect to kilograms (kg).
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
\item{cosgrow}{Returns a data.frame (even if only 1 redshift if requested) with the following parameters evaluated at the respective redshift/s:}
\tabular{ll}{
z       \tab Requested redshift\cr
a       \tab Universe expansion factor, as given by a=1/(1+z)\cr
H       \tab Hubble value in units of (km/s)/Mpc\cr
OmegaM  \tab Omega Matter\cr
OmegaL  \tab Omega Lambda\cr
OmegaR  \tab Omega Radiation\cr
OmegaK  \tab Omega K(c)urvature\cr
Factor  \tab Exact growth factor (g, see cosgrowFactor below for details)\cr
Rate    \tab Exact growth rate (f or f*Sigma8, see cosgrowRate below for details)\cr
Sigma8  \tab Power spectrum fluctuation amplitude on the scale 8 Mpc/z\cr
RhoCrit \tab Critical energy density of the Universe at z [3H(z)^2/8.pi.G], in units of Msol/Mpc^3\cr
RhoMean \tab Mean mass density of the Universe at z [RhoCrit.OmegaM(z)], in units of Msol/Mpc^3\cr
}

The outputs of the standalone functions are:

  \item{cosgrowz}{Returns the input redshift (only included for clarity).}
  
  \item{cosgrowa}{Returns the Universe expansion factor, as given by a=1/(1+z).}
  
  \item{cosgrowH}{Returns the value of the Hubble constant at z, in units of km/s/Mpc.}
  
  \item{cosgrowOmegaM}{Returns the value of Omega Matter at z.}
  
  \item{cosgrowOmegaL}{Returns the value of Omega Lambda at z.}
  
  \item{cosgrowOmegaR}{Returns the value of Omega Radiaion at z.}
  
  \item{cosgrowOmegaK}{Returns the value of Omega K(c)urvature at z.}
  
  \item{cosgrowFactor}{Returns the exact value of the growth factor (typically referred to as 'g' in the astronomy literature), at z.}
  
  \item{cosgrowRate}{Returns either the true (typically referred to as 'f' in the astronomy literature) or RSD type (f*Sigma8) value of the growth rate of structure, at z.}
  
  \item{cosgrowSigma8}{Returns the power spectrum fluctuation amplitude on the scale 8 Mpc/z at z, and is unitless.}
  
  \item{cosgrowFactorApprox}{Returns the approximate value of the growth factor (typically referred to as 'g' in the astronomy literature), at z.}
  
  \item{cosgrowRateApprox}{Returns either the approximate true (typically referred to as 'f' in the astronomy literature) or approximate RSD type (f*Sigma8) value of the growth rate of structure, at z.}
  
  \item{cosgrowSigma8Approx}{Returns the approximate power spectrum fluctuation amplitude on the scale 8 Mpc/z at z, and is unitless.}
  
  \item{cosgrowRhoCrit}{Returns the critical energy density of the Universe at z [3H(z)^2/8.pi.G], in units of Msol/Mpc^3.}
  
  \item{cosgrowRhoMean}{Returns the mean mass density of the Universe at z [RhoCrit.OmegaM(z)], in units of Msol/Mpc^3.}
}
\note{
  The difference between RhoCrit and RhoMean at z=0 is simply RhoMean=RhoCrit*OmegaM. Corrected for 1/(1+z)^3 RhoMean stays constant with redshift (as it should- to first order we do not gain or lose mass within a comovoing volume).
  
  The growth rate and growth factor does not make use of OmegaR in the cosgrow function, hence OmegaR cannot be provided in the individual functions. This is because correctly accounting for the effect of radiation pressure before the surface of last scattering (z~1100) on teh growth rate of structure is highly complex, and beyodn the scope of this calculator. In the case of cosgrow, even if OmegaR is specified it is, in effect, set to zero when making growth factor and rate calculations.
}
\references{
  Based on the equations in:
  
  Hamilton A.J.S., 2001, MNRAS 322 419
  
  Lahav O., et al., 1991, MNRAS, 251, 136
}
\author{
  Aaron Robotham
}

\seealso{
  \code{\link{cosvol}}, \code{\link{cosmap}},  \code{\link{cosdist}}, \code{\link{cosref}}, \code{\link{coshalo}}
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
  cosgrowFactorApprox(0.3) #Approximation better than 1% for reasonable cosmologies.
  cosgrowRate(0.3)
  cosgrowRateApprox(0.3) #Approximation better than 1% for reasonable cosmologies.
  cosgrowRhoCrit(0.3)
  cosgrowRhoMean(0.3)
  cosgrowRhoMean(0)-cosgrowRhoMean(2,Dist='Ang')/(1+2)^3 #Mass is conserved in co-vol
  cosgrowRhoMean(0)-cosgrowRhoMean(10,Dist='Co') #Mass is conserved in co-vol
  
  # Comparison of the various energy densities that make up the Universe for Planck 2013:
  
  plot(cosdistUniAgeAtz(10^seq(-3,4.9,by=0.1), ref='Planck')*1e9,
  cosgrowRhoCrit(10^seq(-3,4.9,by=0.1), ref='Planck', Dist='m', Mass='kg')*
  cosgrowOmegaR(10^seq(-3,4.9,by=0.1), ref='Planck'), type='l',log='xy',
  xlab='Years since Universe formed', ylab=expression('Energy Density'*(kg/m^3)))
  
  lines(cosdistUniAgeAtz(10^seq(-3,4.9,by=0.1), ref='Planck')*1e9,
  cosgrowRhoCrit(10^seq(-3,4.9,by=0.1), ref='Planck', Dist='m', Mass='kg')*
  cosgrowOmegaM(10^seq(-3,4.9,by=0.1), ref='Planck'), col='red')
  
  lines(cosdistUniAgeAtz(10^seq(-3,4.9,by=0.1), ref='Planck')*1e9,
  cosgrowRhoCrit(10^seq(-3,4.9,by=0.1), ref='Planck', Dist='m', Mass='kg')*
  cosgrowOmegaL(10^seq(-3,4.9,by=0.1), ref='Planck'), col='blue')

  abline(v=cosdistUniAgeAtz(0.33,ref='Planck')*1e9,lty=2) # Matter = Vacuum
  abline(v=cosdistUniAgeAtz(3391,ref='Planck')*1e9,lty=2) # Matter = Radiation
  
  legend('topright', legend=c('Radiation Energy Density', 'Matter Energy Density',
  'Vacuum Energy Density'), lty=1, col=c('black','red','blue'))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{cosmo}
\keyword{cosmology}
\keyword{growth}% __ONLY ONE__ keyword per line

