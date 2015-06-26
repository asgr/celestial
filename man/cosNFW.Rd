\name{cosNFW}
\alias{cosNFW}
\alias{cosNFWmass_c}
\alias{cosNFWmass_Rmax}
\alias{cosNFWsigma}
\alias{cosNFWsigma_mean}
\alias{cosNFWgamma}
\title{
Navarro Frenk and White profile
}
\description{
Density and total mass values for Navaro Frenk and White (NFW) profiles
}
\usage{
cosNFW(Rad=0, Rho0=2.412e15, Rs=0.03253)
cosNFWmass_c(Rho0=2.412e15, Rs=0.03253, c=5)
cosNFWmass_Rmax(Rho0=2.412e15, Rs=0.03253, Rmax=0.16265)
cosNFWsigma(Rad=0.03253, Rs=0.03253, c=5, z = 0, H0 = 100, OmegaM = 0.3, OmegaL = 1-OmegaM, Rho = "crit", DeltaVir = 200, ref)
cosNFWsigma_mean(Rad=0.03253, Rs=0.03253, c=5, z = 0, H0 = 100, OmegaM = 0.3, OmegaL = 1-OmegaM, Rho = "crit", DeltaVir = 200, ref)
cosNFWgamma(Rad=0.03253, Rs=0.03253, c=5, SigmaC=1, z = 0, H0 = 100, OmegaM = 0.3, OmegaL = 1-OmegaM, Rho = "crit", DeltaVir = 200, ref)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Rad}{
    Radius at which to calculate output. Either this is a 3D radius (cosNFW) or a projected 2D radius (cosNFWsigma/cosNFWsigma_mean).
  }
  \item{Rho0}{
    The normalising factor.
  }
  \item{Rs}{
    The NFW profile scale radius, where Rs=Rmax/c.
  }
  \item{c}{
    The NFW profile concentration parameter, where c=Rmax/Rs.
  }
  \item{Rmax}{
    The NFW profile Rmax parameter, where Rmax=Rs*c.
  }
  \item{SigmaC}{
    The critical surface mass density (when SigmaC=1 we compute the excess surface density / ESD). See \code{\link{cosdistCrit}} for general computation given source and lens redshifts.
  }
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
    Omega Lambda (default is for a flat Universe with OmegaL = 1-OmegaM = 0.7).
  }
  \item{Rho}{
    Set whether the critical energy density is used (crit) or the mean mass density (mean).
  }
  \item{DeltaVir}{
    Desired overdensity of the halo with respect to Rho.
  }
  \item{ref}{
    The name of a reference cosmology to use, one of 137 / 737 / Planck / WMAP9 / WMAP7 / WMAP5 / WMAP3 / WMAP1 / Millennium / GiggleZ. See \code{\link{cosref}} for details. This overrides any other settings for H0/ OmegaM and OmegaL. If ref=137 or ref=737 no specific Sigma8 is assumed, instead Sigma8 is set to whatever the input value is set to (by default this is 0.8).
  }
}

\details{
These functions calculate various aspects of the NFW profile.
}

\value{
cosNFW outputs the instantaneous NFW profile density.

coshaloNFWmass_c outputs the total mass given Rs and c in Msun/h.

coshaloNFWmass_Rmax outputs the total mass given Rs and Rmax in Msun/h.

cosNFWsigma outputs the line-of-sight surface mass density at Rad (Eqn. 11 of Wright \& Brainerd, 2000).

cosNFWsigma_mean outputs the means surface mass density within Rad (Eqn. 13 of Wright \& Brainerd, 2000).

cosNFWgamma outputs the radial dependence of the weak lensing shear (Eqn. 12 of Wright \& Brainerd, 2000).
}

\references{
Navarro, J.F., Frenk, C.S., White, Simon D.M., 1996, ApJ, 462
Wright, C.Oaxaca., Brainerd, T.G., 2000, ApJ, 534
}

\author{
Aaron Robotham
}

\seealso{
 \code{\link{cosvol}}, \code{\link{cosmap}},  \code{\link{cosdist}}, \code{\link{cosgrow}}, \code{\link{coshalo}}
}

\examples{
#What difference do we see if we use the rad_mean200 radius rather than rad_crit200

rad_crit200=coshaloMvirToRvir(1e12,Lunit=1e6)
rad_mean200=coshaloMvirToRvir(1e12,Lunit=1e6,Rho='mean')
cosNFWmass_Rmax(Rmax=rad_crit200) #By construction we should get ~10^12 Msun/h
cosNFWmass_Rmax(Rmax=rad_mean200) #For the same profile this is a factor 1.31 larger

#Shear checks:

plot(10^seq(-2,2,by=0.1),cosNFWgamma(10^seq(-2,2,by=0.1),Rs=0.2,c=10),type='l',log='xy', xlab='R/Rs', ylab='ESD')
legend('topright', legend=c('Rs=0.2','c=10'))

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{halo}
\keyword{NFW}
\keyword{cosmo}
\keyword{cosmology}

