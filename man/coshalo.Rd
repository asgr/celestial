\name{coshalo}
\alias{coshalo}
\alias{coshaloMvirToSigma}
\alias{coshaloSigmaToMvir}
\alias{coshaloMvirToRvir}
\alias{coshaloRvirToMvir}
\alias{coshaloSigmaToRvir}
\alias{coshaloRvirToSigma}

\title{
Virial halo conversion functions
}
\description{
All 6 Virial parameter conversion functions. Each can map precisely to the other as a one paramter function.
}
\usage{
coshaloMvirToSigma(Mvir=1e+12, z=0, H0=100, OmegaM=0.3, OmegaL=1-OmegaM-OmegaR,
OmegaR=0, Rho='crit', Dist='Co', DeltaVir=200, Munit=1, Lunit=1e6, Vunit=1e3, ref)
coshaloSigmaToMvir(Sigma=230, z=0, H0=100, OmegaM=0.3, OmegaL=1-OmegaM-OmegaR,
OmegaR=0, Rho='crit', Dist='Co', DeltaVir=200, Munit=1, Lunit=1e6, Vunit=1e3, ref)
coshaloMvirToRvir(Mvir=1e12, z=0, H0=100, OmegaM=0.3, OmegaL=1-OmegaM-OmegaR,
OmegaR=0, Rho='crit', Dist='Co', DeltaVir=200, Munit=1, Lunit=1e6, Vunit=1e3, ref)
coshaloRvirToMvir(Rvir=162.635, z=0, H0=100, OmegaM=0.3, OmegaL=1-OmegaM-OmegaR,
OmegaR=0, Rho='crit', Dist='Co', DeltaVir=200, Munit=1, Lunit=1e6, Vunit=1e3, ref)
coshaloSigmaToRvir(Sigma=230, z=0, H0=100, OmegaM=0.3, OmegaL=1-OmegaM-OmegaR,
OmegaR=0, Rho='crit', Dist='Co', DeltaVir=200, Munit=1, Lunit=1e6, Vunit=1e3, ref)
coshaloRvirToSigma(Rvir=162.635, z=0, H0=100, OmegaM=0.3, OmegaL=1-OmegaM-OmegaR,
OmegaR=0, Rho='crit', Dist='Co', DeltaVir=200, Munit=1, Lunit=1e6, Vunit=1e3, ref)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Mvir}{
    Mass within virial radius in units of 'Munit'.
  }
  \item{Sigma}{
    Velocity dispersion (3D) within virial radius in units of 'Vunit'.
  }
  \item{Rvir}{
    Virial radius (3D) in units of 'Lunit'.
  }
  \item{z}{
    Redshift, where z must be > -1 (can be a vector).
  }
  \item{H0}{
    Hubble constant as defined at z=0 (default is H0=100 (km/s)/Mpc).
  }
  \item{OmegaM}{
    Omega Matter (default is 0.3).
  }
  \item{OmegaL}{
    Omega Lambda (default is for a flat Universe with OmegaL = 1-OmegaM = 0.7).
  }
  \item{OmegaR}{
    Omega Radiation (default is 0, but OmegaM/3400 is typical).
  }
  \item{Rho}{
    Set whether the critical energy density is used (crit) or the mean mass density (mean).
  }
  \item{Dist}{
    Determines the distance type, i.e. whether the Rho critical energy or mean mass densities are calculated with respect to angular / physical distances (Ang) or with respect to comoving distances (Co). In effect this means Rvir values are either angular / physical (Ang) or comoving (Co). It does not affect Mvir <-> Sigma conversions, but does affect Mvir <-> Rvir and Rvir <-> Sigma. 
  }
  \item{DeltaVir}{
    Desired overdensity of the halo with respect to Rho.
  }
  \item{Munit}{
    Base mass unit in multiples of Msun.
  }
  \item{Lunit}{
    Base length unit in multiples of parsecs.
  }
  \item{Vunit}{
    Base velocity unit in multiples of m/s.
  }
  \item{ref}{
    The name of a reference cosmology to use, one of 137 / 737 / Planck / WMAP9 / WMAP7 / WMAP5 / WMAP3 / WMAP1 / Millennium / GiggleZ. See \code{\link{cosref}} for details. This overrides any other settings for H0/ OmegaM and OmegaL. If ref=137 or ref=737 no specific Sigma8 is assumed, instead Sigma8 is set to whatever the input value is set to (by default this is 0.8).
  }
}

\details{
These functions allow for various analytic conversions between the 3 major properties related to virial radius: the mass, velocity dispresion and size. The default properties calculate properties for 1e12 Msun halos and assume masses in Msun, velocities in km/s and distances in Kpc.
}

\value{
coshaloMvirToSigma outputs velocity (in units of Vunit) dispersion given mass.

coshaloSigmaToMvir outputs mass (in units of Munit) given velocity dispersion.

coshaloMvirToRvir outputs radius (in units of Lunit) given mass.

coshaloRvirToMvir outputs mass (in units of Munit) given radius.

coshaloSigmaToRvir outputs radius (in units of Lunit) given velocity dispersion.

coshaloRvirToSigma outputs velocity (in units of Vunit) dispersion given radius.
}

\author{
Aaron Robotham, Chris Power
}

\seealso{
 \code{\link{cosvol}}, \code{\link{cosmap}},  \code{\link{cosdist}}, \code{\link{cosgrow}}, \code{\link{cosNFW}}
}

\examples{
coshaloMvirToSigma(1e13) # Velocity in km/s
coshaloMvirToSigma(1e13, Vunit=1) # Velocity in m/s
coshaloSigmaToMvir(coshaloMvirToSigma(1e13, Vunit=1),Vunit=1)
coshaloMvirToRvir(1e13) #Radius in kpc
coshaloSigmaToRvir(coshaloMvirToSigma(1e13, Vunit=1),Vunit=1)

#Some sanity checks

rho_crit200=cosgrowRhoCrit(z=0)*200 #200 times rho critical at z=0
rho_mean200=cosgrowRhoMean(z=0)*200 #200 times rho mean at z=0
#For a 10^12 Msun/h halo, the radius in Mpc/h where the contained density equals rho_crit*200
rad_crit200=(1e12/rho_crit200*3/4/pi)^(1/3)
coshaloMvirToRvir(1e12,Lunit=1e6)-rad_crit200 # ~0 as expected
#For a 10^12 Msun/h halo, the radius in Mpc/h where the contained density equals rho_crit*200
rad_mean200=(1e12/rho_mean200*3/4/pi)^(1/3) # ~0 as expected
coshaloMvirToRvir(1e12,Lunit=1e6,Rho='mean')-rad_mean200
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{halo}
\keyword{cosmo}
\keyword{cosmology}

