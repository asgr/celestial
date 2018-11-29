\name{coshalo}
\alias{coshalo}
\alias{coshaloMvirToSigma}
\alias{coshaloSigmaToMvir}
\alias{coshaloMvirToRvir}
\alias{coshaloRvirToMvir}
\alias{coshaloSigmaToRvir}
\alias{coshaloRvirToSigma}
\alias{coshaloSigmaToTvir}

\title{
Virial halo conversion functions
}
\description{
All 6 Virial parameter conversion functions. Each can map precisely to the other as a one paramter function.
}
\usage{
coshaloMvirToSigma(Mvir=1e+12, z=0, H0=100, OmegaM=0.3, OmegaL=1-OmegaM-OmegaR,
OmegaR=0, Rho='crit', Dist='Co', DeltaVir=200, Munit=1, Lunit=1e6, Vunit=1e3, Dim=3, ref)
coshaloSigmaToMvir(Sigma=230, z=0, H0=100, OmegaM=0.3, OmegaL=1-OmegaM-OmegaR,
OmegaR=0, Rho='crit', Dist='Co', DeltaVir=200, Munit=1, Lunit=1e6, Vunit=1e3, Dim=3, ref)
coshaloMvirToRvir(Mvir=1e12, z=0, H0=100, OmegaM=0.3, OmegaL=1-OmegaM-OmegaR,
OmegaR=0, Rho='crit', Dist='Co', DeltaVir=200, Munit=1, Lunit=1e6, Vunit=1e3, Dim=3, ref)
coshaloRvirToMvir(Rvir=162.635, z=0, H0=100, OmegaM=0.3, OmegaL=1-OmegaM-OmegaR,
OmegaR=0, Rho='crit', Dist='Co', DeltaVir=200, Munit=1, Lunit=1e6, Vunit=1e3, Dim=3, ref)
coshaloSigmaToRvir(Sigma=230, z=0, H0=100, OmegaM=0.3, OmegaL=1-OmegaM-OmegaR,
OmegaR=0, Rho='crit', Dist='Co', DeltaVir=200, Munit=1, Lunit=1e6, Vunit=1e3, Dim=3, ref)
coshaloRvirToSigma(Rvir=162.635, z=0, H0=100, OmegaM=0.3, OmegaL=1-OmegaM-OmegaR,
OmegaR=0, Rho='crit', Dist='Co', DeltaVir=200, Munit=1, Lunit=1e6, Vunit=1e3, Dim=3, ref)
coshaloSigmaToTvir(Sigma=230, Vunit=1e3, Tunit='K', Type='halo', Dim=3) 
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Mvir}{
    Mass within virial radius in units of 'Munit'.
  }
  \item{Sigma}{
    Velocity dispersion (3D) within virial radius in units of 'Vunit'. For coshaloSigmaToTvir the Sigma input should be the virial Sigma which can be found by setting DeltaVir='get' in the the other coshalo functions.
  }
  \item{Rvir}{
    Virial radius (3D) in units of 'Lunit'.
  }
  \item{z}{
    Cosmological redshift, where z must be > -1 (can be a vector).
  }
  \item{H0}{
    Hubble constant as defined at z=0 (default is H0=100 (km/s)/Mpc).
  }
  \item{OmegaM}{
    Omega Matter as defined at z=0 (default is 0.3).
  }
  \item{OmegaL}{
    Omega Lambda as defined at z=0 (default is for a flat Universe with OmegaL = 1-OmegaM = 0.7).
  }
  \item{OmegaR}{
    Omega Radiation as defined at z=0 (default is 0, but OmegaM/3400 is typical).
  }
  \item{Rho}{
    Set whether the critical energy density is used (crit) or the mean mass density (mean).
  }
  \item{Dist}{
    Determines the distance type, i.e. whether the Rho critical energy or mean mass densities are calculated with respect to angular / physical distances (Ang) or with respect to comoving distances (Co). In effect this means Rvir values are either angular / physical (Ang) or comoving (Co). It does not affect Mvir <-> Sigma conversions, but does affect Mvir <-> Rvir and Rvir <-> Sigma. 
  }
  \item{DeltaVir}{
    Desired overdensity of the halo with respect to Rho. If set to 'get' it will estimate the required DeltaVir for a virial collapse using the \code{\link{cosgrowDeltaVir}} function.
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
  \item{Type}{
    Specify the 'halo' or 'gas' temperature to be computed.
  }
  \item{Tunit}{
    Specify the output temperature to be Kelvin ('K'), 'eV' or 'keV'.
  }
  \item{Dim}{
    The dimensional type for the halo (either the 2 or 3). 3 (default) means quantities are intrinsic 3D values. 2 means quantities are for projected systems (i.e. radius and velocity dispersion are compressed). From comparisons to simulations (so NFW, c~5 halos) Rvir[proj]=Rvir[3D]/1.37 and Sigma[proj]=Sigma[3D]/sqrt(3). The former has dependence on the halo profile (so is approximate), whereas the latter is a dimensionality argument that should hold for any virialised system. Note that for projected systems Sigma is measured along one dimension: the line-of-site.
  }
  \item{ref}{
    The name of a reference cosmology to use, one of 137 / 737 / Planck / Planck13 / Planck15 / Planck18 / WMAP / WMAP9 / WMAP7 / WMAP5 / WMAP3 / WMAP1 / Millennium / GiggleZ. Planck = Planck18 and WMAP = WMAP9. The usage is case insensitive, so wmap9 is an allowed input. See \code{\link{cosref}} for details. This overrides any other settings for H0/ OmegaM and OmegaL.
  }
}

\details{
These functions allow for various analytic conversions between the 3 major properties related to virial radius: the mass, velocity dispresion and size. The default properties calculate properties for 1e12 Msun halos and assume masses in Msun, velocities in km/s and distances in Kpc.
}

\value{
\item{coshaloMvirToSigma}{Outputs approximate velocity dispersion (in units of Vunit) given mass (this is exactly the escape velocity at Rvir).}

\item{coshaloSigmaToMvir}{Outputs mass (in units of Munit) given velocity dispersion.}

\item{coshaloMvirToRvir}{Outputs radius (in units of Lunit) given mass.}

\item{coshaloRvirToMvir}{Outputs mass (in units of Munit) given radius.}

\item{coshaloSigmaToRvir}{Outputs radius (in units of Lunit) given velocity dispersion.}

\item{coshaloRvirToSigma}{Outputs approximate velocity dispersion (in units of Vunit) given radius (this is exactly the escape velocity at Rvir).}

\item{coshaloSigmaToTvir}{Output temperture (in units of Tunit) given velocity dispersion. Based on Eqns. 3/7/8/9 of Bryan & Norman (1998).}
}

\references{
  coshaloSigmaToTvir based on the equations in:
  
  Bryan & Norman, 1998, ApJ, 495, 80
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

