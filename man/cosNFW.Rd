\name{cosNFW}
\alias{cosNFW}
\alias{cosNFWmass_c}
\alias{cosNFWmass_Rmax}
\alias{cosNFWvcirc}
\alias{cosNFWvesc}
\alias{cosNFWsigma}
\alias{cosNFWsigma_mean}
\alias{cosNFWgamma}
\alias{cosNFWduffym2c}
\title{
Navarro Frenk and White profile
}
\description{
Density and total mass values for Navaro Frenk and White (NFW) profiles
}
\usage{
cosNFW(Rad=0, Rho0=2.412e15, Rs=0.03253)
cosNFWmass_c(Rho0=2.412e15, Rs=0.03253, c=5, Munit = 1, Lunit = 1e+06)
cosNFWmass_Rmax(Rho0=2.412e15, Rs=0.03253, Rmax=0.16265, Munit = 1, Lunit = 1e+06)
cosNFWvcirc(Rad = 0.16264, Mvir = 1e+12, c = 5, f = Inf, z = 0, H0 = 100, OmegaM = 0.3,
OmegaL = 1 - OmegaM - OmegaR, OmegaR = 0, Rho = "crit", Dist = "Co", DeltaVir = 200,
Munit = 1, Lunit = 1e+06, Vunit = 1000, ref)
cosNFWvesc(Rad = 0.16264, Mvir = 1e+12, c = 5, f = Inf, z = 0, H0 = 100, OmegaM = 0.3,
OmegaL = 1 - OmegaM - OmegaR, OmegaR = 0, Rho = "crit", Dist = "Co", DeltaVir = 200,
Munit = 1, Lunit = 1e+06, Vunit = 1000, ref)
cosNFWsigma(Rad=0.03253, Rs=0.03253, c=5, z = 0, H0 = 100, OmegaM = 0.3,
OmegaL = 1-OmegaM-OmegaR, OmegaR=0, Rho = "crit", DeltaVir = 200, Munit = 1,
Lunit = 1e+06, Vunit = 1000, ref)
cosNFWsigma_mean(Rad=0.03253, Rs=0.03253, c=5, z = 0, H0 = 100, OmegaM = 0.3,
OmegaL = 1-OmegaM-OmegaR, OmegaR=0, Rho = "crit", DeltaVir = 200, Munit = 1,
Lunit = 1e+06, Vunit = 1000, ref)
cosNFWgamma(Rad=0.03253, Rs=0.03253, c=5, SigmaC=1, z = 0, H0 = 100,
OmegaM = 0.3, OmegaL = 1-OmegaM-OmegaR, OmegaR=0, Rho = "crit", DeltaVir = 200,
Munit = 1, Lunit = 1e+06, Vunit = 1000, ref)
cosNFWduffym2c(M=2e12, z = 0, H0 = 100, OmegaM = 0.3, OmegaL = 1-OmegaM-OmegaR,
OmegaR=0, Rho = "crit", A=6.71, B=-0.091, C=-0.44, Munit = 1, ref)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Mvir}{
    Mass within virial radius in units of 'Munit'.
  }
  \item{Rad}{
    Radius at which to calculate output in units of 'Lunit'. Either this is a 3D radius (cosNFW) or a projected 2D radius (cosNFWsigma/cosNFWsigma_mean).
  }
  \item{Rho0}{
    The normalising factor.
  }
  \item{Rs}{
    The NFW profile scale radius, where Rs=Rmax/c, in units of 'Munit'.
  }
  \item{c}{
    The NFW profile concentration parameter, where c=Rmax/Rs.
  }
  \item{f}{
    The NFW profile truncation radius in units of \option{Rmax}.
  }
  \item{Rmax}{
    The NFW profile Rmax parameter, where Rmax=Rs*c, in units of 'Lunit'.
  }
  \item{SigmaC}{
    The critical surface mass density (when SigmaC=1 we compute the excess surface density / ESD). See \code{\link{cosdistCrit}} for general computation given source and lens redshifts.
  }
  \item{M}{
    The halo mass required for computing the Duffy (2008) mass to concentration conversion in units of 'Munit'. Here the halo mass required for input is the 200 times overdense with respect to critical variation.
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
  \item{A}{
    Parameter used for Duffy mass to concentration relation.
  }
  \item{B}{
    Parameter used for Duffy mass to concentration relation.
  }
  \item{C}{
    Parameter used for Duffy mass to concentration relation.
  }
  \item{ref}{
    The name of a reference cosmology to use, one of 137 / 737 / Planck / Planck13 / Planck15 / Planck18 / WMAP / WMAP9 / WMAP7 / WMAP5 / WMAP3 / WMAP1 / Millennium / GiggleZ. Planck = Planck18 and WMAP = WMAP9. The usage is case insensitive, so wmap9 is an allowed input. See \code{\link{cosref}} for details. This overrides any other settings for H0/ OmegaM and OmegaL.
  }
}

\details{
These functions calculate various aspects of the NFW profile.
}

\value{
\item{cosNFW}{Returns the instantaneous NFW profile density.}

\item{cosNFWmass_c}{Returns the total mass given Rs and c in Msun/h.}

\item{cosNFWmass_Rmax}{Returns the total mass given Rs and Rmax in Msun/h.}

\item{cosNFWvcirc}{Returns the circular Keplarian orbit velocity for a given radius assuming an NFW halo potential.}

\item{cosNFWvesc}{Returns the minimum escape (or unbinding) velocity for a given radius assuming an NFW halo potential.}

\item{cosNFWsigma}{Returns the line-of-sight surface mass density at Rad (Eqn. 11 of Wright & Brainerd, 2000).}

\item{cosNFWsigma_mean}{Returns the means surface mass density within Rad (Eqn. 13 of Wright & Brainerd, 2000).}

\item{cosNFWgamma}{Returns the radial dependence of the weak lensing shear (Eqn. 12 of Wright & Brainerd, 2000).}

\item{cosNFWduffym2c}{Returns the Duffy et al (2008) predicted concentration for a given halo mass.}
}

\references{
Duffy A.R., et al., 2008, MNRAS, 390L

Navarro J.F., Frenk C.S., White Simon D.M., 1996, ApJ, 462

Wright C.O. & Brainerd T.G., 2000, ApJ, 534
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

plot(10^seq(-2,2,by=0.1), cosNFWgamma(10^seq(-2,2,by=0.1),Rs=0.2,c=10), type='l',
log='xy', xlab='R/Rs', ylab='ESD')
legend('topright', legend=c('Rs=0.2','c=10'))

#How do critical, mean 200 and 500 masses evolve with redshift? Let's see:

zseq=10^seq(-2, 1, by=0.1)
con=seq(2, 20, by=0.01)
concol=rainbow(length(con), start=0, end=5/6)
rad_crit200=coshaloMvirToRvir(1, z=zseq, Rho='crit', DeltaVir=200, ref='Planck15')
rad_crit500=coshaloMvirToRvir(1, z=zseq, Rho='crit', DeltaVir=500, ref='Planck15')
rad_mean200=coshaloMvirToRvir(1, z=zseq, Rho='mean', DeltaVir=200, ref='Planck15')
rad_mean500=coshaloMvirToRvir(1, z=zseq, Rho='mean', DeltaVir=500, ref='Planck15')
rad_vir=coshaloMvirToRvir(1, z=zseq, Rho='crit', DeltaVir='get', ref='Planck15')

plot(1, 1, type='n', xlim=c(0.01,10), ylim=c(0.8,1.55), xlab='Redshift',
ylab='M200c / M500c', log='x')
for(i in 1:length(con)){
lines(zseq, cosNFWmass_Rmax(Rho0=1, Rs=rad_crit200[1]/con[i], Rmax=rad_crit200)/
cosNFWmass_Rmax(Rho0=1, Rs=rad_crit200[1]/con[i], Rmax=rad_crit500), col=concol[i])
}

plot(1, 1, type='n', xlim=c(0.01,10), ylim=c(0.8,1.55), xlab='Redshift',
ylab='M200m / M500m', log='x')
for(i in 1:length(con)){
lines(zseq, cosNFWmass_Rmax(Rho0=1, Rs=rad_crit200[1]/con[i], Rmax=rad_mean200)/
cosNFWmass_Rmax(Rho0=1, Rs=rad_crit200[1]/con[i], Rmax=rad_mean500), col=concol[i])
}

plot(1, 1, type='n', xlim=c(0.01,10), ylim=c(0.8,1.55), xlab='Redshift',
ylab='M200m / M200c',log='x')
for(i in 1:length(con)){
lines(zseq, cosNFWmass_Rmax(Rho0=1, Rs=rad_crit200[1]/con[i], Rmax=rad_mean200)/
cosNFWmass_Rmax(Rho0=1, Rs=rad_crit200[1]/con[i], Rmax=rad_crit200), col=concol[i])
}

plot(1, 1, type='n', xlim=c(0.01,10), ylim=c(0.8,1.55), xlab='Redshift',
ylab='M500m / M500c', log='x')
for(i in 1:length(con)){
lines(zseq, cosNFWmass_Rmax(Rho0=1, Rs=rad_crit200[1]/con[i], Rmax=rad_mean500)/
cosNFWmass_Rmax(Rho0=1, Rs=rad_crit200[1]/con[i], Rmax=rad_crit500), col=concol[i])
}

plot(1, 1, type='n', xlim=c(0.01,10), ylim=c(0.8,1.55), xlab='Redshift',
ylab='Mvir / M200c',log='x')
for(i in 1:length(con)){
lines(zseq, cosNFWmass_Rmax(Rho0=1, Rs=rad_crit200[1]/con[i], Rmax=rad_vir)/
cosNFWmass_Rmax(Rho0=1, Rs=rad_crit200[1]/con[i], Rmax=rad_crit200), col=concol[i])
}

plot(1, 1, type='n', xlim=c(0.01,10), ylim=c(0.8,1.55), xlab='Redshift',
ylab='Mvir / M200m',log='x')
for(i in 1:length(con)){
lines(zseq, cosNFWmass_Rmax(Rho0=1, Rs=rad_crit200[1]/con[i], Rmax=rad_vir)/
cosNFWmass_Rmax(Rho0=1, Rs=rad_crit200[1]/con[i], Rmax=rad_mean200), col=concol[i])
}

plot(zseq, rad_crit200/rad_crit500, type='l', xlim=c(0.01,10), ylim=c(0.8,1.55),
xlab='Redshift', ylab='R200 / R500', log='x')

plot(zseq, rad_mean200/rad_crit200, type='l', xlim=c(0.01,10), ylim=c(0.8,1.55), 
xlab='Redshift', ylab='Rm / Rc', log='x')

plot(zseq, rad_vir/rad_crit200, type='l', xlim=c(0.01,10), ylim=c(0.8,1.55), 
xlab='Redshift', ylab='Rvir / R200c', log='x')

plot(zseq, rad_vir/rad_mean200, type='l', xlim=c(0.01,10), ylim=c(0.8,1.55), 
xlab='Redshift', ylab='Rvir / R200c', log='x')

#R200m and R200c go either side of Rvir, so by cosmic conspiracy the mean is nearly flat:

plot(zseq, 2*rad_vir/(rad_mean200+rad_crit200), type='l', xlim=c(0.01,10),
ylim=c(0.8,1.55), xlab='Redshift', ylab='2Rvir / (R200c+R200m)', log='x')

#To check Vcirc and Vesc for a 10^12 Msun halo:

plot(0:400, cosNFWvcirc(0:400,f=1,Lunit=1e3), type='l', lty=1, xlab='R / kpc',
ylab='V / km/s', ylim=c(0,500))
lines(0:400, cosNFWvesc(0:400,f=1,Lunit=1e3), lty=2)
legend('topright', legend=c('Vel-Circ','Vel-Escape'), lty=c(1,2))
abline(v=coshaloMvirToRvir(Lunit=1e3), lty=3)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{halo}
\keyword{NFW}
\keyword{cosmo}
\keyword{cosmology}

