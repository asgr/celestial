\name{cosNFW}
\alias{cosNFW}
\alias{cosNFWmass_c}
\alias{cosNFWmass_Rmax}
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
cosNFWmass_c(Rho0=2.412e15, Rs=0.03253, c=5)
cosNFWmass_Rmax(Rho0=2.412e15, Rs=0.03253, Rmax=0.16265)
cosNFWsigma(Rad=0.03253, Rs=0.03253, c=5, z = 0, H0 = 100, OmegaM = 0.3,
OmegaL = 1-OmegaM-OmegaR, OmegaR=0, Rho = "crit", DeltaVir = 200, ref)
cosNFWsigma_mean(Rad=0.03253, Rs=0.03253, c=5, z = 0, H0 = 100, OmegaM = 0.3,
OmegaL = 1-OmegaM-OmegaR, OmegaR=0, Rho = "crit", DeltaVir = 200, ref)
cosNFWgamma(Rad=0.03253, Rs=0.03253, c=5, SigmaC=1, z = 0, H0 = 100,
OmegaM = 0.3, OmegaL = 1-OmegaM-OmegaR, OmegaR=0, Rho = "crit", DeltaVir = 200, ref)
cosNFWduffym2c(M=2e12, z = 0, H0 = 100, OmegaM = 0.3, OmegaL = 1-OmegaM-OmegaR,
OmegaR=0, Rho = "crit", A=6.71, B=-0.091, C=-0.44, ref) 
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
  \item{M}{
    The halo mass required for computing the Duffy (2008) mass to concentration conversion. Here the halo mass required for input is the 200 times overdense with respect to critical variation.
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
  \item{OmegaR}{
    Omega Radiation (default is 0, but OmegaM/3400 is typical).
  }
  \item{Rho}{
    Set whether the critical energy density is used (crit) or the mean mass density (mean).
  }
  \item{DeltaVir}{
    Desired overdensity of the halo with respect to Rho.
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
    The name of a reference cosmology to use, one of 137 / 737 / Planck / WMAP9 / WMAP7 / WMAP5 / WMAP3 / WMAP1 / Millennium / GiggleZ. See \code{\link{cosref}} for details. This overrides any other settings for H0/ OmegaM and OmegaL. If ref=137 or ref=737 no specific Sigma8 is assumed, instead Sigma8 is set to whatever the input value is set to (by default this is 0.8).
  }
}

\details{
These functions calculate various aspects of the NFW profile.
}

\value{
\item{cosNFW}{Returns the instantaneous NFW profile density.}

\item{cosNFWmass_c}{Returns the total mass given Rs and c in Msun/h.}

\item{cosNFWmass_Rmax}{Returns the total mass given Rs and Rmax in Msun/h.}

\item{cosNFWsigma}{Returns the line-of-sight surface mass density at Rad (Eqn. 11 of Wright \& Brainerd, 2000).}

\item{cosNFWsigma_mean}{Returns the means surface mass density within Rad (Eqn. 13 of Wright \& Brainerd, 2000).}

\item{cosNFWgamma}{Returns the radial dependence of the weak lensing shear (Eqn. 12 of Wright \& Brainerd, 2000).}

\item{cosNFWduffym2c}{Returns the Duffy (2008) predicted concentration for a given halo mass given .}
}

\references{
Duffy A.R., et al., 2008, MNRAS, 390L

Navarro J.F., Frenk C.S., White Simon D.M., 1996, ApJ, 462

Wright C.O., Brainerd T.G., 2000, ApJ, 534
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
crit=seq(2, 20, by=0.01)
critcol=rainbow(length(crit), start=0, end=5/6)
rad_crit200=coshaloMvirToRvir(1, z=zseq, Rho='crit', DeltaVir=200, ref='Planck15')
rad_crit500=coshaloMvirToRvir(1, z=zseq, Rho='crit', DeltaVir=500, ref='Planck15')
rad_mean200=coshaloMvirToRvir(1, z=zseq, Rho='mean', DeltaVir=200, ref='Planck15')
rad_mean500=coshaloMvirToRvir(1, z=zseq, Rho='mean', DeltaVir=500, ref='Planck15')

plot(1, 1, type='n', xlim=c(0.01,10), ylim=c(1,1.55), xlab='Redshift',
ylab='M200c/M500c', log='x')
for(i in 1:length(crit)){
lines(zseq, cosNFWmass_Rmax(Rho0=1, Rs=rad_crit200[1]/crit[i], Rmax=rad_crit200)/
cosNFWmass_Rmax(Rho0=1, Rs=rad_crit200[1]/crit[i], Rmax=rad_crit500), col=critcol[i])
}

plot(1, 1, type='n', xlim=c(0.01,10), ylim=c(1,1.55), xlab='Redshift',
ylab='M200m/M500m', log='x')
for(i in 1:length(crit)){
lines(zseq, cosNFWmass_Rmax(Rho0=1, Rs=rad_crit200[1]/crit[i], Rmax=rad_mean200)/
cosNFWmass_Rmax(Rho0=1, Rs=rad_crit200[1]/crit[i], Rmax=rad_mean500), col=critcol[i])
}

plot(1, 1, type='n', xlim=c(0.01,10), ylim=c(1,1.55), xlab='Redshift',
ylab='M200m/M200c',log='x')
for(i in 1:length(crit)){
lines(zseq, cosNFWmass_Rmax(Rho0=1, Rs=rad_crit200[1]/crit[i], Rmax=rad_mean200)/
cosNFWmass_Rmax(Rho0=1, Rs=rad_crit200[1]/crit[i], Rmax=rad_crit200), col=critcol[i])
}

plot(1, 1, type='n', xlim=c(0.01,10), ylim=c(1,1.55), xlab='Redshift',
ylab='M500m/M500c', log='x')
for(i in 1:length(crit)){
lines(zseq, cosNFWmass_Rmax(Rho0=1, Rs=rad_crit200[1]/crit[i], Rmax=rad_mean500)/
cosNFWmass_Rmax(Rho0=1, Rs=rad_crit200[1]/crit[i], Rmax=rad_crit500), col=critcol[i])
}

plot(zseq, rad_crit200/rad_crit500, type='l', xlim=c(0.01,10), ylim=c(1,1.55),
xlab='Redshift', ylab='R200/R500', log='x')

plot(zseq, rad_mean200/rad_crit200, type='l', xlim=c(0.01,10), ylim=c(1,1.55), 
xlab='Redshift', ylab='Rmean/Rcrit', log='x')
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{halo}
\keyword{NFW}
\keyword{cosmo}
\keyword{cosmology}

