\name{cosNFW}
\alias{cosNFW}
\alias{cosNFWmass_c}
\alias{cosNFWmass_Rmax}
\title{
Navaro Frenk and White profile
}
\description{
Density and total mass values for Navaro Frenk and White (NFW) profiles
}
\usage{
cosNFW(Rad=0, Rho0=2.412e5, Rs=0.03253)
cosNFWmass_c(Rho0=2.412e5, Rs=0.03253, c=5)
cosNFWmass_Rmax(Rho0=2.412e5, Rs=0.03253, Rmax=0.1626)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Rad}{
Radius at which the 3D profile denisty should be calculated.
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
The NFW profile Rmax parameter, where Rmax=Rs*c
}
}

\details{
These functions calculate various aspects of the NFW profile.
}

\value{
cosNFW outputs the instantaneous NFW profile density.

coshaloNFWmass_c outputs the total mass given Rs and c in Msun/h.

coshaloNFWmass_Rmax outputs the total mass given Rs and Rmax in Msun/h.
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
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{halo}
\keyword{NFW}
\keyword{cosmo}
\keyword{cosmology}

