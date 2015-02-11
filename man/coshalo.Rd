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
All 6 Virial paramter conversion functions. Each can map precisely to the other as a one paramter function.
}
\usage{
coshaloMvirToSigma(Mvir = 1e+12, Munit = 1, Lunit = 1000, Vunit = 1, DeltaVir = 200, H0=100)
coshaloSigmaToMvir(Sigma = 230, Munit = 1, Lunit = 1000, Vunit = 1, DeltaVir = 200, H0=100)
coshaloMvirToRvir(Mvir=1e12, Munit=1, Lunit=1e3, Vunit=1, DeltaVir=200, H0=100)
coshaloRvirToMvir(Rvir=162.635, Munit=1, Lunit=1e3, Vunit=1, DeltaVir=200, H0=100)
coshaloSigmaToRvir(Sigma=230, Munit=1, Lunit=1e3, Vunit=1, DeltaVir=200, H0=100)
coshaloRvirToSigma(Rvir=162.635, Munit=1, Lunit=1e3, Vunit=1, DeltaVir=200, H0=100)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Mvir}{
Mass within virial radius in units of 'Munit'.
}
  \item{Sigma}{
Velocity dispersion within virial radius in units of 'Vunit'.
}
  \item{Rvir}{
Size of virial radius in units of 'Lunit'.
}
  \item{Munit}{
Base mass unit in multiples of Msun.
}
  \item{Lunit}{
Base length unit in multiples of parsecs.
}
  \item{Vunit}{
Base velocity unit in multiples of km/s.
}
  \item{DeltaVir}{
Relative over-density virial radius.
}
  \item{H0}{
Hubble constant as defined at z=0 (default is H0=100 (km/s)/Mpc).
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
 \code{\link{cosvol}}, \code{\link{cosmap}},  \code{\link{cosdist}}, \code{\link{cosgrow}}
}

\examples{
coshaloMvirToSigma(1e13) # Velocity in km/s
coshaloMvirToSigma(1e13, Vunit=1e-3) # Velocity in m/s
coshaloSigmaToMvir(coshaloMvirToSigma(1e13, Vunit=1e-3),Vunit=1e-3)
coshaloMvirToRvir(1e13) #Radius in kpc
coshaloSigmaToRvir(coshaloMvirToSigma(1e13, Vunit=1e-3),Vunit=1e-3)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{halo}
\keyword{cosmo}
\keyword{cosmology}

