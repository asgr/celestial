\name{cosvol}
\alias{cosvol}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Cosmological volume calculator
}
\description{
Given the sky area, two redshifts and the cosmology, this function calculates the comoving volume.
}
\usage{
cosvol(area = 60, zmax = 1, zmin = 0, H0 = 100, OmegaM = 0.3, OmegaL = 1 - OmegaM,
inunit = "deg2")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{area}{
Sky area in units of innunit (default is square degrees)
}
  \item{zmax}{
Maximum redshift of comoving cone.
}
  \item{zmin}{
Minimum redshift of comoving cone.
}
  \item{H0}{
Hubble constant (default is H0=100 (km/s)/Mpc)
}
  \item{OmegaM}{
Omega matter (default is 0.3)
}
  \item{OmegaL}{
Omega Lambda (default is for a flat Universe with OmegaL=1-OmegaM)
}
  \item{inunit}{
The units of angular area provided. Allowed options are deg2 for square degrees, amin2 for square arc minutes, asec2 for square arc seconds and rad2 or sr for steradians.
}
}

\value{
A single number specifying the comoving volume of the specified cone segment in Gpc^3.
}
\references{
Based on the equations in Hogg D.W., 1999, arXiv, 9905116 and Wright E.L., 2006, PASP, 118, 1711
}
\author{
Aaron Robotham
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{cosdist}},\code{\link{skyarea}}, \code{\link{cosmap}}
}
\examples{
#Approximate volume of the GAMA survey (area given in skyarea example, zmax is approx
#limit of main galaxy sample):
TotalGAMAvol=cosvol(293.82,0.6)
print(paste('The GAMA survey volume is ~',round(TotalGAMAvol,2),'Gpc^3'))

#Approximate volume of SDSS (area given for DR7, zmax is approx limit of main galaxy sample):
TotalSDSSvol=cosvol(8423,0.3)
print(paste('The SDSS survey volume is ~',round(TotalSDSSvol,2),'Gpc^3'))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{cosmo}
\keyword{cosmology}
\keyword{volume}% __ONLY ONE__ keyword per line
