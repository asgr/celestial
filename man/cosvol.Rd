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
cosvol(area=60, zmax=1, zmin=0, H0=100, OmegaM=0.3, OmegaL=1-OmegaM-OmegaR, OmegaR=0,
w0=-1, wprime=0, inunit = "deg2", ref)
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
    Hubble constant as defined at z=0 (default is H0=100 (km/s)/Mpc)
  }
  \item{OmegaM}{
    Omega Matter (default is 0.3).
  }
  \item{OmegaL}{
    Omega Lambda (default is for a flat Universe with OmegaL = 1-OmegaM-OmegaR = 0.7).
  }
  \item{OmegaR}{
    Omega Radiation (default is 0, but OmegaM/3400 is typical).
  }
  \item{w0}{
    The value of dark energy equation of state at z=0. See \code{\link{cosgrow}} for more details.
  }
  \item{wprime}{
    The evolution term that governs how the dark energy equation of state evolves with redshift. See \code{\link{cosgrow}} for more details.
  }
  \item{inunit}{
    The units of angular area provided. Allowed options are deg2 for square degrees, amin2 for square arc minutes, asec2 for square arc seconds and rad2 or sr for steradians.
}
  \item{ref}{
    The name of a reference cosmology to use, one of 137 / 737 / Planck / Planck13 / Planck15 / WMAP / WMAP9 / WMAP7 / WMAP5 / WMAP3 / WMAP1 / Millennium / GiggleZ. Planck=Planck15 and WMAP=WMAP9. The usage is case insensitive, so wmap9 is an allowed input. See \code{\link{cosref}} for details. This overrides any other settings for H0/ OmegaM and OmegaL.
  }
}

\value{
A 3 element vector. The first element (voltot) specifies the comoving volume of the requested cone segment in Gpc^3, the second element (volmeanz) specifies the mean redshift when mass is uniformly distributed in the volume, the third element (volmedz) specifies the median redshift when mass is uniformly distributed in the volume.
}
\references{
Based on the equations in:

Hogg D.W., 1999, arXiv, 9905116

Wright E.L., 2006, PASP, 118, 1711
}
\author{
Aaron Robotham
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{cosdist}},\code{\link{skyarea}}, \code{\link{cosmap}}, \code{\link{cosgrow}}
}
\examples{
#Approximate volume of the GAMA survey (area given in skyarea example, zmax is approx
#limit of main galaxy sample):
TotalGAMAvol=cosvol(293.82,0.6)[1]
print(paste('The GAMA survey volume is ~',round(TotalGAMAvol,2),'Gpc^3'))

#Approximate volume of SDSS (area given for DR7, zmax is approx limit of main galaxy sample):
TotalSDSSvol=cosvol(8423,0.3)[1]
print(paste('The SDSS survey volume is ~',round(TotalSDSSvol,2),'Gpc^3'))

#Change of reference cosmology
cosvol(293.82,0.6,ref='Planck')
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{cosmo}
\keyword{cosmology}
\keyword{volume}% __ONLY ONE__ keyword per line
