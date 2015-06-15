\name{cosvar}
\alias{cosvar}
\alias{cosvarcar}
\alias{cosvarsph}
\alias{cosvararea}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Driver & Robotham (2010) cosmic variance calculator
}
\description{
The main cosmic variance calculator function taken from Driver & Robotham (2010).
cosvarcar is an interface to the Cartesian coordinate version, whilst cosvarsph is
a utility interface to give approximate cosmic variance for astronomy survey regions
(usually defined by RA, Dec and redshift limits).
}
\usage{
cosvarcar(aside = 50, bside = 50, cside = 50, regions = 1)
cosvarsph(long = c(129, 141), lat = c(-2, 3), zmax = 1, zmin = 0, regions = 1,
inunit='deg', sep=":")
cosvararea(area=60, zmax=1, zmin=0, regions=1, inunit='deg2')
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{aside}{
The aside (shortest projected side) of the Cartesian box, must be defined using 737 cosmology.
}
  \item{bside}{
The bside (longest projects side) of the Cartesian box, must be defined using 737 cosmology.
}
  \item{cside}{
The cside (radial side) of the Cartesian box, must be defined using 737 cosmology.
}
  \item{regions}{
How many well separated regions of this size will there be?
}
  \item{long}{
Upper and lower longitude (RA) limits of interest in units of inunit. If of length 1 then the number specified is assumed to be the upper limit and the lower limit is set to 0.
}
  \item{lat}{
Upper and lower latitude (Dec) limits of interest in units of inunit. If of length 1 then the number specified is assumed to be the upper limit and the lower limit is set to 0.
}
  \item{zmax}{
Maximum redshift of comoving cone.
}
  \item{zmin}{
Minimum redshift of comoving cone.
}
  \item{inunit}{
    \tabular{ll}{
cosvarsph \tab The units of angular coordinate provided for long and lat (see \code{\link{skyarea}}). \cr
cosvararea \tab The units of angular area provided (see \code{\link{cosvol}}). \cr
    }
}
  \item{sep}{
When inunit='sex', sep defines the type of separator used for the HMS and DMS strings (i.e. H:M:S and D:M:S would be sep=':', which is the default). See \code{\link{hms2deg}} and \code{\link{dms2deg}} for more details.
}
  \item{area}{
    Sky area in units of innunit (default is square degrees)
}
}
\details{
These functions use the empircally motivated cosmic variance percentage formula provided in Driver & Robotham (2010) Eqn 4.

cosvarsph is a 'best effort' approximation of the comoving box subtended by the specified spherical coordinates using the following conversions:

CoDistLow = cosdistCoDist(z=zmin,H0=70,OmegaM=0.3)     
CoDistHigh = cosdistCoDist(z=zmax,H0=70,OmegaM=0.3)
cside=CoDistHigh-CoDistLow
area=skyarea(long = long, lat = lat, inunit = inunit, outunit='deg2')[1]
volume=cosvol(area=area, zmax = zmax, zmin=zmin, H0 = 70, OmegaM = 0.3, inunit='deg2')[1]
aside=cos(mean(lat)*pi/180)*(abs(diff(long))/360)*(CoDistLow+cside/2)
bside=(abs(diff(long))/180)*(CoDistLow+cside/2)
scale=sqrt(volume*1e9/(aside*bside*cside))
aside=aside*scale
bside=bside*scale
return(cosvarcar(aside=aside, bside=bside, cside=cside, subsets=subsets))

cosvararea is a simplifed version of cosvarsph, where the assumption is that aside=bside (so the aspect ratio on the sky is 1:1).
}
\value{
The output is the approximate percentage cosmic (or sample) variance that is expected for the volume specified.
}
\references{
Driver & Robotham, 2010, MNRAS, 407, 2131
}
\author{
Aaron Robotham and Simon Driver
}
\note{
Many people get upset at the term 'cosmic variance' and prefer 'sample variance'. Whilst I am sympathetic to the argument, more astronomers are familiar with the former term.
	
These cosmic variance estimates are defined using SDSS at z~0.1, caveats abound at higher redshifts, but these numbers should serve as a reasonably conservative (i.e. pessimistic) upper limit.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{cosvol}}, \code{\link{skyarea}}
}
\examples{
#Approximate CV of the GAMA equatorial regions:
cosvarsph(long=12, lat=5, zmax=0.5)*1/sqrt(3)
#Or using the GAMA sexigesimal coordinates (should be the same):
cosvarsph(long = c('11:36:0','12:24:0'), lat = c('-2:0:0','3:0:0'), zmax=0.5,
inunit='sex')*1/sqrt(3)
#Approximate CV of the SDSS:
cosvarsph(long=150, lat=100, zmax=0.3)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{cosmic}
\keyword{variance}% __ONLY ONE__ keyword per line
\keyword{sample}
