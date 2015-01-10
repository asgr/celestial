\name{cosdist}
\alias{cosdist}
\alias{Hz}
\alias{rhocrit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Cosmological distance calculator
}
\description{
This function allows comoving, angular size and luminosity distances to be calculated for a given redshift, it can also return look back time. It makes use of curvature properly using the relation OmegaM+OmegaL+OmegaK=1, but by default it assumes a flat Universe where only OmegaM needs to be specified.
}
\usage{
cosdist(z = 1, H0 = 100, OmegaM = 0.3, OmegaL = 1 - OmegaM, age = FALSE)
Hz(z = 1, H0 = 100, OmegaM = 0.3, OmegaL = 1 - OmegaM)
rhocrit(z = 0, H0 = 100, OmegaM = 0.3, OmegaL = 1 - OmegaM)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{z}{
Redshift (can be a vector).
}
  \item{H0}{
Hubble constant (default is H0=100 (km/s)/Mpc).
}
  \item{OmegaM}{
Omega matter (default is 0.3).
}
  \item{OmegaL}{
Omega Lambda (default is for a flat Universe with OmegaL=1-OmegaM).
}
  \item{age}{
Flag to tell program to return age or not- this slows calculation, so is by default turned off.
}
}
\details{
Function is based on D. W. Hogg et all 1999 (arXiv 9905116)
}
\value{
cosdist function returns a data.frame (even if only 1 redshift if requested) with the following columns:
\tabular{ll}{
z       \tab Requested redshift\cr
CoDist  \tab Line-of-sight (i.e. radial) comoving distance in units of Mpc\cr
LumDist \tab Luminosity distance in units of Mpc\cr
AngDist \tab Angular diameter distance in units of Mpc\cr
CoDistTran  \tab Transverse comoving distance in units of Mpc\cr
DistMod \tab The distance modulus used where AbsMag = ApMag - DistMod, and DistMod = 5log10(LumDist)+25\cr
AngSize \tab Physical projected size of an object at z in units of kpc/arcsec\cr
CoVolGpc3 \tab Comoving volume of Universe within z in units of Gpc^3\cr
}
If age=TRUE is set then three extra pieces of age information are calculated for each z as extra columns:
\tabular{ll}{
HubTime \tab Approximate Hubble age of the Universe in units of years \cr
UniAgeNow  \tab Age of the Universe now in units of years\cr
UniAgeAtz  \tab Age of the Universe at the specified redshift (z) in units of years\cr
TravelTime \tab Light travel time from the specified redshift (AKA look back time) in units of years\cr
}

The Hz function returns the Hubble constant for the specified redshift and cosmology, in units of km/s/Mpc.

The rhocrit function returns the critical density of the Universe for the specified redshift and cosmology, in units of Msol/Mpc^3.
}
\references{
Based on the equations in Hogg D.W., 1999, arXiv, 9905116 and Wright E.L., 2006, PASP, 118, 1711
}
\author{
Aaron Robotham
}

\seealso{
\code{\link{cosvol}}
}

\examples{
print(cosdist(0.3,70))

redshifts=seq(0,3,by=0.01)

plot(redshifts, cosdist(redshifts, 70)[,'CoDist'], type='l', col='darkgreen', xlab='Redshift / z',
ylab='Distance / Mpc')
lines(redshifts, cosdist(redshifts, 70)[,'LumDist'], col='red')
lines(redshifts, cosdist(redshifts, 70)[,'AngDist'], col='blue')
legend('topleft', legend=c('Comoving Distance', 'Luminosity Distance', 'Angular Diameter Distance'),
col=c('darkgreen', 'red', 'blue'),lty=1)

plot(redshifts, cosdist(redshifts, 70, age=TRUE)[,'TravelTime'], type='l', xlab='Redshift / z',
ylab='Light travel time / Yrs')
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{cosmo}
\keyword{cosmology}
\keyword{distance}% __ONLY ONE__ keyword per line

