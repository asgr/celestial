\name{cosdist}
\alias{cosdist}
\alias{cosdistz}
\alias{cosdista}
\alias{cosdistCoDist}
\alias{cosdistLumDist}
\alias{cosdistAngDist}
\alias{cosdistCoDistTran}
\alias{cosdistDistMod}
\alias{cosdistAngSize}
\alias{cosdistCoVol}
\alias{cosdistHubTime}
\alias{cosdistUniAgeNow}
\alias{cosdistUniAgeAtz}
\alias{cosdistTravelTime}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Cosmological Distance Calculations
}
\description{
These functions allows comoving, angular size and luminosity distances to be calculated for a given redshift, it can also return look back time. They use curvature correctly, calculated internally using the relation OmegaM+OmegaL+OmegaK=1, but by default they assume a flat Universe where only OmegaM needs to be specified.
}
\usage{
cosdist(z = 1, H0 = 100, OmegaM = 0.3, OmegaL = 1 - OmegaM, age = FALSE, error = FALSE)
cosdistz(z = 1)
cosdista(z = 1)
cosdistCoDist(z = 1, H0 = 100, OmegaM = 0.3, OmegaL = 1 - OmegaM)
cosdistLumDist(z = 1, H0 = 100, OmegaM = 0.3, OmegaL = 1 - OmegaM)
cosdistAngDist(z = 1, H0 = 100, OmegaM = 0.3, OmegaL = 1 - OmegaM)
cosdistCoDistTran(z = 1, H0 = 100, OmegaM = 0.3, OmegaL = 1 - OmegaM)
cosdistDistMod(z = 1, H0 = 100, OmegaM = 0.3, OmegaL = 1 - OmegaM)
cosdistAngSize(z = 1, H0 = 100, OmegaM = 0.3, OmegaL = 1 - OmegaM)
cosdistCoVol(z = 1, H0 = 100, OmegaM = 0.3, OmegaL = 1 - OmegaM)
cosdistHubTime(H0 = 100)
cosdistUniAgeNow(z = 1, H0 = 100, OmegaM = 0.3, OmegaL = 1 - OmegaM)
cosdistUniAgeAtz(z = 1, H0 = 100, OmegaM = 0.3, OmegaL = 1 - OmegaM)
cosdistTravelTime(z = 1, H0 = 100, OmegaM = 0.3, OmegaL = 1 - OmegaM)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{z}{
Redshift (can be a vector).
}
  \item{H0}{
Hubble constant as defined at z=0 (default is H0=100 (km/s)/Mpc).
}
  \item{OmegaM}{
Omega matter (default is 0.3).
}
  \item{OmegaL}{
Omega Lambda (default is for a flat Universe with OmegaL=1-OmegaM).
}
  \item{age}{
Flag for cosdist function to return age or not- this slows calculation, so is by default turned off.
}
  \item{error}{
Flag for cosdist to calculate the relative error for distance/age values.
}
}
\details{
Function is based on D. W. Hogg et all 1999 (arXiv 9905116)
}
\value{
cosdist function returns a data.frame (even if only 1 redshift if requested) with the following columns:
\tabular{ll}{
z       \tab Requested redshift\cr
a       \tab Universe expansion factor, as given by a=1/(1+z)\cr
CoDist  \tab Line-of-sight (i.e. radial) comoving distance in units of Mpc\cr
LumDist \tab Luminosity distance in units of Mpc\cr
AngDist \tab Angular diameter distance in units of Mpc\cr
CoDistTran  \tab Transverse comoving distance in units of Mpc\cr
DistMod \tab The distance modulus used where AbsMag = ApMag - DistMod, and DistMod = 5log10(LumDist)+25 in units of mag\cr
AngSize \tab Physical projected size of an object at z in units of kpc/arcsec\cr
CoVol \tab Comoving volume of Universe within z in units of Gpc^3\cr
}
If age=TRUE is set then additional age-related information is calculated for each z as extra columns:
\tabular{ll}{
HubTime \tab Approximate Hubble age of the Universe in units of Gyrs\cr
UniAgeNow  \tab Age of the Universe now in units of Gyrs\cr
UniAgeAtz  \tab Age of the Universe at the specified redshift (z) in units of Gyrs\cr
TravelTime \tab Light travel time from the specified redshift (AKA look back time) in units of Gyrs\cr
}
If error=TRUE is set then the relative error for distance/age values is calculated for each z as an extra column:
\tabular{ll}{
RelError \tab Relative error of the distance/age integrals (this is the main source of error in the calculations)\cr
}

cosdistz function returns the input redshift (only included for clarity).

cosdista function returns the Universe expansion factor, as given by a=1/(1+z).

cosdistCoDist function returns the Line-of-sight (i.e. radial) comoving distance in units of Mpc.

cosdistLumDist function returns the Luminosity distance in units of Mpc.

cosdistAngDist function returns the Angular diameter distance in units of Mpc.

cosdistCoDistTran function returns the Transverse comoving distance in units of Mpc.

cosdistDistMod function returns the The distance modulus used where AbsMag = ApMag - DistMod, and DistMod = 5log10(LumDist)+25 in units of mag.

cosdistAngSize function returns the Physical projected size of an object at z in units of kpc/arcsec.

cosdistCoVol function returns the Comoving volume of Universe within z in units of Gpc^3.

cosdistHubTime function returns the Approximate Hubble age of the Universe in units of Gyrs.

cosdistUniAgeNow function returns the Age of the Universe now in units of Gyrs.

cosdistUniAgeAtz function returns the Age of the Universe at the specified redshift (z) in units of Gyrs.

cosdistTravelTime function returns the Light travel time from the specified redshift (AKA look back time) in units of Gyrs.
}
\references{
Based on the equations in Hogg D.W., 1999, arXiv, 9905116 and Wright E.L., 2006, PASP, 118, 1711
}
\author{
Aaron Robotham
}

\seealso{
\code{\link{cosvol}}, \code{\link{cosmap}}, \code{\link{cosgrow}}
}

\examples{
cosdist(0.3,70,age=TRUE)
cosdistz(0.3)
cosdista(0.3)
cosdistCoDist(0.3,70)
cosdistLumDist(0.3,70)
cosdistAngDist(0.3,70)
cosdistCoDistTran(0.3,70)
cosdistDistMod(0.3,70)
cosdistAngSize(0.3,70)
cosdistCoVol(0.3,70)
cosdistHubTime(70)
cosdistUniAgeNow(0.3,70)
cosdistUniAgeAtz(0.3,70)
cosdistTravelTime(0.3,70)

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

