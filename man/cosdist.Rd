\name{cosdist}
\alias{cosdist}
\alias{cosdistz}
\alias{cosdista}
\alias{cosdistCoDist}
\alias{cosdistLumDist}
\alias{cosdistAngDist}
\alias{cosdistAngDist12}
\alias{cosdistCoDistTran}
\alias{cosdistCoDist12ang}
\alias{cosdistDistMod}
\alias{cosdistAngScale}
\alias{cosdistAngSize}
\alias{cosdistAngArea}
\alias{cosdistCoVol}
\alias{cosdistHubTime}
\alias{cosdistUniAgeNow}
\alias{cosdistUniAgeAtz}
\alias{cosdistTravelTime}
\alias{cosdistRelError}
\alias{cosdistCrit}
\alias{cosdistCoDist12ang}
\alias{cosdistLumDist12ang}
\alias{cosdistAngDist12ang}
\alias{cosdistzem12ang}
\alias{cosdistzeff12ang}
\alias{cosdistzeff}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Cosmological distance calculations
}
\description{
These functions allow comoving, angular size and luminosity distances to be calculated for a given redshift, it can also return look back time. They use curvature correctly, calculated internally using the relation OmegaM+OmegaL+OmegaR+OmegaK=1, but by default they assume a flat Universe where only OmegaM needs to be specified and OmegaR=0 (so no radiation pressure at any epoch).
}
\usage{
cosdist(z=1, H0=100, OmegaM=0.3, OmegaL=1-OmegaM-OmegaR, OmegaR=0, w0 = -1, wprime = 0,
age=FALSE, ref, error=FALSE)
cosdistz(z=1)
cosdistzeff(zref = 1, zem = 2)
cosdista(z=1)
cosdistCoDist(z=1, H0=100, OmegaM=0.3, OmegaL=1-OmegaM-OmegaR, OmegaR=0, w0 = -1,
wprime = 0, ref)
cosdistLumDist(z=1, H0=100, OmegaM=0.3, OmegaL=1-OmegaM-OmegaR, OmegaR=0, w0 = -1,
wprime = 0, ref)
cosdistAngDist(z=1, H0=100, OmegaM=0.3, OmegaL=1-OmegaM-OmegaR, OmegaR=0, w0 = -1,
wprime = 0, ref)
cosdistAngDist12(z1=1, z2=2, H0=100, OmegaM=0.3, OmegaL=1-OmegaM-OmegaR, OmegaR=0,
w0 = -1, wprime = 0, ref)
cosdistCoDistTran(z=1, H0=100, OmegaM=0.3, OmegaL=1-OmegaM-OmegaR, OmegaR=0, w0 = -1,
wprime = 0, ref)
cosdistCoDist12ang(z1=1, z2=2, ang=0, H0=100, OmegaM=0.3, OmegaL=1-OmegaM-OmegaR,
OmegaR=0, w0 = -1, wprime = 0, inunit='deg', ref) 
cosdistLumDist12ang(z1=1, z2=2, ang=0, H0=100, OmegaM=0.3, OmegaL=1-OmegaM-OmegaR,
OmegaR=0, w0 = -1, wprime = 0, inunit='deg', ref)
cosdistAngDist12ang(z1=1, z2=2, ang=0, H0=100, OmegaM=0.3, OmegaL=1-OmegaM-OmegaR,
OmegaR=0, w0 = -1, wprime = 0, inunit='deg', ref)
cosdistzem12ang(z1=1, z2=2, ang=0, H0=100, OmegaM=0.3, OmegaL=1-OmegaM-OmegaR,
OmegaR=0, w0 = -1, wprime = 0, inunit='deg', ref)
cosdistzeff12ang(z1=1, z2=2, ang=0, H0=100, OmegaM=0.3, OmegaL=1-OmegaM-OmegaR,
OmegaR=0, w0 = -1, wprime = 0, inunit='deg', ref)
cosdistDistMod(z=1, H0=100, OmegaM=0.3, OmegaL=1-OmegaM-OmegaR, OmegaR=0, w0 = -1,
wprime = 0, ref)
cosdistAngScale(z=1, H0=100, OmegaM=0.3, OmegaL=1-OmegaM-OmegaR, OmegaR=0, w0 = -1,
wprime = 0, ref)
cosdistAngSize(z=1, Size=1, H0=100, OmegaM=0.3, OmegaL=1-OmegaM-OmegaR, OmegaR=0, w0=-1,
wprime=0, Dim=1, Dist='Co',  outunit='deg', ref)
cosdistAngArea(z=1, Size=1, H0=100, OmegaM=0.3, OmegaL=1-OmegaM-OmegaR, OmegaR=0, w0=-1,
wprime=0, Dim=2, Dist='Co',  outunit='deg2', ref)
cosdistCoVol(z=1, H0=100, OmegaM=0.3, OmegaL=1-OmegaM-OmegaR, OmegaR=0, w0 = -1,
wprime = 0, ref)
cosdistHubTime(H0=100)
cosdistUniAgeNow(z=1, H0=100, OmegaM=0.3, OmegaL=1-OmegaM-OmegaR, OmegaR=0, w0 = -1,
wprime = 0, ref)
cosdistUniAgeAtz(z=1, H0=100, OmegaM=0.3, OmegaL=1-OmegaM-OmegaR, OmegaR=0, w0 = -1,
wprime = 0, ref)
cosdistTravelTime(z=1, H0=100, OmegaM=0.3, OmegaL=1-OmegaM-OmegaR, OmegaR=0, w0 = -1,
wprime = 0, ref)
cosdistRelError(z=1, OmegaM=0.3, OmegaL=1-OmegaM-OmegaR, OmegaR=0, w0 = -1, wprime = 0,
ref)
cosdistCrit(z_lens=1, z_source=2, H0=100, OmegaM=0.3, OmegaL=1-OmegaM-OmegaR, OmegaR=0,
w0 = -1, wprime = 0, ref)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{z}{
    Cosmological redshift, where z must be > -1 (can be a vector).
  }
  \item{H0}{
    Hubble constant as defined at z=0 (default is H0=100 (km/s)/Mpc).
  }
  \item{OmegaM}{
    Omega matter (default is 0.3).
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
  \item{age}{
    Flag for cosdist function to return age or not- this slows calculation, so is by default turned off.
  }
  \item{ref}{
    The name of a reference cosmology to use, one of 137 / 737 / Planck / Planck13 / Planck15 / WMAP / WMAP9 / WMAP7 / WMAP5 / WMAP3 / WMAP1 / Millennium / GiggleZ. Planck=Planck15 and WMAP=WMAP9. The usage is case insensitive, so wmap9 is an allowed input. This overrides any other settings for H0, OmegaM and OmegaL. If OmegaR is missing from the reference set then it is inherited from the function input (0 by default). See \code{\link{cosref}} for details.
  }
  \item{error}{
    Flag for cosdist to calculate the relative error for distance/age values.
  }
  \item{z1}{
    Redshift for object 1, where z1 must be > -1 (can be a vector) and less than z2.
  }
  \item{z2}{
    Redshift for object 2, where z2 must be > -1 (can be a vector) and greater than z1.
  }
  \item{zref}{
    Redshift for the reference object, i.e. the object that we caste as the observer of another object at zem.
  }
  \item{zem}{
    Redshift for the emitting  object, i.e. the object that we caste as being observed by another object at zref.
  }
  \item{z_lens}{
    Redshift, where z_lens must be > -1 (can be a vector) and z_lens < z_source.
  }
  \item{z_source}{
    Redshift, where z_source must be > -1 (can be a vector) and z_lens < z_source.
  }
  \item{ang}{
    The observed angular separation between object 1 and object 2 in degrees.
  }
  \item{Size}{
    The 1D size of the object (i.e. diameter or total length) in Mpc. Either comoving or angular, as specified by \option{Dist}. For cosdistAngArea this is always taken to be the diameter of either projected or 3D object.
  }
  \item{Dim}{
    Specifies whether the object being considered is 1D (a line) 2D (e.g. a face on galaxy) or 3D (e.g. dark matter halo). This makes a very small modification to the geometry used (tan of 1D/2D and sin for 3D), but is only noticeable for large structures at low redshifts.
  }
  \item{Dist}{
    Determines the distance type, i.e. angular / physical distances (Ang) or with respect to comoving distances (Co).
  }
  \item{inunit}{
    The units of angular coordinate provided for \option{ang}. Allowed options are deg for degrees, amin for arc minutes, asec for arc seconds, and rad for radians.
  }
  \item{outunit}{
    For \code{cosdistAngSize} units of angular coordinate output. Allowed options are deg for degress (default), amin for arc minutes, asec for arc seconds, and rad for radians.
    
    For \code{cosdistAngArea} units of angular area output. Allowed options are deg2 for square degrees (default), amin2 for square arc minutes, asec2 for square arc seconds and rad2 or sr for steradians.
  }
}
\details{
Functions are largely based on D. W. Hogg et al. 1999 and Wright et al. 2006.

Negative value of z> -1 are allowed, which produces future predictions based on present day cosmology.

cosdistAngDist12 is only available for OmegaK>=0.
}
\value{
\item{cosdist}{Returns a data.frame (even if only 1 redshift if requested) with the following columns:}
\tabular{ll}{
z       \tab Requested redshift\cr
a       \tab Universe expansion factor, as given by a=1/(1+z)\cr
CoDist  \tab Line-of-sight (i.e. radial) comoving distance in units of Mpc\cr
LumDist \tab Luminosity distance in units of Mpc\cr
AngDist \tab Angular diameter distance in units of Mpc\cr
CoDistTran  \tab Transverse comoving distance in units of Mpc\cr
DistMod \tab The distance modulus used where AbsMag = ApMag - DistMod, and DistMod = 5log10(LumDist)+25 in units of mag\cr
AngScale \tab Physical projected scale of an object at z in units of kpc/arcsec\cr
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

The outputs of the standalone functions are:

\item{cosdistz}{Returns the input redshift (only included for clarity).}

\item{cosdistzeff}{Returns the apparent redshift that the object at zref will observe the object at zem for the Universe age that zref is observed to have now. This is given by \eqn{(1+zem)/(1+zref)}{(1+zem)/(1+zref)}.}

\item{cosdista}{Returns the Universe expansion factor, as given by \eqn{a=1/(1+z)}{a=1/(1+z)}.}

\item{cosdistCoDist}{Returns the line-of-sight (i.e. radial) comoving distance in units of Mpc. For a flat Universe (OmegaK=0) this is exactly the samething as the transverse comoving distance, and by extension it is also the proper motion distance.}

\item{cosdistLumDist}{Returns the luminosity distance in units of Mpc.}

\item{cosdistAngDist}{Returns the angular diameter distance in units of Mpc.}

\item{cosdistAngDist12}{Returns the radial angular diameter distance separation in units of Mpc between objects at \option{z1} and \option{z2} that have small angular separations on sky.}

\item{cosdistCoDistTran}{Returns the transverse comoving distance in units of Mpc. This is equivilant to the proper motion distance for all values of Universe curvature (OmegaK !=0), and is the same thing as the line-of-sight comoving distance for a flat Universe (OmegaK=0).}

\item{cosdistCoDist12ang}{Returns the total comoving distance in units of Mpc between objects at \option{z1} and \option{z2} with a separation \option{ang}. This works for curved cosmologies (i.e. OmegaK!=0) and for large radial and tangential separations. For small separations at a certain value of z for both objects the result is very similar to cosdistCoDistTran(z)*sin(ang*pi/180). This function was mostly extracted from Eqn 3.19 in Peacock (1999).}

\item{cosdistLumDist12ang}{Returns the total luminosity distance in units of Mpc between objects at \option{z1} and \option{z2} with a separation \option{ang}. This is equal to cosdistCoDist12ang*(1+zeff), where zeff is the apparent redshift that the object at z1 will observe the object at z2 for the Universe age that z1 is observed to have now. See \code{cosdistCoDist12ang} for details.}

\item{cosdistAngDist12ang}{Returns the total angular diameter distance in units of Mpc between objects at \option{z1} and \option{z2} with a separation \option{ang}. This is equal to cosdistCoDist12ang/(1+zeff), where zeff is the apparent redshift that the object at z1 will observe the object at z2 for the Universe age that z1 is observed to have now. See \code{cosdistCoDist12ang} for details.}

\item{cosdistzem12ang}{Returns the apparent redshift that the object at \option{z1} would observe the object at \option{z2} to be for our current Universe age. See \code{cosdistCoDist12ang} for details.}

\item{cosdistzeff12ang}{Returns the apparent redshift that the object at \option{z1} would observe the object at \option{z2} to be for the Universe age that z1 is observed to have now. See \code{cosdistCoDist12ang} for details.}

\item{cosdistDistMod}{Returns the distance modulus used where AbsMag = ApMag - DistMod, and DistMod = 5log10(LumDist)+25 in units of mag.}

\item{cosdistAngScale}{Returns the physical projected scale of an object at z in units of kpc/arcsec.}

\item{cosdistAngSize}{Returns the angular size (length or diameter) of an object (by default in degrees).}

\item{cosdistAngArea}{Returns the angular area of an object (by default degrees^2), taking the specified Size to be the diameter.}

\item{cosdistCoVol}{Returns the comoving volume of Universe within z in units of Gpc^3.}

\item{cosdistHubTime}{Returns the approximate Hubble age of the Universe in units of Gyrs.}

\item{cosdistUniAgeNow}{Returns the age of the Universe now in units of Gyrs.}

\item{cosdistUniAgeAtz}{Returns the age of the Universe at the specified redshift (z) in units of Gyrs.}

\item{cosdistTravelTime}{Returns the light travel time from the specified redshift (AKA look back time) in units of Gyrs.}

\item{cosdistRelError}{Returns the relative error of the distance/age integrals (this is the main source of error in the calculations).}

\item{cosdistCrit}{Returns the critical surface mass density, SigmaC (see also \code{\link{cosNFW}}).}
}
\references{
Based on the equations in:

Davis T.M. & Lineweaver, Charles H., 2004, PASA, 21, 97
  
Hogg D.W., 1999, arXiv, 9905116

Liske J., 2000, MNRAS, 319, 557L

Peacock J.A., 1999, Cosmological Physics, Cambridge University Press

Wright E.L., 2006, PASP, 118, 1711
}
\author{
Aaron Robotham
}

\seealso{
\code{\link{cosvol}}, \code{\link{cosmap}}, \code{\link{cosgrow}}, \code{\link{cosref}}, \code{\link{cosNFW}}
}

\examples{
\dontrun{
cosdist(0.3,70,age=TRUE)
cosdist(0.3,70,age=TRUE,ref='Planck')
cosdistz(0.3)
cosdista(0.3)
cosdistCoDist(0.3,70)
cosdistLumDist(0.3,70)
cosdistAngDist(0.3,70)
cosdistAngDist12(0.3,0.5,70)
cosdistCoDistTran(0.3,70)
cosdistCoDist12ang(0,2,10)
cosdistDistMod(0.3,70)
cosdistAngScale(0.3,70)
cosdistAngSize(0.3,1,70)
cosdistCoVol(0.3,70)
cosdistHubTime(70)
cosdistUniAgeNow(0.3,70)
cosdistUniAgeAtz(0.3,70)
cosdistTravelTime(0.3,70)
cosdistRelError(0.3)
cosdistCrit(0.3,0.5,70)
cosdistzeff(1,2)
cosdistzem12ang(1,2)
cosdistzeff12ang(1,2)

#A check of the comoving separation between objects function:

cosdistCoDistTran(2,OmegaM = 0.3, OmegaL=1)*sin(pi/180)
cosdistCoDist12ang(2,2,ang=1,OmegaM=0.3,OmegaL=1)

#Very close, however cosdistCoDist12ang lets us go further:

cosdistCoDist12ang(1,2,ang=10,OmegaM=0.3,OmegaL=1)
cosdistCoDist12ang(2,2,ang=180,OmegaM=0.3,OmegaL=1)

#The second number should be be the same as:

cosdistCoDist(2,OmegaM=0.3,OmegaL=1)*2

#Example 1 by John Peacock for EDS Universe (answer should be nearly 3):

cosdistzem12ang(3,4,56.4,H0=100,OmegaM=1,OmegaL=0)

#Example 2 by John Peacock for EDS Universe (answer should be nearly 2995 Mpc/h):

cosdistCoDist12ang(3,4,56.4,H0=100,OmegaM=1,OmegaL=0)

#Example 3 by John Peacock for Milne Universe (answer should be nearly 5294 Mpc/h):

cosdistCoDist12ang(3,4,56,H0=100,OmegaM=0,OmegaL=0)

#Example 4 by John Peacock for Milne Universe (answer should be nearly 4.846):

cosdistzeff12ang(3,4,56,H0=100,OmegaM=0,OmegaL=0)

#Example 5 by John Peacock for Milne Universe (answer should be nearly 364 Mpc/h):

cosdistAngDist12ang(3,4,56,H0=100,OmegaM=0,OmegaL=0)

#Nice plot of distance estimates:

redshifts=seq(0,3,by=0.01)
plot(redshifts, cosdistCoDist(redshifts, ref='planck'), type='l', col='darkgreen',
xlab='Redshift / z', ylab='Distance / Mpc')
lines(redshifts, cosdistLumDist(redshifts, ref='planck'), col='red')
lines(redshifts, cosdistAngDist(redshifts, ref='planck'), col='blue')
legend('topleft', legend=c('Comoving Distance', 'Luminosity Distance', 'Angular Diameter Distance'),
col=c('darkgreen', 'red', 'blue'),lty=1)

plot(redshifts, cosdistTravelTime(redshifts, ref='planck'), type='l',
xlab='Redshift / z', ylab='Light travel time / Yrs')

#Actual time example (Figure 1 of Davis & Lineweaver 2004)
zseq=10^seq(-2,6,len=1e3)-1
dists=cosdistCoDist(zseq, ref='737')*0.00326
times=cosdistTravelTime(zseq, ref='737')
plot(dists, times, type='l', xlab='Comoving Distance / Glyr',
ylab='Time / Gyr')
abline(v=0, h=0, lty=1)
abline(h=c(min(times), max(times)), lty=2)
abline(v=c(min(dists), max(dists)), lty=2)

#Conformal time example (Figure 1 of Davis & Lineweaver 2004):
#Mpc to Glyr conversion is 0.00326

zseq=10^seq(-2,6,len=1e3)-1
dists=cosdistCoDist(zseq, ref='737')*0.00326
plot(dists, dists, type='l',
xlab='Comoving Distance / Glyr', ylab='Conformal Time / Gyr')
abline(v=0, h=0, lty=1)
abline(h=c(min(dists), max(dists)), lty=2)
abline(v=c(min(dists), max(dists)), lty=2)
}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{cosmo}
\keyword{cosmology}
\keyword{distance}% __ONLY ONE__ keyword per line

