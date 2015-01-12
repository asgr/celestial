\name{cosmap}
\alias{cosmap}
\alias{cosmapval}
\alias{cosmapfunc}
\title{
Cosmological Mapping Functions
}
\description{
Functions for mapping from one arbitrary csomological parameter to another. This includes the provision of a generic interpolation function and another exact value lookup.
}
\usage{
cosmapval(val = 1, cosparamx = "z", cosparamy = "CoDist", H0 = 100, OmegaM = 0.3, OmegaL = 1 - OmegaM, zrange = c(0, 20), step='z', res = 100)
cosmapfunc(cosparamx = "z", cosparamy = "CoDist", H0 = 100, OmegaM = 0.3, OmegaL = 1 - OmegaM, zrange = c(0, 20), step='z', res = 100)
}

\arguments{
  \item{val}{
The value/s to be mapped from parameter cosparamx to parameter cosparamy (this can be a vector or a single number).
}
  \item{cosparamx}{
Cosmological parameter, must be one of: z, a, CoDist, LumDist, CoDistTran, DistMod, CoVol, UniAgeAtz, TravelTime (see cosdist help for further description of these). Note that AngDist and AngSize are not an option for cosparamx due to degenerate redshift solutions causing problems approxfun mapping.
}
  \item{cosparamy}{
Cosmological parameter, must be one of: z, a, CoDist, LumDist, AngDist, CoDistTran, DistMod, AngSize, CoVol, UniAgeAtz, TravelTime (see cosdist help for further description of these).
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
  \item{zrange}{
Lower and upper z limits that the approxfun mapping is generated over (increase range if default is not sufficient, and decrease if it is wasteful, i.e. the possible redshift window is known to be quite narrow).
}
  \item{step}{
The type of stepping used. Allowed values are 'z' (uniform stepping in z), 'logz' (uniform stepping in log10(1+z) and expansion factor 'a' (uniform stepping in a=1/(1+z)). Default is z. For mappings using time (UniAgeNow, UniAgeAtz, TravelTime) or comoving quantities (CoDist, CoDistTran, CoVol) or distance modulus (DistMod) 'a' or 'logz' map the numeric range more uniformly. This is because a and log10(1+z) are approximately linear in light travel time (positive and negative correlation respectively), and typically they have better behaviour than stepping uniformly in z directly.
}
  \item{res}{
The resolution of steps. Larger numbers will be more accurate, but will be slower to compute.
}
}
\details{
The default zrange and res should be sufficient for most reasonable cosmologies if the approximate redshift location of the region to be mapped is entirely unknown.
}
\value{
cosmapval contains the value/s at cosparamx mapped onto cosparamy, in the natural units of cosparamx: unitless (z), unitless (a), Mpc (CoDist), Mpc (LumDist), Mpc (AngDist), Mpc (CoDistTran), mag (DistMod), kpc/arcsec (AngSize), Gpc^3 (CoVol), Gyrs (UniAgeAtz), Gyrs (TravelTime)

cosmapfunc uses base R approxfun to map cosparamx onto cosparamy between zrange[1] and zrange[2] in uniform steps of expansion factor (a=1/(1+z)). cosmofunc returns the output function created by approxfun.
}
\references{
Based on the equations in Hogg D.W., 1999, arXiv, 9905116 and Wright E.L., 2006, PASP, 118, 1711
}
\author{
Aaron Robotham
}
\seealso{
\code{\link{cosdist}}, \code{\link{cosvol}}
}

\examples{
tempfunc=cosmapfunc('CoVol', 'UniAgeAtz')
tempfunc(50)

cosmapval(50, 'CoVol', 'UniAgeAtz')
}

\keyword{ cosmo }
\keyword{ cosmology }
\keyword{ mapping }
