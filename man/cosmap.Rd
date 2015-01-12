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
cosmapval(val = 50, cosparam = "CoVol", H0 = 100, OmegaM = 0.3, OmegaL = 1 - OmegaM, zrange = c(0, 100), res = 10, iter = 12 ,age = FALSE)
cosmapfunc(cosparamx = "CoVol", cosparamy = "z", H0 = 100, OmegaM = 0.3, OmegaL = 1 - OmegaM, zrange = c(0, 20), step='z', res = 100)
}

\arguments{
  \item{val}{
The value/s to be mapped from parameter cosparamx to parameter cosparamy (this can be a vector or a single number).
}
  \item{cosparam}{
Cosmological parameter, must be one of: z, a, CoDist, LumDist, CoDistTran, DistMod, CoVol, UniAgeAtz, TravelTime (see cosdist help for further description of these). Note that AngDist and AngSize are not an option for cosparamx due to degenerate redshift solutions causing problems approxfun mapping.
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
  \item{iter}{
The number of iterations to make when calculating the exact location of a given cosmological parameter when using cosmapval.  
}
  \item{age}{
Flag to tell cosdist function to return age or not- this slows calculation, so is by default turned off.
}
}
\details{
The default zrange and res should be sufficient for most reasonable cosmologies if the approximate redshift location of the region to be mapped is entirely unknown.

The default res and iter for cosmapval is appropriate for most mappings with z<20 using a fiducial 737 cosmology. If this proves insufficient (this should be obvious from error column) then increase both of these. Overall accuracy goes as res^iter.
}
\value{
cosmapval contains the cosdist output for parameter 'cosparam' at value 'val'. See \code{\link{cosdist}} for information on the cosdist output. It includes an additional final column named 'error' which gives the approximate relative error of the values returned compared to the desired lookup location. Smaller is obviously better, but at the cost of computational time.

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

cosmapval(50:60, 'CoVol')
}

\keyword{ cosmo }
\keyword{ cosmology }
\keyword{ mapping }
