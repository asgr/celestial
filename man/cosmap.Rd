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
cosmapval(val=50, cosparam="CoVol", H0=100, OmegaM=0.3, OmegaL=1-OmegaM, zrange=c(0,100),
res=10, iter=12, out='cos')
cosmapfunc(cosparamx="CoVol", cosparamy="z", H0=100, OmegaM=0.3, OmegaL=1-OmegaM,
zrange=c(0,20), step='z', res=100)
}

\arguments{
  \item{val}{
The value/s to be mapped from parameter cosparamx to parameter cosparamy (this can be a vector or a single number).
}
  \item{cosparam}{
Cosmological parameter, must be one of: z, a, CoDist, LumDist, CoDistTran, DistMod, CoVol, UniAgeAtz, TravelTime (see \code{\link{cosdist}} help for further description of these) H, OmegaM, OmegaL, OmegaK, Factor, Rate, RhoCrit (see \code{\link{cosgrow}} help for further description of these). Note that AngDist and AngSize are not an option for cosparam due to degenerate redshift solutions causing problems with the approxfun mapping.
}
  \item{cosparamx}{
Cosmological parameter, must be one of: z, a, CoDist, LumDist, CoDistTran, DistMod, CoVol, UniAgeAtz, TravelTime (see \code{\link{cosdist}} help for further description of these) H, OmegaM, OmegaL, OmegaK, Factor, Rate, RhoCrit (see \code{\link{cosgrow}} help for further description of these). Note that AngDist and AngSize are not an option for cosparamx due to degenerate redshift solutions causing problems with the approxfun mapping.
}
  \item{cosparamy}{
Cosmological parameter, must be one of: z, a, CoDist, LumDist, AngDist, CoDistTran, DistMod, AngSize CoVol, UniAgeAtz, TravelTime (see \code{\link{cosdist}} help for further description of these) H, OmegaM, OmegaL, OmegaK, Factor, Rate, RhoCrit (see \code{\link{cosgrow}} help for further description of these).
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
  \item{out}{
Either out='cos', in which case the output is a data.frame containing the output of \code{\link{cosdist}} and \code{\link{cosgrow}} for the specified input vaue, or out='z', in which case the output is a vector of the corresponding redshift (z) values.
}
}
\details{
The default zrange and res should be sufficient for most reasonable cosmologies if the approximate redshift location of the region to be mapped is entirely unknown.

The default res and iter for cosmapval is appropriate for most mappings with z<20 using a fiducial 737 cosmology. If this proves insufficient (this should be obvious from error column) then increase both of these. Overall accuracy goes as res^iter.
}
\value{
If out='cos', cosmapval contains the concatenation of the cosdist (with age=TRUE and error=TRUE) and cosgrow functions for parameter 'cosparam' at value 'val'. The 'z' and 'a' columns are only included once (from the output of cosdist). See \code{\link{cosdist}} and \code{\link{cosgrow}} for information on the cosdist and cosgrow outputs. If out='z', then cosmapval merely returns the corresponding redshifts.

The cosmapval output (when out='cos') includes an additional final column named 'MapError' which gives the approximate relative error of the values returned compared to the desired lookup location. Smaller is obviously better, but at the cost of computational time.

cosmapfunc uses base R approxfun to map cosparamx onto cosparamy between zrange[1] and zrange[2] in uniform steps of expansion factor (a=1/(1+z)). cosmofunc returns the output function created by approxfun.
}
\references{
Based on the equations in Hogg D.W., 1999, arXiv, 9905116 and Wright E.L., 2006, PASP, 118, 1711
}
\author{
Aaron Robotham
}
\seealso{
\code{\link{cosdist}}, \code{\link{cosvol}}, \code{\link{cosgrow}}
}

\examples{
tempfunc=cosmapfunc('CoVol', 'UniAgeAtz')
tempfunc(50)

cosmapval(50:60, 'CoVol')
}

\keyword{ cosmo }
\keyword{ cosmology }
\keyword{ mapping }
