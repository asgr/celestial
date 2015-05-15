\name{Sky Coordinate Matching}
\alias{coordmatch}
\alias{coordmatchsing}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Sky matching
}
\description{
These functions allows the user to match a reference set of sky coordinates against a comparison set of sky coordinates. The match radius can be varied per source (all matches per source are given within this radius), and mutual best matches are also extracted. coordmatch should be used for finding multiple matches and coordmatchsing should be used when trying to find matches around a single source.
}
\usage{
coordmatch(coordref, coordcompare, rad = 2, inunitref = "deg", inunitcompare = "deg", radunit = "asec", sep = ":", kstart = 10, ignoreexact = FALSE)
coordmatchsing(RAref,Decref, coordcompare, rad=2, inunitref = "deg", inunitcompare="deg", radunit='asec', sep = ":", ignoreexact=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{coordref}{
For coordmatch this is the reference dataset, i.e. you want to find matches for each object in this catalogue. A two column matrix or data.frame, where column one is the RA and column two the Dec.
}
  \item{RAref}{
For coordmatch sing this is the reference RA for the sinlge object of interest.
}
  \item{Decref}{
For coordmatch sing this is the reference Dec for the sinlge object of interest.
}
  \item{coordcompare}{
The comparison dataset, i.e. you want to find objects in this catalogue that match locations in coordref. A two column matrix or data.frame, where column one is the RA and column two the Dec.
}
  \item{rad}{
The matching radius to use. If this is length one then the same radius is used for all objects, otherwise it must be the same length as the number of rows in coordref.
}
  \item{inunitref}{
The units of angular coordinate provided for coordref / RAref / Decref. Allowed options are deg for degress, rad for radians and sex for sexigesimal (i.e. HMS for RA and DMS for Deg).
}
  \item{inunitcompare}{
The units of angular coordinate provided for coordcompare. Allowed options are deg for degress, rad for radians and sex for sexigesimal (i.e. HMS for RA and DMS for Deg).
}
  \item{radunit}{
The unit type for the radius specified. Allowed options are deg for degress, amin for arc minutes, asec for arc seconds and rad for radians.
}
  \item{sep}{
If inunitref or inunitcompare is set to 'sex' then sep defines the separation type as detailed in \code{\link{hms2deg}} and \code{\link{dms2deg}}.
}
  \item{kstart}{
The number of matching nodes to attempt initial. The code iterates until all matches within the specified radius (rad) have been found, but it works faster if the kstart is close to the maximum number of matches for any coordref object.
}
  \item{ignoreexact}{
Shoudl exact matches be ignored in the output? This might be helpful when matching the same table against itself, where you have no interest in finding object matches with respect to themselves.
}
}
\details{
For coordmatch the main matching is done using nn2 that comes as part of the RANN package. coordmatch adds a large amount of sky coordinate oriented functionality beyond the simple implementation of nn2. For single object matches coordmatchsing should be used since it is substantially faster in this regime (making use of direct dot products).
}
\value{
The output of coordmatch is a list containing:

\item{ID}{The full matrix of matching IDs. The rows are ordered identically to coordmatch, and the ID value is the row position in coordcompare for the match.}
\item{sep}{The full matrix of matching separations in the same units as radunit. The rows are ordered identically to coordmatch, and the sep value is the separation for each matrix location in the ID list object.}
\item{Nmatch}{Nmatch is a vector giving the total number of matches for each coordmatch row.}
\item{bestmatch}{A two column matrix giving the best matching IDs. Only objects with at least one match are listed. Column 1 gives the row position from coordref and column 2 gives the corresponding best matching position in coordcompare.}

The output of coordmatchsing is a list containing:

\item{ID}{The full vector of matching IDs. The ID values are the row positions in coordcompare for the match.}
\item{sep}{The full vector of matching separations in the same units as radunit. The sep value is the separation for each vector location in the ID list object.}
\item{Nmatch}{Total number of matches within the specified radius.}
\item{bestmatch}{The best matching ID, where the ID value is the row position in coordcompare for the match.}

}
\author{
Aaron Robotham
}

\seealso{
\code{\link{hms2deg}}, \code{\link{dms2deg}}, \code{\link{sph2car}}
}

\examples{
set.seed(666)

#Here we make objects in a virtual 1 square degree region

mocksky=cbind(runif(1e3), runif(1e3))

#Now we match to find all objects within an arc minute, ignoring self matches

mockmatches=coordmatch(mocksky, mocksky, ignoreexact=TRUE, rad=1, radunit='amin')

#Now we match to find all objects with varying match radii, ignoring self matches

mockmatchesvary=coordmatch(mocksky, mocksky,ignoreexact=TRUE, rad=seq(0,1,length=1e3), radunit='amin')
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{match}
\keyword{matching}
\keyword{coord}
\keyword{coordinates}
