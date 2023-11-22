\name{Sky Coordinate Matching}
\alias{coordmatch}
\alias{coordmatchsing}
\alias{internalclean}
\alias{group_links}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Sky matching
}
\description{
These functions allows the user to match a reference set of sky coordinates against a comparison set of sky coordinates. The match radius can be varied per source (all matches per source are given within this radius), and mutual best matches are also extracted. \code{coordmatch} should be used for finding multiple matches and \code{coordmatchsing} should be used when trying to find matches around a single source. \code{internalclean} is a utility function that will remove closely duplicated objects via some \option{tiebreak} criterion, and is probably only of interest to advanced users trying to clean catalogues that were produced from overlapping frames.
}
\usage{
coordmatch(coordref, coordcompare, rad = 2, inunitref = "deg", inunitcompare = "deg",
  radunit = "asec", sep = ":", kstart = 10, ignoreexact = FALSE, ignoreinternal = FALSE,
  matchextra = FALSE, smallapprox = FALSE, jitter = FALSE, jamount = 1e-12, jseed = 666)

coordmatchsing(RAref, Decref, coordcompare, rad = 2, inunitref = "deg",
  inunitcompare = "deg", radunit = 'asec', sep = ":", ignoreexact = FALSE,
  smallapprox = FALSE)

internalclean(RA, Dec, rad = 2, tiebreak, decreasing = FALSE, Nmatch = 'all',
  iter = FALSE, group = FALSE, ...)
  
group_links(links, selfgroup = FALSE, return_groupinfo = FALSE, return_linkslist = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{coordref}{
For coordmatch this is the reference dataset, i.e. you want to find matches for each object in this catalogue. A minimum two column matrix or data.frame, where column one is the RA and column two the Dec. See \option{matchextra}.
}
  \item{coordcompare}{
The comparison dataset, i.e. you want to find objects in this catalogue that match locations in coordref. A minimum two column matrix or data.frame, where column one is the RA and column two the Dec. If \option{coordcompare} is not provided then it is set to \option{coordref} automatically. Since this means the user is doing a single table internal match \option{ignoreinternal} is automatically set to TRUE (but this can be overridden). See \option{matchextra}.
}
  \item{RAref}{
For \code{coordmatchsing} this is the reference RA for the single object of interest.
}
  \item{Decref}{
For \code{coordmatchsing} this is the reference Dec for the single object of interest.
}
  \item{RA}{
For \code{internalclean} this is a vector of right ascensions for internal cleaning. If \option{RA} is a two column structure then the second column is taken to be \option{Dec}.
}
  \item{Dec}{
For \code{internalclean} this is a vector of declinations for internal cleaning. If \option{RA} is a two column structure then the second column is taken to be \option{Dec}. 
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
If inunitref, inunitcompare or inunit is set to 'sex' then sep defines the separation type as detailed in \code{\link{hms2deg}} and \code{\link{dms2deg}}.
}
  \item{kstart}{
The number of matching nodes to attempt initial. The code iterates until all matches within the specified radius (rad) have been found, but it works faster if the kstart is close to the maximum number of matches for any coordref object.
}
  \item{ignoreexact}{
Should exact matches be ignored in the output? If TRUE then 0 separation ID matches are set to 0 and the separation is NA. This might be helpful when matching the same table against itself, where you have no interest in finding object matches with respect to themselves.
}
  \item{ignoreinternal}{
Should identical row matches be ignored in the output? If TRUE then exact row ID matches are set to 0 and the separation is NA. The bestmatch output will ignore these trivial matchesw also. This only makes sense if \option{coordref} and \option{coordcompare} are the same table and you are trying to do an internal table match where you do not want the trivial result of rows matching to themselves. Automatically switches to TRUE if \option{coordcompare} is not provided. 
}
  \item{matchextra}{
Should extra columns in \option{coordref} and \option{coordcompare} be used as part of the N-D match? Extra columns beyond the requried RA and Dec can be provided and these will be used as part of the N-D match. The meaning of \option{rad} in this case is not trivial of course since the match is done within a hyper-sphere. When the extra columns have the same value \option{rad} can still be interpretted as an angular coordinate match. These extra columns should be appropriately scaled, e.g. you might want to make a 2 arcsec match with an extra magnitude column. In this case even if two objects sit on top of each other on sky, they cannot differ by more than 2 mag in flux to be a match.
}
  \item{smallapprox}{
Should the small angle approximation of asin(a/b) = a/b be used? If TRUE then some computations may be much faster, since asin is an expensive computation to make for lots of near matches.
}
  \item{jitter}{
Logical; should a very small jitter be applied to the \option{coordcompare} data? This should be negligibly small, but does make the output slightly variable depending on the seed. If \option{jitter} is not explicitly set and \option{ignoreinternal} is TRUE (which is the case if only one catalogue is supplied) and \option{ignoreexact} is FALSE then \option{jitter} is automatically switched to TRUE. See \code{\link{jitter}}.
}
  \item{jamount}{
Numeric scaler; max amount of jitter in radians (passed to \option{amount} in \code{\link{jitter}}). The default is less than a millionth of an asec, so there should be no circumstances where this affects real results (that is massively sub Event Horizon Telescope resolution even).
}
  \item{jseed}{
Integer scalar; random seed to use for jitter.
}
  \item{rad}{
The matching radius to use. If this is length one then the same radius is used for all objects, otherwise it must be the same length as the number of rows in coordref.
}
  \item{tiebreak}{
For \code{internalclean} this is a vector of values to determine the preferred source, e.g. something like magnitude of distance to the centre of the origin frame. By default smaller values are considered better, but this can be flipped by setting \option{decreasing}=TRUE. If \option{tiebreak} is not provided then the first source that appears is considered the better object in the cleaned catalogue.
}
  \item{decreasing}{
Determines whether smaller (\option{decreasing}=FALSE) or larger (\option{decreasing}=TRUE) \option{tiebreak} values are considered preferable.
}
  \item{Nmatch}{
Character scalar or integer vector; if default 'all' then cleaning is done for any number of matches and no match objects (unambiguous) are passed through. If a vector then the number of matches has to exist in the \option{Nmatch} vector to get passed out. This is useful if we only want cleaned objects in overlap regions where we expect a certain number of matches. Note \option{Nmatch}=1 means two objects match- the base reference and 1 other etc.
}
  \item{iter}{
Logical; should \code{internalclean} by run until no matches within the \option{rad} exist? The reason you might want to do this is to break degeneracies in chains where all links are within the search \option{rad} but the extremes are outside of this. In that situation you might want to keep sub-groups based on their preferred \option{tiebreak} (in which case use \option{iter} = FALSE) or you might want to fully resolve the whole complexes down to guarantee nothing matches within the \option{rad} (in which case use \option{iter} = TRUE). The latter is the only solution that guarantees no sources are left that match within \option{rad}. See Examples.
}
  \item{group}{
Logical; this goes a step further than \option{iter} and cretes full FoF groups. See Examples.
}
  \item{links}{
Integer Matrix; the two column link associations that we wish to group together. These do not need to be symmetric and can be non-contiguous, they must be integers though.
}
  \item{selfgroup}{
Logical; should 1-1 self-groupings be ignored (i.e. 2 matching to 2)?  
}
  \item{return_groupinfo}{
Logical; should basic group info (Ngroup, Nlist, IDmin, IDmax) be returned as a data.frame? Note Nlist is only present if \option{return_linkslist} is also TRUE.
}
  \item{return_linkslist}{
Logical; should the links that form the groups be returned?
}
  \item{\dots}{
Other arguments to be passed to \code{coordmatch}.  
}
}
\details{
For coordmatch the main matching is done using nn2 that comes as part of the RANN package. coordmatch adds a large amount of sky coordinate oriented functionality beyond the simple implementation of nn2. For single object matches coordmatchsing should be used since it is substantially faster in this regime (making use of direct dot products).

\option{ignoreexact} is more strict in a sense since all objects exactly matching are ignored, whereas with \option{ignoreinternal} only identical row IDs are interpretted as being the same object. This setting can be useful in \code{internalclean} (via \dots) too, since having exact matches and internal matches can create ambiguous cleaning.
}
\value{
The output of \code{coordmatch} is a list containing:

\item{ID}{The full matrix of matching IDs. The rows are ordered identically to \option{coordref}, and the ID value is the row position in \option{coordcompare} for the match.}
\item{sep}{The full matrix of matching separations in the same units as \option{radunit}. The rows are ordered identically to \option{coordref}, and the sep value is the separation for each matrix location in the ID list object.}
\item{Nmatch}{Nmatch is a vector giving the total number of matches for each \option{coordref} row.}
\item{bestmatch}{A three column data.frame giving the best matching IDs. Only objects with at least one match are listed. Column 1 (refID) gives the row position from \option{coordref} and column 2 (compareID) gives the corresponding best matching row position in \option{coordcompare}. Column 3 (sep) gives the separation between the matched ref and compare positions in the same units as radunit.}

The output of \code{coordmatchsing} is a list containing:

\item{ID}{The full vector of matching IDs. The ID values are the row positions in \option{coordcompare} for the match.}
\item{sep}{The full vector of matching separations in the same units as \option{radunit}. The sep value is the separation for each vector location in the ID list object.}
\item{Nmatch}{Total number of matches within the specified radius.}
\item{bestmatch}{The best matching ID, where the ID value is the row position in \option{coordcompare} for the match.}

The output of \code{internalclean} is a vector containing IDs of the rows to keep to achieve the cleaned catalogue.

The output of \code{group_links} depends on the return parameters. Minimally it is a list of the FoF groups. Can also be a list structure containing:

\item{grouplist}{The list of the FoF groups.}
\item{groupinfo}{The table of basic group information when \option{return_groupinfo} is TRUE. The data.frame will contain Ngroup, Nlist, IDmin, IDmax, though note Nlist is only present if \option{return_linkslist} is also TRUE.}
\item{linkslist}{The list of the links that form the respective groups.}
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

mockmatchesvary=coordmatch(mocksky, mocksky, ignoreexact=TRUE, rad=seq(0,1,length=1e3),
radunit='amin')

#We can do this also by using the internal table match mode:

mockmatchesvary2=coordmatch(mocksky, rad=seq(0,1,length=1e3), radunit='amin')

#Check that this looks the same (should be identical with all zeroes):

summary(mockmatchesvary$bestmatch-mockmatchesvary2$bestmatch)

#Internal matching can be complicated:

coords = cbind(c(-8:8, 11:14),c(-8:8, 11:14))*0.0001
plot(coords)

#Note the search radius is such two adjacent objects can be associated:

library(plotrix)
points(0,0,pch=4)
draw.circle(0,0,radius=1.2/3600)

(internalclean(coords, rad=1.2))
(internalclean(coords, rad=1.2, iter=TRUE))
(internalclean(coords, rad=1.2, group=TRUE))

#Basic Fof groups:
set.seed(666)
links = cbind(sample(40,20,replace=TRUE), sample(40,20,replace=TRUE))
print(links)

group_links(links, return_groupinfo=TRUE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{match}
\keyword{matching}
\keyword{coord}
\keyword{coordinates}
