\name{car2sph}
\alias{car2sph}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Transforms 3D cartesian coordinates to spherical coordinates
}
\description{
Transforms 3D cartesian coordinates to spherical coordinates. The user can choose to return the spherical coordinates in degrees or radians.
}
\usage{
car2sph(x, y, z, deg = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
x values, can also contain a matrix of x, y and z (in that order).
}
  \item{y}{
y values.
}
  \item{z}{
z values
}
  \item{deg}{
Should degrees be returned (default) or radians.
}
}
\details{
This is a low level function that is used for plot transformations.
}
\value{
A data.frame is returned containing the columns long (longitude), lat (latitude) and radius.
}

\author{
Aaron Robotham
}

\seealso{
\code{\link{sph2car}}
}

\examples{
print(car2sph(x=1,y=1,z=0,deg=TRUE))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{transform}

