\name{deg2dms}
\alias{deg2dms}

\title{
Convert decimal degrees to dms format.
}

\description{
Convert decimal degrees to dms (degrees, minutes, seconds) format. This is probably most useful for declination conversion, since dms is fairly standard method of presenting declination coordinates. The decimal degrees=d+m/60+s/3600. Degrees should range from -90 to +90.
}

\usage{
deg2dms(deg, type='mat', sep=':', digits=2)
}

\arguments{
  \item{deg}{
The decimal degrees you are converting. All deg values should be -90<=deg<=90
}
  \item{type}{
The output type desired. If 'mat' then the output is a 3 column data.frame where column 1 is the degree, column 2 is the minutes and column 3 is the seconds. If 'cat' then the output is a single vector of strings where the separator is defined by the 'sep' argument.
}
  \item{sep}{
Defines the type of separator used when type='cat'. Any value other than 'DMS' and 'dms' is used for all separations, so the default ':' would produce an output like 3:34:45.5. If set to 'dms' or 'DMS' then the output is of the format 3d34m45.5s and 3D34M45.5s resepctively.
}
  \item{digits}{
The digits to print for angular seconds. See \code{\link[base]{formatC}} for details on how digits is parsed.
}
}

\value{
A data.frame with the columns degrees, minutes and seconds if type='mat'. If type='cat' then a vector of strings with separators defined by the 'sep' argument.
}

\author{
Aaron Robotham
}

\seealso{
\code{\link{dms2deg}}
}

\examples{
print(deg2dms(12.345))
print(deg2dms(12.345,type='cat',sep=':'))
print(deg2dms(12.345,type='cat',sep='dms'))
print(deg2dms(12.345,type='cat',sep='DMS'))
}

\keyword{convert}

