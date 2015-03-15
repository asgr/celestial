\name{deg2hms}
\alias{deg2hms}

\title{
Convert decimal degrees to hms format.
}

\description{
Convert decimal degrees to hms (hours, minutes, seconds) format. This is probably most useful for right-ascension (RA) conversion, since hms is fairly standard method of presenting RA coordinates. The decimal degrees=15*h+15*m/60+15*s/3600 (i.e. there are 24 hours in 360 degrees). Degrees should range from 0 to 360.
}

\usage{
deg2hms(deg,type='mat',sep=':')
}

\arguments{
  \item{deg}{
The decimal degrees you are converting. All deg values should be 0<=d<=360.
}
  \item{type}{
The output type desired. If 'mat' then the output is a 3 column data.frame where column 1 is the degree, column 2 is the minutes and column 3 is the seconds. If 'cat' then the output is a single vector of strings where the separator is defined by the 'sep' argument.
}
  \item{sep}{
Defines the type of separator used when type='cat'. Any value other than 'DMS' and 'dms' is used for all separations, so the default ':' would produce an output like 3:34:45.5. If set to 'hms' or 'HMS' then the output is of the format 3h34m45.5s and 3H34M45.5s resepctively.
}
}

\value{
A data.frame with the columns degrees, minutes and seconds if type='mat'. If type='cat' then a vector of strings with separators defined by the 'sep' argument.
}

\author{
Aaron Robotham
}

\seealso{
\code{\link{hms2deg}}
}

\examples{
deg2hms(123.456)
deg2hms(123.456,type='cat',sep=':')
deg2hms(123.456,type='cat',sep='hms')
deg2hms(123.456,type='cat',sep='HMS')
}

\keyword{convert}

