\name{hms2deg}
\alias{hms2deg}

\title{
Convert hms to degrees format.
}
\description{
Convert hms (hours, minutes, seconds) to degrees format. This is probably most useful for right ascension (RA) conversion, since hms is fairly standard method of presenting RA coordinates. The decimal degrees=15*h+15*m/60+15*s/3600. Should range between 0 and 24 hours. Hours and minutes should be integer and seconds can be decimal.
}
\usage{
hms2deg(h,m,s,sep=':')
}

\arguments{
  \item{h}{
The integer number of hours you are converting. If it is not integer then the floor of the number is taken. All m values should be 0<=h<=24.
}
  \item{m}{
The integer number of minutes you are converting. If it is not integer then the floor of the number is taken. All m values should be 0<=m<60.
}
  \item{s}{
The decimal number of seconds you are converting. All s values should be 0<=s<60.
}
  \item{sep}{
Defines the type of separator used when 'h' is a vector of strings. Any value other than 'HMS' and 'hms' is used for all separations, so the default ':' would be for an input like 3:34:45.5. If set to 'hms' or 'HMS' then the output is of the format 3h34m45.5s and 3H34M45.5s resepctively.
}
}

\value{
A value of decimal degrees.
}

\author{
Aaron Robotham
}

\seealso{
\code{\link{deg2hms}}
}

\examples{
print(hms2deg(12,10,36))
print(hms2deg('12:10:36'))
print(hms2deg('12h10m36s',sep='hms'))
print(hms2deg(c('12H10M36S','3H4M10S'),sep='HMS'))
}

\keyword{convert}

