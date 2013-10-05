\name{dms2deg}
\alias{dms2deg}

\title{
Convert DMS to degrees format.
}
\description{
Convert DMS (degrees, minutes, seconds) to degrees format. This is probably most useful for declination conversion, since dms is fairly standard method of presenting declination coordinates. The decimal degrees=d+m/60+s/3600. Degrees should range from -90 to +90. Degrees and minutes should be integer and seconds can be decimal.
}
\usage{
dms2deg(d,m,s,sign='d',sep=':')
}

\arguments{
  \item{d}{
The integer number of degrees you are converting. If it is not integer then the floor of the number is taken. This can contain the sign of the declination when sign='d', but must be all positive if the sign argument is specified (this is required if d contains any 0s, see below). If sign is specified, all d values should be 0<=d<=90, otherwise d values should be 0<=d<=90.
}
  \item{m}{
The integer number of minutes you are converting. If it is not integer then the floor of the number is taken. All m values should be 0<=m<60.
}
  \item{s}{
The decimal number of seconds you are converting. All s values should be 0<=s<60.
}
  \item{sign}{
The sign of the declination. The default 'd' inherits the sign of the d argument. This is ambiguous when d is 0 since the sign of +/-0 is taken to be 0. If d contains any 0s, you must supply a vector of the same length as d with +ve or -ve values (e.g. +/- 1), the sign of these value will be taken as the sign for the declination.
}
  \item{sep}{
Defines the type of separator used when 'd' is a vector of strings. Any value other than 'DMS' and 'dms' is used for all separations, so the default ':' would be for an input like 3:34:45.5. If set to 'dms' or 'DMS' then the output is of the format 3d34m45.5s and 3D34M45.5s resepctively.
}
}

\value{
A value of decimal degrees.
}

\author{
Aaron Robotham
}

\seealso{
\code{\link{deg2dms}}
}

\examples{
print(dms2deg(70,45,19,-1))
print(dms2deg('-70:45:19'))
print(dms2deg('-70d45m19s',sep='dms'))
print(dms2deg(c('-70D45M19S','3D5M15S'),sep='DMS'))
}

\keyword{convert}

