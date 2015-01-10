\name{skyarea}
\alias{skyarea}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Exact angular area calculator
}
\description{
This function takes a survey geometry defined by RA (long) and Dec (latitude) limits and calculates the exact angular area covered.
}
\usage{
skyarea(long = c(129, 141), lat = c(-2, 3), inunit = "deg", outunit = "deg2")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{long}{
Upper and lower longitude (RA) limits of interest in units of inunit. If of length 1 then the number specified is assumed to be the upper limit and the lower limit is set to 0.
}
  \item{lat}{
Upper and lower latitude (Dec) limits of interest in units of inunit. If of length 1 then the number specified is assumed to be the upper limit and the lower limit is set to 0.
}
  \item{inunit}{
The units of angular coordinate provided. Allowed options are deg for degress, amin for arc minutes, asec for arc seconds and rad for radians.
}
  \item{outunit}{
The units of angular area desired. Allowed options are deg2 for square degrees, amin2 for square arc minutes, asec2 for square arc seconds and rad2 or sr for steradians.
}
}

\value{
Two value vector. First value is the sky area covered in units of outunit (names area), second value is the fraction of the celestial sphere covered by the specified geometry (named areafrac).
}
\references{
Based on the equations in Hogg D.W., 1999, arXiv, 9905116 and Wright E.L., 2006, PASP, 118, 1711
}
\author{
Aaron Robotham
}

\seealso{
\code{\link{cosvol}}
}
\examples{
#The GAMA survey areas:
G02area=skyarea(c(30.2,38.8),c(-10.25,-3.72))
G09area=skyarea(c(129,141),c(-2,3))
G12area=skyarea(c(174,186),c(-3,2))
G15area=skyarea(c(211.5,223.5),c(-2,3))
G23area=skyarea(c(338.1,351.9),c(-35,-30))

#Total GAMA survey area:
TotalGAMAarea=G02area+G09area+G12area+G15area+G23area
print(paste('The GAMA survey area is',round(TotalGAMAarea['area'],2),'sq. deg.'))

#Future TACs note: this is less than 1% of the sky ;-)
print(paste('The GAMA survey area is',round(TotalGAMAarea['areafrac']*100,2),'\% of the sky'))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{sky}
\keyword{area}% __ONLY ONE__ keyword per line
