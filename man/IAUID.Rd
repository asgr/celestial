\name{IAUID}
\alias{IAUID}

\title{
IAU name creator.
}

\description{
Creates IAU legal names for objects given coordinates, name and epoch.
}

\usage{
IAUID(ra, dec, name = "GAMA", epoch = "J")
}

\arguments{
  \item{ra}{
The decimal degrees right ascension.
}
  \item{dec}{
The decimal degrees declination.
}
  \item{name}{
Name to be appended to IAU designation as a string.
}
  \item{epoch}{
Epoch, i.e. 'J' (default) or 'B'. Enter as a string.
}
}

\value{
Text string that outputs an IAU legal name for an object.
}

\author{
Aaron Robotham
}

\examples{
print(IAUID(123.45,67.89,'GAMA','J'))
}

\keyword{convert}

