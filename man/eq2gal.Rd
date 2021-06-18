\name{eq2gal}
\alias{eq2gal}
\alias{gal2eq}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Convert Equtorial to Galactic Coordinates
}
\description{
Simple conversions accurate to about 0.1 arcseconds between J2000 Equatorial and MW Galactic.
}
\usage{
eq2gal(RA, Dec, pole_RA = 192.859508, pole_Dec = 27.128336, eta = 32.932)
gal2eq(gal_long, gal_lat, pole_RA = 192.859508, pole_Dec = 27.128336, eta = 32.932)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{RA}{
Numeric vector; Right Ascension in degrees (J2000 system). For convenience this can also be a two column matrix / data.frame containing \option{RA} and \option{Dec} columns.
}
  \item{Dec}{
Numeric vector; Declination in degrees (J2000 system).
}
  \item{gal_long}{
Numeric vector; Galactic Longitude in degrees (J2000 system). For convenience this can also be a two column matrix / data.frame containing \option{gal_long} and \option{gal_lat} columns.
}
  \item{gal_lat}{
Numeric vector; Galactic Latitude in degrees (J2000 system).
}
  \item{pole_RA}{
Numeric scalar; Right Ascension of Galactic Pole in degrees (J2000 system).
}
  \item{pole_Dec}{
Numeric scalar; Declination of Galactic Pole in degrees (J2000 system).
}
  \item{eta}{
Numeric scalar; Tilt of the Earth's axis in degrees (J2000 system).
}
}
\details{
Approximate simple conversions, accurate usully to 0.1 arcseconds
}
\value{
\code{eq2gal} returns the Galactic longitude and latitude in degrees (two column output).
\code{gal2eq} returns the Equatorial RA and Dec in degress (two column output).
}
\references{
"Practical astronomy with your calculator" (Peter Duffett-Smith)
}
\author{
Aaron Robotham
}

\seealso{
\code{\link{sph2car}},\code{\link{car2sph}}
}
\examples{
# Notice the conversion anbiguity at the poles

gal2eq(eq2gal(seq(0,360, len=19), seq(-90,90, len=19)))
}
\concept{ coordinate }

