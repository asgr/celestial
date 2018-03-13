\name{skyproj}
\alias{tanproj}
\alias{sinproj}
\alias{radec2xy}
\alias{xy2radec}
\alias{gnomonic}
\alias{orthographic}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Tan Gnomonic and Sine Orthographic Projection System WCS Solver Functions
}
\description{
Converts RA/Dec (degrees) to x/y (pixels) position using the Tan Gnomonic or Sine Orthographic projection systems, and vice-versa.
Translations adapted from: http://mathworld.wolfram.com/GnomonicProjection.html and http://mathworld.wolfram.com/OrthographicProjection.html.
}
\usage{
radec2xy(RA, Dec, header, CRVAL1 = 0, CRVAL2 = 0, CRPIX1 = 0, CRPIX2 = 0, CD1_1 = 1,
CD1_2 = 0, CD2_1 = 0, CD2_2 = 1, proj='tan')
xy2radec(x, y, header, CRVAL1 = 0, CRVAL2 = 0, CRPIX1 = 0, CRPIX2 = 0, CD1_1 = 1,
CD1_2 = 0, CD2_1 = 0, CD2_2 = 1, proj='tan')
}

\arguments{
  \item{RA}{
Vector or matrix; target right ascension in degrees. If matrix then the first column will be used as RA and the second column as Dec.
}
  \item{Dec}{
Vector; target declination in degrees. Ignored if \option{RA} is a matrix.
}
  \item{x}{
Vector or matrix; target x-pixel. If Matrix then the first column will be used as the x-axis and the second column as y-axis.
}
  \item{y}{
Vector; target y-pixel. Ignored if \option{x} is a matrix.
}
  \item{CRVAL1}{
FITS header CRVAL1 for the Tan Gnomonic projection system. This is the RA in degrees at the location of \option{CRPIX1}.
}
  \item{CRVAL2}{
FITS header CRVAL2 for the Tan Gnomonic projection system. This is the Dec in degrees at the location of \option{CRPIX2}.
}
  \item{CRPIX1}{
FITS header CRPIX1 for the Tan Gnomonic projection system. This is the x pixel value at the location of \option{CRVAL1}.
}
  \item{CRPIX2}{
FITS header CRPIX2 for the Tan Gnomonic projection system. This is the y pixel value at the location of \option{CRVAL2}.
}
  \item{CD1_1}{
FITS header CD1_1 for the Tan Gnomonic projection system. Change in RA-Tan in degrees along x-Axis.
}
  \item{CD1_2}{
FITS header CD1_2 for the Tan Gnomonic projection system. Change in RA-Tan in degrees along y-Axis.
}
  \item{CD2_1}{
FITS header CD2_1 for the Tan Gnomonic projection system. Change in Dec-Tan in degrees along x-Axis.
}
  \item{CD2_2}{
FITS header CD2_2 for the Tan Gnomonic projection system. Change in Dec-Tan in degrees along y-Axis.
}
  \item{CTYPE1}{
The RA projection system type. Either 'RA--TAN' for Tan Gnomonic, or 'RA--SIN' for Sine Orthographic. 'RA--NCP' is approximated by Sine Orthographic with a warning. Over-ridden by the FITS header.  
}
  \item{CTYPE2}{
The DEC projection system type. Either 'DEC--TAN' for Tan Gnomonic, or 'DEC--SIN' for Sine Orthographic. 'DEC--NCP' is approximated by Sine Orthographic with a warning. Over-ridden by the FITS header.  
}
  \item{header}{
Full FITS header in table or vector format. Legal table format headers are provided by the \code{read.fitshdr} function or the \option{hdr} list output of \code{read.fits} in the \code{astro} package). Also the \option{hdr} output of \code{readFITS} in the \code{FITSio} package provides legal vector format inputs. If a header is provided then key words will be taken from here as a priority. Missing header keywords are printed out and other header option arguments are used in these cases.
}
}
\details{
These functions encode the standard FITS Tan Gnomonic and Sine Orthographic projection systems for solving an image WCS (covering most moden imaging and radio data). They do not deal with higher order polynomial distortion terms.
}
\value{
  \item{radec2xy}{
Returns a two column matrix with columns x and y.
}
  \item{xy2radec}{
Returns a two column matrix with columns RA and Dec (in degrees).
}
}
\references{
http://mathworld.wolfram.com/GnomonicProjection.html
http://mathworld.wolfram.com/OrthographicProjection.html
}

\author{
Aaron Robotham
}

\seealso{
\code{\link{deg2dms}}, \code{\link{deg2hms}}, \code{\link{dms2deg}}, \code{\link{hms2deg}}
}
\examples{
#A simple example:

radec2xy(10,20)
xy2radec(radec2xy(10,20))

#A more complicated example, where we transform and rotate large amounts:

exdata_start=expand.grid(1:10,21:30)
plot(exdata_start)
exradec=radec2xy(exdata_start, CRVAL1=20, CRPIX1=100, CRVAL2=30, CRPIX2=130, CD1_1=0.1,
CD1_2=-0.05, CD2_1=0.05, CD2_2=0.1)
plot(exradec)
exdata_end=xy2radec(exradec, CRVAL1=20, CRPIX1=100, CRVAL2=30, CRPIX2=130, CD1_1=0.1,
CD1_2=-0.05, CD2_1=0.05, CD2_2=0.1)
plot(exdata_start,cex=2)
points(exdata_end,col='red')

#The residuals should be very small (in the noice of double precision arithmetic):

plot(density(exdata_start[,1]-exdata_end[,1]))
lines(density(exdata_start[,2]-exdata_end[,2]),col='red')
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ tan }% use one of  RShowDoc("KEYWORDS")
\keyword{ projection }% __ONLY ONE__ keyword per line
\keyword{ gnomonic }
