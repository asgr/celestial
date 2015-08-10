\name{planck}
\alias{planck}
\alias{cosplanck}
\alias{cosplanckLawRadFreq}
\alias{cosplanckLawRadWave}
\alias{cosplanckLawEnFreq}
\alias{cosplanckLawEnWave}
\alias{cosplanckLawRadFreqN}
\alias{cosplanckLawRadWaveN}
\alias{cosplanckPeakFreq}
\alias{cosplanckPeakWave}
\alias{cosplanckSBLawRad}
\alias{cosplanckSBLawEn}
\alias{cosplanckLawRadPhotEnAv}
\alias{cosplanckLawRadPhotN}
\alias{cosplanckCMBTemp}
\title{
Planck's Law and Related Functions
}
\description{
Functions related to Planck's Law of thermal radiation.
}
\usage{
cosplanckLawRadFreq(nu,Temp=2.725)
cosplanckLawRadWave(lambda,Temp=2.725)
cosplanckLawEnFreq(nu,Temp=2.725)
cosplanckLawEnWave(lambda,Temp=2.725)
cosplanckLawRadFreqN(nu,Temp=2.725)
cosplanckLawRadWaveN(lambda,Temp=2.725)
cosplanckPeakFreq(Temp=2.725)
cosplanckPeakWave(Temp=2.725)
cosplanckSBLawRad(Temp=2.725)
cosplanckSBLawEn(Temp=2.725)
cosplanckLawRadPhotEnAv(Temp=2.725)
cosplanckLawRadPhotN(Temp=2.725)
cosplanckCMBTemp(z)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{nu}{
The frequency of radiation in Hertz (Hz).
}
  \item{lambda}{
The wavelength of radiation in metres (m).
}
  \item{Temp}{
The absolute temperature of the system in Kelvin (K).
}
  \item{z}{
Redshift, where z must be > -1 (can be a vector).
}
}
\details{
The functions with \code{Rad} in the name are related the spectral radiance form of Planck's Law (typically designated I or B), whilst those with \code{En} are related to the spectral energy density form of Planck's Law (u), where u=I.4.pi/c
}
\value{
Planck's Law in terms of spectral radiance:

\item{cosplanckLawRadFreq}{
The power per steradian per metre squared per unit frequency for a black body (W.sr\eqn{^{-1}}{^-1}.m\eqn{^{-2}}{^-2}.Hz\eqn{^{-1}}{^-1}).
}

\item{cosplanckLawRadWave}{
The power per steradian per metre squared per unit wavelength for a black body (W.sr\eqn{^{-1}}{^-1}.m\eqn{^{-2}}{^-2}.m\eqn{^{-1}}{^-1}).
}

Planck's Law in terms of spectral energy density:

\item{cosplanckLawEnFreq}{
The energy per metre cubed per unit frequency for a black body (J.m\eqn{^{-3}}{^-3}.Hz\eqn{^{-1}}{^-1}).
}

\item{cosplanckLawEnWave}{
The energy per metre cubed per unit wavelength for a black body (J.m\eqn{^{-3}}{^-3}.m\eqn{^{-1}}{^-1}).
}

Photon counts:

\item{cosplanckLawRadFreqN}{
The number of photons per unit frequency for a black body (Hz\eqn{^{-1}}{^-1}).
}

\item{cosplanckLawRadWaveN}{
The number of photons per unit wavelength for a black body (m\eqn{^{-1}}{^-1}).
}

Peak locations (via Wien's displacement law):

\item{cosplanckPeakFreq}{
The frequency location of the radiation peak for a black body as found in \code{cosplanckLawRadFreq}.
}

\item{cosplanckPeakWave}{
The wavelength location of the radiation peak for a black body as found in \code{cosplanckLawRadWave}.
}

Stefan-Boltzmann Law:

\item{cosplanckSBLawRad}{
Total energy radiated per unit surface area per unit time per steradian across all wavelengths for a black body (W.m\eqn{^{-2}}{^-2}). This is the radiance version of the relation.
}

\item{cosplanckSBLawEn}{
Total energy radiated per unit surface area per unit time across all wavelengths for a black body (W.m\eqn{^{-2}}{^-2}.sr\eqn{^{-1}}{^-1}). This is the power version of the relation.
}

Photon properties:

\item{cosplanckLawRadPhotEnAv}{
Average black body photon energy (J).
}

\item{cosplanckLawRadPhotN}{
Total number of photons produced by black body per metre squared per second per steradian (m\eqn{^{-2}}{^-2}.s\eqn{^{-1}}{^-1}.sr\eqn{^{-1}}{^-1}).
}

Cosmic Microwave Background:

\item{cosplanckCMBTemp}{
The temperaure of the CMB at redshift z.
}
}
\references{
Marr J.M., Wilkin F.P., 2012, AmJPh, 80, 399
}
\author{
Aaron Robotham
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{cosgrow}}
}
\examples{
#CMB now:

plot(10^seq(9,12,by=0.01), cosplanckLawRadFreq(10^seq(9,12,by=0.01)),
log='x', type='l', xlab=expression(Frequency / Hz),
ylab=expression('Spectral Radiance' / W*sr^{-1}*m^{-2}*Hz^{-1}))
abline(v=cosplanckPeakFreq(),lty=2)

plot(10^seq(-4,-1,by=0.01), cosplanckLawRadWave(10^seq(-4,-1,by=0.01)),
log='x', type='l', xlab=expression(Wavelength / m),
ylab=expression('Spectral Radiance' / W*sr^{-1}*m^{-2}*m^{-1}))
abline(v=cosplanckPeakWave(),lty=2)

#CMB at surface of last scattering:

TempLastScat=cosplanckCMBTemp(1100) #Note this is still much cooler than our Sun!

plot(10^seq(12,15,by=0.01), cosplanckLawRadFreq(10^seq(12,15,by=0.01),TempLastScat),
log='x', type='l', xlab=expression(Frequency / Hz),
ylab=expression('Spectral Radiance' / W*sr^{-1}*m^{-2}*Hz^{-1}))
abline(v=cosplanckPeakFreq(TempLastScat),lty=2)

plot(10^seq(-7,-4,by=0.01), cosplanckLawRadWave(10^seq(-7,-4,by=0.01),TempLastScat),
log='x', type='l', xlab=expression(Wavelength / m),
ylab=expression('Spectral Radiance' / W*sr^{-1}*m^{-2}*m^{-1}))
abline(v=cosplanckPeakWave(TempLastScat),lty=2)

#Exact number of photons produced by black body:

cosplanckLawRadPhotN()

#We can get pretty much the correct answer through direct integration, i.e.:

integrate(cosplanckLawRadFreqN,1e8,1e12)
integrate(cosplanckLawRadWaveN,1e-4,1e-1)

#Stefan-Boltzmann Law:

cosplanckSBLawRad()

#We can get (almost, some rounding is off) the same answer by multiplying
#the total number of photons produced by a black body per metre squared per
#second per steradian by the average photon energy:

cosplanckLawRadPhotEnAv()*cosplanckLawRadPhotN()

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{planck}
\keyword{blackbody}% __ONLY ONE__ keyword per line
