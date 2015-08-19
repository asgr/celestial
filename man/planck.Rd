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
\alias{cosplanckSBLawRad_sr}
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
cosplanckSBLawRad_sr(Temp=2.725)
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
The functions with \code{Rad} in the name are related the spectral radiance form of Planck's Law (typically designated I or B), whilst those with \code{En} are related to the spectral energy density form of Planck's Law (u), where \eqn{u=4\pi I/c}{u=4.pi.I/c}.

To calculate the number of photons in a mode we simply use \eqn{E=h\nu=h c / \lambda}{E=h.nu=h.c/lambda}.

Below h is the Planck constant, \eqn{k_B}{k[B]} is the Boltzmann constant, c is the speed-of-light in a vacuum and \eqn{\sigma}{sigma} is the Stefan-Boltzmann constant.

\code{cosplanckLawRadFreq} is the spectral radiance per unit frequency version of Planck's Law, defined as:

\deqn{B_\nu(\nu,T) = I_\nu(\nu,T) = \frac{2 h \nu^3}{c^2} \frac{1}{e^{h \nu / k_B T}-1}}{%
      B[nu](nu,T) = I[nu](nu,T) = (2.h.nu^3/c^2).(1/(exp(h.nu/k[B].T)-1))}
      
\code{cosplanckLawRadWave} is the spectral radiance per unit wavelength version of Planck's Law, defined as:

\deqn{B_\lambda(\lambda,T) = I_\lambda(\lambda,T) = \frac{2 h c^2}{\lambda^5} \frac{1}{e^{h c / \lambda k_B T}-1}}{%
      B[lambda](lambda,T) = I[lambda](lambda,T) = (2.h.c^2/lambda^5).(1/(exp(h.c/lambda.k[B].T)-1))}

\code{cosplanckLawRadFreqN} is the number of photons per unit frequency, defined as:

\deqn{B_\nu(\nu,T) = I_\nu(\nu,T) = \frac{2 \nu^2}{c^2} \frac{1}{e^{h \nu / k_B T}-1}}{%
      B[nu](nu,T) = I[nu](nu,T) = (2.nu^2/c^2).(1/(exp(h.nu/k[B].T)-1))}
      
\code{cosplanckLawRadWaveN} is the number of photons per unit wavelength, defined as:

\deqn{B_\lambda(\lambda,T) = I_\lambda(\lambda,T) = \frac{2 c}{\lambda^4} \frac{1}{e^{h c / \lambda k_B T}-1}}{%
      B[lambda](lambda,T) = I[lambda](lambda,T) = (2.c/lambda^4).(1/(exp(h.c/lambda.k[B].T)-1))}
      
\code{cosplanckLawEnFreq} is the spectral energy density per unit frequency version of Planck's Law, defined as:

\deqn{u_\nu(\nu,T) = \frac{8 \pi h \nu^3}{c^3} \frac{1}{e^{h\nu/k_B T}-1}}{%
      u[nu](nu,T) = (8.pi.h.nu^3/c^3).(1/(exp(h.nu/k[B].T)-1))}
      
\code{cosplanckLawEnWave} is the spectral energy density per unit wavelength version of Planck's Law, defined as:

\deqn{u_\lambda(\lambda,T) = \frac{8 \pi h c}{\lambda^5} \frac{1}{e^{h c / \lambda k_B T}-1}}{%
      u[lambda](lambda,T) = (8.pi.h.c/lambda^5).(1/(exp(h.c/lambda.k[B].T)-1))}

\code{cosplanckPeakFreq} gives the location in frequency of the peak of \eqn{I_\nu(\nu,T)}{I[nu](nu,T)}, defined as:

\deqn{\nu_{peak} = 2.821 k_B T}{%
      nu[peak] = 2.821.k[B].T}
      
\code{cosplanckPeakWave} gives the location in wavelength of the peak of \eqn{I_\lambda(\lambda,T)}{I[lambda](lambda,T)}, defined as:

\deqn{\lambda_{peak} = 4.965 k_B T}{%
      lambda[peak] = 4.965.k[B].T}

\code{cosplanckSBLawRad} gives the emissive power (or radiant exitance) version of the Stefan-Boltzmann Law, defined as:

\deqn{j^* = \sigma T^4}{%
      j^* = sigma.T^4}
      
\code{cosplanckSBLawRad_sr} gives the spectral radiance version of the Stefan-Boltzmann Law, defined as:

\deqn{L = \sigma T^4/\pi}{%
      L = sigma.T^4/pi}
      
\code{cosplanckSBLawEn} gives the energy density version of the Stefan-Boltzmann Law, defined as:

\deqn{\epsilon = 4 \sigma T^4 / c}{%
      epsilon = 4.sigma.T^4/c}
      
Notice that \eqn{j^*}{J^*} and L merely differ by a factor of \eqn{\pi}{pi}, i.e. L is per steradian.

\code{cosplanckLawRadPhotEnAv} gives the average energy of the emitted black body photon, defined as:

\deqn{<E_{phot}> = 3.729282 \times 10^{-23} T}{%
      <E[phot]> = 3.729282e-23 T}

\code{cosplanckLawRadPhotN} gives the total number of photons produced by black body per metre squared per second per steradian, defined as:

\deqn{N_{phot} = 1.5205 \times 10^{15} T^3 / \pi}{%
      N[phot] = 1.5205e+15.T^3/pi}

Various confidence building sanity checks of how to use these functions are given in the Examples below.
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
The number of photons per steradian per metre squared per second per unit frequency for a black body (photons.sr\eqn{^{-1}}{^-1}.m\eqn{^{-2}}{^-2}.s\eqn{^{-1}}{^-1}.Hz\eqn{^{-1}}{^-1}).
}

\item{cosplanckLawRadWaveN}{
The number of photonsper steradian per metre squared per second per unit wavelength for a black body (photons.sr\eqn{^{-1}}{^-1}.m\eqn{^{-2}}{^-2}.s\eqn{^{-1}}{^-1}.m\eqn{^{-1}}{^-1}).
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
Total energy radiated per metre squared per second across all wavelengths for a black body (W.m\eqn{^{-2}}{^-2}.sr\eqn{^{-1}}{^-1}). This is the power version of the Stefan-Boltzmann Law.
}

\item{cosplanckSBLawRad_sr}{
Total energy radiated per metre squared per second per steradian across all wavelengths for a black body (W.m\eqn{^{-2}}{^-2}). This is the radiance version of the Stefan-Boltzmann Law.
}

\item{cosplanckSBLawEn}{
Total energy per metre cubed across all wavelengths for a black body (J.m\eqn{^{-3}}{^-3}). This is the energy density version of the Stefan-Boltzmann Law.
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
#Classic example for different temperature stars:

waveseq=10^seq(-7,-5,by=0.01)
plot(waveseq, cosplanckLawRadWave(waveseq,5000),
log='x', type='l', xlab=expression(Wavelength / m),
ylab=expression('Spectral Radiance' / W*sr^{-1}*m^{-2}*m^{-1}), col='blue')
lines(waveseq, cosplanckLawRadWave(waveseq,4000), col='green')
lines(waveseq, cosplanckLawRadWave(waveseq,3000), col='red')
legend('topright', legend=c('3000K','4000K','5000K'), col=c('red','green','blue'), lty=1)

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

cosplanckSBLawRad_sr()

#We can get (almost, some rounding is off) the same answer by multiplying
#the total number of photons produced by a black body per metre squared per
#second per steradian by the average photon energy:

cosplanckLawRadPhotEnAv()*cosplanckLawRadPhotN()

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{planck}
\keyword{blackbody}% __ONLY ONE__ keyword per line
