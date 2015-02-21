\name{Cosmology Reference Sets}
\alias{cosref}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Cosmology parameter data included in celestial package.
}
\description{
cosref:   Cosmology H0 / OmegaM / OmegaL / Sigma8 parameters taken from Planck, WMAP1-9, Millennium Simulation and GiggleZ.
}
\usage{
data(cosref)
}
\details{
The included data is a table of the following cosmological parameters:
\tabular{lllll}{
Ref         \tab  H0    \tab  OmegaM  \tab  OmegaL  \tab  Sigma8  \cr
737         \tab  70.0  \tab  0.3000  \tab  0.7000  \tab  Test    \cr
Planck      \tab  67.0  \tab  0.3183  \tab  0.6817  \tab  0.8347  \cr
WMAP9       \tab  69.7  \tab  0.2819  \tab  0.7181  \tab  0.8170  \cr
WMAP7       \tab  70.4  \tab  0.2750  \tab  0.7250  \tab  0.8160  \cr
WMAP5       \tab  70.5  \tab  0.2740  \tab  0.7260  \tab  0.8120  \cr
WMAP3       \tab  70.4  \tab  0.2680  \tab  0.7320  \tab  0.7760  \cr
WMAP1       \tab  72.0  \tab  0.2900  \tab  0.7100  \tab  0.9000  \cr
Millennium  \tab  73.0  \tab  0.2500  \tab  0.7500  \tab  0.9000  \cr
GiggleZ     \tab  70.5  \tab  0.2730  \tab  0.7270  \tab  0.8120  \cr
}
}
\references{
\tabular{lll}{
Name    \tab  Full Reference                              \tab  arXiv Refence   \cr
737     \tab  Simplified concordance cosomology           \tab  NA              \cr
Planck  \tab  Planck Collaboration, 2014, A&A, 571, 16    \tab  arXiv:1303.5076 \cr
WMAP9	  \tab  Hinshaw G., et al., 2013, ApJS, 208, 19     \tab  arXiv:1212.5226 \cr
WMAP7	  \tab  Komatsu E., et al., 2010, ApJS, 192, 18     \tab  arXiv:1001.4538 \cr
WMAP5	  \tab  Komatsu E., et al., 2009, ApJS, 180, 306    \tab  arXiv:0803.0547 \cr
WMAP3	  \tab  Spergel D. N., et al., 2007, ApJS, 170, 377 \tab  arXiv:astro-ph/0603449 \cr
WMAP1	  \tab  Spergel D. N., et al., 2003, ApJS, 148, 175 \tab  arXiv:astro-ph/0302209 \cr
Millennium \tab	Springel V., et al., 2005, Nature, 435, 629 \tab arXiv:astro-ph/0504097\cr
GiggleZ	\tab  Poole G. B., et al., 201?, MNRAS, ?, ?      \tab  arXiv:1407.0390\cr
}
}
\author{
Aaron Robotham
}
\seealso{
  \code{\link{cosvol}}, \code{\link{cosmap}},  \code{\link{cosdist}}, \code{\link{cosgrow}}
}
\examples{
data(cosref)
cosref[cosref[,'Ref']=='Planck',]
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ fit }
\keyword{ data }

