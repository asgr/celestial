\name{Cosmology Reference Sets}
\alias{cosref}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Cosmology parameter data included in celestial package.
}
\description{
cosref:   Cosmology H0 / OmegaM / OmegaL / OmegaR (via OmegaM/zeq) and Sigma8 parameters taken from Planck, WMAP1-9, Millennium Simulation and GiggleZ. Not all of these exist for each source, so NA values are used in these cases.
}
\usage{
data(cosref)
}
\details{
The included data is a table of the following cosmological parameters:
\tabular{llllll}{
Ref         \tab  H0    \tab  OmegaM \tab  OmegaL \tab  OmegaR  \tab  Sigma8  \cr
737         \tab  70.0  \tab  0.300  \tab  0.700  \tab  NA      \tab  NA      \cr
137         \tab  100.0 \tab  0.300  \tab  0.700  \tab  NA      \tab  NA      \cr
Planck13    \tab  67.3  \tab  0.315  \tab  0.685  \tab  9.289295e-05  0.829   \cr 
Planck15    \tab  67.8  \tab  0.308  \tab  0.692  \tab  9.150327e-05  0.815   \cr
Planck      \tab  67.8  \tab  0.308  \tab  0.692  \tab  9.150327e-05  0.815   \cr
WMAP9       \tab  69.7  \tab  0.288  \tab  0.712  \tab  8.780488e-05  \tab  0.817  \cr
WMAP7       \tab  70.4  \tab  0.275  \tab  0.725  \tab  8.569648e-05  \tab  0.816  \cr
WMAP5       \tab  70.5  \tab  0.274  \tab  0.726  \tab  8.45679e-05   \tab  0.812  \cr
WMAP3       \tab  70.4  \tab  0.268  \tab  0.732  \tab  NA      \tab  0.776   \cr
WMAP1       \tab  72.0  \tab  0.290  \tab  0.710  \tab  NA      \tab  0.900   \cr
Millennium  \tab  73.0  \tab  0.250  \tab  0.750  \tab  NA      \tab  0.900   \cr
GiggleZ     \tab  70.5  \tab  0.273  \tab  0.727  \tab  NA      \tab  0.812   \cr
}
}
\references{
\tabular{lll}{
Name    \tab  Full Reference                              \tab  arXiv Refence   \cr
737     \tab  Simplified concordance cosomology           \tab  NA              \cr
7378    \tab  Simplified concordance cosomology           \tab  NA              \cr
Planck  \tab  Planck Collaboration, 2014, A&A, 571, 16    \tab  arXiv:1303.5076v3 \cr
WMAP9	  \tab  Hinshaw G., et al., 2013, ApJS, 208, 19     \tab  arXiv:1212.5226v3 \cr
WMAP7	  \tab  Komatsu E., et al., 2010, ApJS, 192, 18     \tab  arXiv:1001.4538v3 \cr
WMAP5	  \tab  Komatsu E., et al., 2009, ApJS, 180, 306    \tab  arXiv:0803.0547v2 \cr
WMAP3	  \tab  Spergel D. N., et al., 2007, ApJS, 170, 377 \tab  arXiv:astro-ph/0603449v2 \cr
WMAP1	  \tab  Spergel D. N., et al., 2003, ApJS, 148, 175 \tab  arXiv:astro-ph/0302209v3 \cr
Millennium \tab	Springel V., et al., 2005, Nature, 435, 629 \tab arXiv:astro-ph/0504097v2 \cr
GiggleZ	\tab  Poole G. B., et al., 201?, MNRAS, ?, ?      \tab  arXiv:1407.0390v1 \cr
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

