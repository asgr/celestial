\name{cosorb}
\alias{cosorbVisViva}
\alias{cosorbFreeFall}
\alias{cosorbRocheRad}
\alias{cosorbRocheSize}
\alias{potential_part}
\alias{kinetic_part}

\title{
  Orbital functions
}
\description{
  A variety of obital analysis functions. These are useful for setting up initial conditions for merging systems etc.
}
\usage{
cosorbVisViva(mass=1e12, Rad=162.635, SemiMajRad=162.635, Munit=1, Lunit=1e3, Vunit=1)
cosorbFreeFall(M1=1e12, M2=1, Rad=162.635, Munit=1, Lunit=1e3, Vunit=1, Tunit=1e9)
cosorbRocheRad(M1=1e12, M2=1e10, Size=35.03865, Rfac=2.44)
cosorbRocheSize(M1=1e12, M2=1e10, Rad=396.8294, Rfac=2.44)
  
potential_part(part, eval = NULL, mass = 1, soft = 1, Munit = 1, Lunit = 1000,
  Vunit = 1) 
kinetic_part(vel, mass=NULL, Vunit = 1)
  
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mass}{
    Mass in units of 'Munit'.
  }
  \item{M1}{
    Mass of primary body in units of 'Munit'.
  }
  \item{M2}{
    Mass of secondary body in units of 'Munit'.
  }
  \item{Rad}{
    Separation between bodies in units of 'Lunit' (for cosorbRocheSize this is in arbitrary units).
  }
  \item{part}{
    Cartesian particle positions [x,y,z] in units of 'Lunit'.
  }
  \item{eval}{
    Cartesian evalution positions [x,y,z] in units of 'Lunit'. If this is not provided then the potentials are evaluated at the positions of \option{part}.
  }
  \item{vel}{
    Cartesian particle velocities [vx,vy,vz] in units of 'Vunit'.
  }
  \item{soft}{
    Numeric softening in units of 'Lunit'.
  }
  \item{SemiMajRad}{
    The semi major radius of the orbit (a > 0 for ellipses, a = Rad for circles, 1/a = 0 for parabolas, and a < 0 for hyperbolas).
  }
  \item{Size}{
    The size radius of the secondary object. Inside of this radius the object is bound to the secondary, outside of this radius the object is stripped by the primary.
  }
  \item{Rfac}{
    The Roche factor. Approximately taken to be 2.44, but in reality it varies depending on the shape of the potentials etc.
  }
  \item{Munit}{
    Base mass unit in multiples of Msun.
  }
  \item{Lunit}{
    Base length unit in multiples of parsecs.
  }
  \item{Vunit}{
    Base velocity unit in multiples of km/s.
  }
  \item{Tunit}{
    Base time unit in multiples of years.
  }
}

\details{
  These functions allow for various analytic conversions between the 3 major properties related to virial radius: the mass, velocity dispresion and size. The default properties calculate properties for 1e12 Msun halos and assume masses in Msun, velocities in km/s and distances in Kpc.
}

\value{
  cosorbVisViva function gives the required velocity in units of Vunit to create the specified orbit.
  
  cosorbFreeFall function gives the free fall time to static initial velocity separated bodies.
  
  cosorbRocheRad function gives the orbital radius at which the secondary will become stripped within a specified bound radius.
  
  cosorbRocheSize function gives the limiting bound radius of the secondary for a given system.
  
  potential_part functions gives the total gravitational potential either at the positions of particles provided (\option{part}) or at the requested evaluation positions (\option{eval}).
  
  kinetic_part functions gives the total kinetic energy of particles provided (\option{vel})
}

\author{
  Aaron Robotham, Chris Power
}

\seealso{
  \code{\link{cosvol}}, \code{\link{cosmap}},  \code{\link{cosdist}}, \code{\link{cosgrow}}
}

\examples{
  cosorbVisViva(M=1e15, Rad=1, Lunit=1e6)
  cosorbFreeFall(M1=1e15, M2=1, Rad=1, Lunit=1e6)
  cosorbRocheRad(M1=1e12, M2=1e12, Size=162.635, Rfac=2.44)
  cosorbRocheSize(M1=1e12, M2=1e12, Rad=396.8294, Rfac=2.44)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{halo}
\keyword{cosmo}
\keyword{cosmology}


