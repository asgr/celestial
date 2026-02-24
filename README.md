# celestial (R package)

<!-- badges: start -->
![R-CMD-check](https://github.com/asgr/celestial/workflows/R-CMD-check/badge.svg)
<!-- badges: end -->

## Synopsis

Core package containing a collection of common astronomical conversion routines and functions. It provides tools for coordinate conversions (decimal degrees to/from HMS/DMS formats), cosmological calculations (distances, volumes, ages, structure growth), sky coordinate matching, WCS projections, halo profiles, and more.

## Installation

### Getting R

First things first, you will probably want to install a recent version of **R** that lets you build packages from source. The advantage of choosing this route is you can then update bleeding edge versions directly from GitHub. If you rely on the pre-built binaries on CRAN you might be waiting much longer.

#### Mac

For Mac just get the latest binaries from the **R** project pages:

<https://cloud.r-project.org/bin/macosx/>

#### Windows

For Windows just get the latest binaries from the **R** project pages:

<https://cloud.r-project.org/bin/windows/>

#### Linux

Debian:	`sudo apt-get install r-base r-base-dev`

Fedora:	`sudo yum install R`

Suse:	More of a pain, see here <https://cloud.r-project.org/bin/linux/suse/README.html>

Ubuntu:	`sudo apt-get install r-base-dev`

All the info on binaries is here: <https://cloud.r-project.org/bin/linux/>

Now you have **R** installed (hopefully) I would also suggest you get yourself **R-Studio**. It is a very popular and well maintained **R** IDE that gives you a lot of helpful shortcuts to scripting and analysing with **R**. The latest version can be grabbed from <https://www.rstudio.com/products/rstudio/> where you almost certainly want the free Desktop version.

### Getting celestial

Source installation from GitHub should be easy:

```R
install.packages('remotes')
remotes::install_github("asgr/celestial")
library(celestial)
```

A few Mac people seem to have issues with the above due to the backend used to download files. A work around seems to be to either use devtools:

```R
install.packages('devtools')
devtools::install_github("asgr/celestial")
library(celestial)
```

Or try the following:

```R
Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS="true")
remotes::install_github("asgr/celestial")
```

If all of these do not work then the nuclear option is to download (or clone) the GitHub repo, cd to where the tar.gz file is and run in the **console** (or **Terminal** on Mac):

```console
R CMD install celestial_X.Y.Z.tar.gz
```

where X, Y and Z should be set as appropriate for the version downloaded (check the name of the file basically).

Failing all of the above, please raise an Issue here.

#### Package Dependencies

The above should also install the required packages. If you have trouble with this you can try installing the required packages manually first and then retry the installation for **celestial**:

```R
install.packages(c('RANN', 'NISTunits', 'pracma')) # Required packages
install.packages(c('plotrix', 'igraph')) # Suggested packages
install.packages('remotes')
remotes::install_github("asgr/celestial")
```

Assuming you have installed all of the packages that you need/want, you should now be able to load **celestial** within **R** with the usual:

```R
library(celestial)
```

## Code Example

### Coordinate Conversions

Convert decimal degrees to HMS (Hours, Minutes, Seconds) format, as commonly used for right ascension:

```R
deg2hms(182.5, type='cat', sep=':')
# [1] "12:10:00.00"

deg2hms(182.5, type='cat', sep='hms')
# [1] "12h10m00.00s"
```

Convert decimal degrees to DMS (Degrees, Minutes, Seconds) format, as commonly used for declination:

```R
deg2dms(-27.345, type='cat', sep=':')
# [1] "-27:20:42.00"

deg2dms(-27.345, type='cat', sep='dms')
# [1] "-27d20m42.00s"
```

### Cosmological Distance Calculations

Compute cosmological distances and related quantities for a given redshift:

```R
cosdist(z=1, H0=70, OmegaM=0.3)
```

This returns a data.frame with columns including comoving distance, luminosity distance, angular diameter distance, look-back time, and more.

For quick lookups you can use the convenient single-output wrapper functions:

```R
cosdistLumDist(z=1, H0=70, OmegaM=0.3)   # Luminosity distance in Mpc
cosdistAngDist(z=1, H0=70, OmegaM=0.3)   # Angular diameter distance in Mpc
cosdistAngScale(z=1, H0=70, OmegaM=0.3)  # Angular scale in kpc/arcsec
```

You can also use built-in reference cosmologies:

```R
cosdist(z=1, ref='Planck18')
```

### Cosmological Volume Calculations

```R
cosvol(zlo=0, zhi=0.1, H0=70, OmegaM=0.3)
```

### Sky Coordinate Matching

Match two catalogues of sky coordinates:

```R
# Create two sets of RA/Dec coordinates
coordref = cbind(RA=c(10.1, 10.5, 11.0), Dec=c(-20.1, -20.5, -21.0))
coordcompare = cbind(RA=c(10.11, 10.6, 12.0), Dec=c(-20.09, -20.4, -22.0))

# Match within 10 arcsec
match = coordmatch(coordref, coordcompare, rad=10)
```

### WCS Projections

Convert between RA/Dec (degrees) and x/y pixel coordinates:

```R
radec2xy(RA=180.0, Dec=-30.0, CRVAL1=180.0, CRVAL2=-30.0, CRPIX1=512, CRPIX2=512,
         CD1_1=-2.78e-4, CD2_2=2.78e-4)
```

## Motivation

**celestial** was developed to provide a self-contained collection of common astronomical utility functions in **R**, particularly for cosmological calculations and sky coordinate handling. It is used as a dependency by several other **ProTools** packages including **ProFound** and **ProFit**.

## Contributors

Aaron Robotham

## References

Driver S.P. & Robotham A.S.G., 2010, MNRAS, 407, 2131 (cosmic variance calculator)

## Resources

<https://github.com/asgr/celestial>

## Forums

Please sign up to <http://profit.freeforums.net/> if you want to ask a question (or browse the questions asked).

## License

GPL-3
