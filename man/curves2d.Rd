\name{curves2d}
\alias{curves2d}
\title{Select points to "slide" along two-dimensional curves.}
\usage{
  curves2d(file, nsliders)
}
\arguments{
  \item{file}{Name of sample file in working directory (can
  include path) containing two-dimensional landmark data.
  File names can be written in manually, including paths,
  or obtained using directory/file manipulation functions
  e.g., \code{\link{list.files}}}

  \item{nsliders}{Number of landmarks to slide along
  curves.}
}
\value{
  Function returns an n-x-3 .nts file containing the
  positions along which each of n chosen semilandmarks will
  "slide". e.g., "4 3 2", semilandmark 3 will slide between
  4 and 2.
}
\description{
  An interactive function to select 2d "sliding"
  semilandmarks along curves.
}
\details{
  This function selects points on the created template to
  "slide" along 2d curves lacking known landmarks (see
  Bookstein 1991:376-382, 1997 for algorithm details). Each
  sliding semi-landmark (sliders) will slide between two
  designated points, along a line tangent to the specified
  curvature. Using the right mouse button, users: 1. Select
  sample ".nts" file digitized using function digitize2d 2.
  Select the first point between which semi-landmark will
  "slide" 3. Select sliding semi-landmark, 4.  Select point
  along which sliding trajectory will end. Screen plot
  shows arrows connecting the three points, and direction
  of "slider" semilandmark. The LEFT mouse button is used
  to select points.
}
\author{
  Dean Adams, Erik Otarola-Castillo, and Emma Sherratt
}
\references{
  Bookstein, F. J. 1991 Morphometric Tools for Landmark
  Data: Geometry and Biology. Cambridge University Press,
  New York.

  Bookstein, F. J. 1997 Landmark Methods for Forms without
  Landmarks: Morphometrics of Group Differences in Outline
  Shape. Medical Image Analysis 1(3):225-243.
}
\seealso{
  \code{\link{list.files}}
}
\keyword{digicurves}

