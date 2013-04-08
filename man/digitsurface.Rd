\name{digitsurface}
\alias{digitsurface}
\title{Digitize 3D fixed landmarks and surface semilandmarks.}
\usage{
  digitsurface(specimen, fixed)
}
\arguments{
  \item{specimen}{Name of data matrix in working directory
  containing three-dimensional landmark coordinates.}

  \item{fixed}{numeric: the number of fixed template
  landmarks}
}
\description{
  A function to digitize three dimensional fixed landmarks
  and surface semilandmarks.
}
\details{
  Function to digitize landmarks on 3D pointclouds. "n"
  Landmark points are selected by user akin to landmarks
  selected to construct template using function
  \code{\link{buildtemplate}}. Following selection of
  points, function digitsurface finds surface semilandmarks
  following algorithm outlined in Gunz et al. (2005:90-92)
  and Mitteroecker and Gunz (2009:242). digitsurface finds
  the same number of surface semi-landmarks as the template
  (created by buildtemplate) by downsampling scanned mesh,
  registering template with current specimen via GPA. A
  nearest neighbor algorithm is used to match template
  surface landmarks to current specimen's.
  \subsection{Digitizing}{Digitizing using buildtemplate is
  interactive between landmark selection using a mouse (see
  below for instructions), and the R console. Once a point
  is selected, the user is asked if the system should keep
  or discard the selection #'(y/n). If "y", the user is
  asked to continue to select the next landmark.If "n" the
  removes the last chosen landmark, and the user is askesd
  to select it again. This can be repeated until the user
  is comfortable with the landmark chosen.

  To digitize with a standard 3-button (PC) buildtemplate
  uses: \enumerate{ \item the RIGHT mouse button (primary)
  to select points to be digitized First, \item the LEFT
  mouse button (secondary) is used to rotate mesh, \item
  the mouse SCROLLER (third/middle) is used to zoom in and
  out. } Note: Digitizing functions on MACINTOSH computers
  using a standard 3-button mice works as specified. Macs
  using platform specific single button mice: \enumerate{
  \item press button to rotate 3D mesh, \item press button
  while pressing CONTROL key to select points to be
  digitized, \item press button while pressing OPTION key
  to adjust mesh perspective. } Mac mouse settings need
  adjustment to acquire zooming functions of the
  "middle/third" mouse button/scroller. Some instructions
  can be found at \url{http://blog.boastr.net/}. In
  addition, freeware providing "Middle Click" functionality
  is available at \url{http://magicprefs.com/} for
  "magicmice" now standard on many Macintosh machines. }
  When completed, an ".nts" file is created in working
  directory using the specimen name, adding "coords.nts" as
  a suffix. This file contains the specimen coordinates to
  be used by GPA.
}
\author{
  \href{http://www.people.fas.harvard.edu/~eotarolacastillo}{Erik
  Otarola-Castillo} and
  \href{http://www.public.iastate.edu/~dcadams}{Dean Adams}
}
\references{
  Gunz P, Mitteroecker P, & Bookstein FJ (2005)
  Semilandmarks in Three Dimensions. Modern Morphometrics
  in Physical Anthropology, ed Slice DE (Springer-Verlag,
  New York), pp 73-98.

  Mitteroecker P & Gunz P (2009) Advances in Geometric
  Morphometrics. Evolutionary Biology 36(2):235-247.
}
\seealso{
  \code{\link{buildtemplate}}
}
\keyword{digitsurface}

