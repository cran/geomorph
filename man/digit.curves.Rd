\name{digit.curves}
\alias{digit.curves}
\title{Select points to "slide" along three-dimensional curves.}
\usage{
  digit.curves(n, curves)
}
\arguments{
  \item{n}{Number of total "fixed" and curve sliding
  landmarks.}

  \item{curves}{Number of landmarks to slide along curves.}
}
\value{
  Function returns a matrix containing the landmark adress
  of the curve sliders, indicating the points between which
  the selected point will "slide". In addition, the
  function returns a .csv file to be used by
  \code{\link{gpagen}} during GPA.
}
\description{
  A function to select points which will "slide" along
  three-dimensional curves.
}
\details{
  This function helps user select points on the created
  template to "slide" along curves lacking known landmarks
  (see Bookstein 1991:376-382, 1997 for algorithm details).
  Each sliding semi-landmark (sliders) will slide between
  two designated points, along a line tangent to the
  specified curvature. template.txt file must be in current
  working directory. \enumerate{ \item Select the first
  point between which semi-landmark will "slide" \item
  Select sliding point, \item Select point along which
  sliding trajectory will end. Screen will show lines
  connecting the three points, and will highlight the
  sliding semilandmark in red. } \subsection{Digitizing}{
  Digitizing using buildtemplate is interactive between
  landmark selection using a mouse (see below for
  instructions), and the R console. Once a point is
  selected, the user is asked if the system should keep or
  discard the selection #'(y/n). If "y", the user is asked
  to continue to select the next landmark.If "n" the
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
  to adjust mesh perspective.}

  Mac mouse settings need adjustment to acquire zooming
  functions of the "middle/third" mouse button/scroller.
  Some instructions can be found at
  \url{http://blog.boastr.net/}. In addition, freeware
  providing "Middle Click" functionality is available at
  \url{http://magicprefs.com/} for "magicmice" now standard
  on many Macintosh machines. }
}
\author{
  \href{http://www.people.fas.harvard.edu/~eotarolacastillo}{Erik
  Otarola-Castillo} and
  \href{http://www.public.iastate.edu/~dcadams}{Dean
  Adams}.
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
  \code{\link{digitsurface}}, \code{\link{gpagen}}
}
\keyword{digicurves}

