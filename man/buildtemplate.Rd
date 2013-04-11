\name{buildtemplate}
\alias{buildtemplate}
\title{Build 3D surface template}
\usage{
  buildtemplate(specimen, fixed, surface.sliders)
}
\arguments{
  \item{specimen}{Name of matrix containing
  three-dimensional coordinates of a surface scan}

  \item{fixed}{numeric: the number of fixed template
  landmarks}

  \item{surface.sliders}{numeric: the number of template
  surface sliders desired}
}
\value{
  Function returns a matrix containing the x,y,z
  coordinates of the down sampled points, which can be used
  as a template for digitizing other surface scans using
  the function \code{\link{digitsurface}}. Additionally,
  the file 'template.txt' is generated, and an NTS file
  with the name of the specimen containing the digitized
  points for the template specimen (for use in subsequent
  morphometric analsyes).
}
\description{
  A function to build 3D template to extract 3D surface
  sliding semilandmarks from all specimens
}
\details{
  Function buildtemplate constructs a template surface with
  which to down sample point clouds of specimens to be used
  in 3d shape analysis. builddtemplate allows users to
  choose a predetermined number of points with which to
  represent the structure of interest as sliding surface
  semilandmarks. Template surface is used in analyses of
  surface semi-landmarks outline in Gunz et al.
  (2005:90-92) and Mitteroecker and Gunz (2009:242).
  Landmark points are first digitized by the user (see
  Digitizing subsection), then a roughly equidistant set of
  predetermined number of points are automatically chosen.
  The now down-sampled mesh is used as a "template" to
  extract a mesh of similarly numbered points for analysis.
  \subsection{Digitizing}{Digitizing using buildtemplate is
  interactive between landmark selection using a mouse (see
  below for instructions), and the R console. Once a point
  is selected, the user is asked if the system should keep
  or discard the selection #'(y/n). If "y", the user is
  asked to continue to select the next landmark.If "n" the
  removes the last chosen landmark, and the user is askesd
  to select it again. This can be repeated until the user
  is comfortable with the landmark chosen.

  NOTE: To ensure a strong match between the scan and the
  template, it is recommended that a reasonable number of
  fixed points be used.

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
}
\author{
  \href{http://www.people.fas.harvard.edu/~eotarolacastillo}{Erik
  Otarola-Castillo} and
  \href{http://www.public.iastate.edu/~dcadams}{Dean
  Adams}.
}
\references{
  Gunz P, Mitteroecker P, & Bookstein FJ (2005)
  Semilandmarks in Three Dimensions. Modern Morphometrics
  in Physical Anthropology, ed Slice DE (Springer-Verlag,
  New York), pp 73-98.

  Mitteroecker P & Gunz P (2009) Advances in Geometric
  Morphometrics. Evolutionary Biology 36(2):235-247.
}
\keyword{buildtemplate}
\keyword{template}

