\name{digit.fixed}
\alias{digit.fixed}
\title{Digitize fixed 3D landmarks only.}
\usage{
  digit.fixed(spec, fixed, index = FALSE, ptsize = 1)
}
\arguments{
  \item{spec}{Name of surface file, as either an object of
  class shape3d/mesh3d, or matrix of three-dimensional
  vertex coordinates.}

  \item{fixed}{numeric: The number of fixed template
  landmarks}

  \item{index}{logical: Whether selected landmark addresses
  should be returned}

  \item{ptsize}{numeric: Size to plot the mesh points
  (vertices), e.g. 0.1 for dense meshes, 3 for sparse
  meshes}
}
\description{
  A function to digitize only fixed landmarks. Input for
  the function is either a matrix of 3D coordinates or a
  mesh3d object as obtained from read.ply().
}
\details{
  Function to digitize landmarks on 3D surfaces or
  pointclouds. "n" Landmark points are selected by user. No
  template is used to select surface sliding semi
  landmarks. Select points to be digitized, using the RIGHT
  mouse button. The LEFT mouse button is used to ROTATE
  mesh, and the mouse SCROLLER is used to zoom in and out.
  When selection of n landmarks is completed, an ".nts"
  file is created in working directory using the specimen
  name, adding "fixedlmcoords.nts" as a suffix.

  NOTE: To ensure a strong match between the scan and the
  template, it is recommended that a reasonable number of
  fixed points be used.

  \subsection{Digitizing}{ Digitizing using buildtemplate
  is interactive between landmark selection using a mouse
  (see below for instructions), and the R console. Once a
  point is selected, the user is asked if the system should
  keep or discard the selection #'(y/n). If "y", the user
  is asked to continue to select the next landmark. If "n"
  the removes the last chosen landmark, and the user is
  asked to select it again. This can be repeated until the
  user is comfortable with the landmark chosen.

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
  Erik Otarola-Castillo & Emma Sherratt
}
\keyword{digifix}

