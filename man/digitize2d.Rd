\name{digitize2d}
\alias{digitize2d}
\title{Digitize 2D landmarks.}
\usage{
digitize2d(spec, nlandmarks, scale)
}
\arguments{
  \item{spec}{Name of jpeg file to be digitized. File names
  can be written in manually, including paths, or obtained
  using directory/file manipulation functions e.g.,
  \code{\link{list.files}}}

  \item{nlandmarks}{Number of landmarks to be digitized.}

  \item{scale}{Length of scale placed in image.}
}
\value{
Function writes to the working directory an NTS file with
the name of the specimen and .nts suffix containing the
landmark coordinates. Function also returns to console an n
x 2 matrix containing the x,y coordinates of the digitized
landmarks.
}
\description{
An interactive function to digitize two-dimensional(2D)
landmarks from .jpg files.
}
\details{
Function for digitizing 2D landmarks on specimen images
(.jpg). "nlandmarks" are digitized by the user (see below).
The landmarks are "fixed" (traditional landmarks).  They
can be later designated as "curve sliders" (semilandmarks,
that will "slide" along curves lacking known landmarks (see
Bookstein 1991:376-382, 1997 for algorithm details) if
required. These landmarks are defined as "sliders" using
function \code{\link{define.sliders.2d}} or with similar
format matrix made outside R. Digitizing sequence should
proceed by selecting all "fixed" landmark points first,
followed by selection of "curve sliders".

\subsection{Digitizing}{ Digitizing landmarks involves
landmark selection using a mouse in the plot window, using
the LEFT mouse button (or regular button for Mac users):
\enumerate{ \item Digitize the scale bar by selecting the
two end points (single click for start and end), NOTE: no
dot is printed on the image on the first click; only when
both ends have been digitized.  \item Digitize each
landmark with single click and the landmark is shown in
red, When selection of n landmarks is completed, an ".nts"
file is created in working directory using the specimen
name, adding ".nts" as a suffix. } }
}
\author{
Erik Otarola-Castillo & Emma Sherratt
}
\seealso{
\code{\link{list.files}},
}
\keyword{digitize2d}

