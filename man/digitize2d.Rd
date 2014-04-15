\name{digitize2d}
\alias{digitize2d}
\title{Digitize 2D landmarks.}
\usage{
digitize2d(file, nlandmarks, scale, verbose = TRUE)
}
\arguments{
  \item{file}{Name of jpeg file to be digitized. File names
  can be written in manually, including paths, or obtained
  using directory/file manipulation functions e.g.,
  \code{\link{list.files}}}

  \item{nlandmarks}{Number of landmarks to be digitized.}

  \item{scale}{Length of scale placed in image.}

  \item{verbose}{logical. User decides whether to digitize
  in verbose or silent format (see details), default is
  verbose}
}
\value{
Function returns an n-x-2 .nts file containing the 2d
coordinates of digitized landmarks
}
\description{
An interactive function to digitize two-dimensional(2D)
landmarks from .jpg files.
}
\details{
Function for digitizing 2D landmarks on specimen images
(.jpg). "nlandmarks" is the number of landmark points to be
digitized by the user. Landmarks should include "true"
landmarks and semi-landmarks to be "sliders". For best
results, digitizing sequence should proceed by selecting
all true landmark points first, followed by selection of
sliding semi-landmarks. Use function
\code{\link{define.sliders.2d}} to select sliding
semi-landmarks.

\subsection{Digitizing}{ Digitizing landmarks involves
landmark selection using a mouse in the plot window, using
the LEFT mouse button (or regular button for Mac users):
\enumerate{ \item Digitize the scale bar by selecting the
two end points (single click for start and end), \item
Digitize each landmark with single click and the landmark
is shown in red, When selection of n landmarks is
completed, an ".nts" file is created in working directory
using the specimen name, adding ".nts" as a suffix. } If
verbose = TRUE, digitizing is interactive between landmark
selection using a mouse and the R console. Once a landmark
is selected, the user is asked if the system should keep or
discard the selection (y/n). If "y", the user is asked to
continue to select the next landmark. If "n", the user is
asked to select it again. This can be repeated until the
user is comfortable with the landmark chosen. verbose =
FALSE is silent, and digitizing of each landmark is
continuous and uninterupted until all landmarks are chosen.
Landmark coordinates are returned scaled.}
}
\author{
Erik Otarola-Castillo and Emma Sherratt
}
\seealso{
\code{\link{list.files}}
}
\keyword{digitizing}

