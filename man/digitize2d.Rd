\name{digitize2d}
\alias{digitize2d}
\title{Digitize 2d landmarks.}
\usage{
  digitize2d(file, nlandmarks, scale)
}
\arguments{
  \item{file}{Name of jpeg file in working directory (can
  include path) to be digitized. File names can be written
  in manually, including paths, or obtained using
  directory/file manipulation functions e.g.,
  \code{\link{list.files}}}

  \item{nlandmarks}{Number of landmarks to be digitized.}

  \item{scale}{Desired length of scale placed on image.}
}
\value{
  Function returns an n-x-2 .nts file containing the 2d
  coordinates of digitized landmarks

  Returns text file with the distance in scale
}
\description{
  A function to digitize 2d landmarks from .jpg files.
}
\details{
  digitize2d is a function to digitize 2d landmarks on
  specimen images (.jpg). "nlandmarks" is the number of
  landmark points to be digitized by the user. Landmarks
  should include "true" landmarks and semilandmarks to be
  "slid". For best results, digitizing sequence should
  proceed by selecting all true landmark points first,
  followed by selection of sliding semi-landmarks. Use
  function "curves2d" to select sliding semilandmarks.
  After choosing image to digitize, users digitize scale
  within image (LEFT mouse button). Then landmark points
  can be digitized using the LEFT mouse button. When
  selection of n landmarks is completed, an ".nts" file is
  created in working directory using the specimen name,
  adding "2dcoords.nts" as a suffix.
}
\author{
  Erik Otarola-Castillo and Dean Adams
}
\seealso{
  \code{\link{list.files}}
}
\keyword{digitize2d}

