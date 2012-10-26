\name{read.ply}
\alias{read.ply}
\title{Read landmark data from ply files}
\usage{
  read.ply(file)
}
\arguments{
  \item{file}{A ply file in NextEngine or David scanner
  format}
}
\value{
  Function returns a list with the following components:
  \item{coords}{The x,y,z coordinates of the ply surface}
  \item{polygons}{A set of polygons connecting triplets of
  coordinates}
}
\description{
  Read ply files to obtain landmark coordinates
}
\details{
  This function reads three-dimensional surface data in the
  form of a single ply file (in either NextEngine or David
  scanner format). The landmarks of this surface may then
  be used to digitize three-dimensional points, and
  semilandmarks on curves and surfaces.
}
\author{
  Dean Adams
}
\keyword{read.ply}

