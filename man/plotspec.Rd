\name{plotspec}
\alias{plotspec}
\title{Plot 3D specimen, fixed landmarks and surface semilandmarks.}
\usage{
plotspec(spec, digitspec, fixed, ptsize = 1)
}
\arguments{
  \item{spec}{Name of surface file, as either an object of
  class shape3d/mesh3d, or matrix of three-dimensional
  vertex coordinates.}

  \item{digitspec}{Name of data matrix containing
  three-dimensional fixed and/or surface sliding
  coordinates.}

  \item{fixed}{numeric: The number of fixed template
  landmarks}

  \item{ptsize}{numeric: Size to plot the mesh points
  (vertices), e.g. 0.1 for dense meshes, 3 for sparse
  meshes}
}
\description{
A function to plot three dimensional (3D) specimen along
with its landmarks.
}
\details{
Function to plot 3D specimens along with their digitized
"fixed" landmarks and semilandmarks "surface sliders" and
"curve sliders". If specimen is a 3D surface (ply file),
mesh is plotted.
}
\examples{
data(scallopPLY)
ply <- scallopPLY$ply
digitdat <- scallopPLY$coords
plotspec(spec=ply,digitspec=digitdat,fixed=16)
}
\author{
Erik Otarola-Castillo & Emma Sherratt
}
\keyword{plotspec}

