\name{plotTangentSpace}
\alias{plotTangentSpace}
\title{Plot specimens in tangent space}
\usage{
  plotTangentSpace(A, warpgrids = TRUE)
}
\arguments{
  \item{A}{An array (p x k x n) containing landmark
  coordinates for a set of specimens}

  \item{warpgrids}{A logical value indicating whether
  deformation grids for shapes along PC1 should be
  displayed}
}
\description{
  Function plots a set of Procrustes-aligned specimens in
  tangent space along their principal axes
}
\details{
  The function performs a principal compoments analysis of
  shape variation and plots the first two dimensions of
  tangent space for a set of Procrustes-aligned specimens.
  The percent variation along each PC-axis is returned.
  Additionally (and optionally), deformation grids can be
  requested, which display the shape of specimens at the
  ends of the range of variability along PC1.
}
\examples{
data(plethodon)
Y.gpa<-gpagen(plethodon$land)    #GPA-alignment
ref<-mshape(Y.gpa$coords)

plotTangentSpace(Y.gpa$coords)
}
\author{
  Dean Adams
}
\keyword{plotTangentSpace}

