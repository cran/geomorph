\name{plotTangentSpace}
\alias{plotTangentSpace}
\title{Plot specimens in tangent space}
\usage{
plotTangentSpace(A, axis1 = 1, axis2 = 2, warpgrids = TRUE,
  label = FALSE, groups = NULL)
}
\arguments{
  \item{A}{An array (p x k x n) containing landmark
  coordinates for a set of specimens}

  \item{warpgrids}{A logical value indicating whether
  deformation grids for shapes along PC1 should be
  displayed}

  \item{axis1}{A value indicating which PC axis should be
  displayed as the X-axis (default = PC1)}

  \item{axis2}{A value indicating which PC axis should be
  displayed as the X-axis (default = PC2)}

  \item{label}{A logical value indicating whether labels
  for each specimen should be displayed}

  \item{groups}{An optional factor vector specifying group
  identity for each specimen}
}
\description{
Function plots a set of Procrustes-aligned specimens in
tangent space along their principal axes
}
\details{
The function performs a principal components analysis of
shape variation and plots two dimensions of tangent space
for a set of Procrustes-aligned specimens (default is PC1
vs. PC2). The percent variation along each PC-axis is
returned. Additionally (and optionally), deformation grids
can be requested, which display the shape of specimens at
the ends of the range of variability along PC1. The
function returns a table summarizing the percent variation
explained by each pc axis, and the set of principal
component scores for all specimens. If groups are provided,
specimens from each group are plotted using distinct colors
based on the order in which the groups are found in the
dataset, and using R's standard color palette: black, red,
green, blue, cyan, magenta, yellow, and gray.
}
\examples{
data(plethodon)
Y.gpa<-gpagen(plethodon$land)    #GPA-alignment
ref<-mshape(Y.gpa$coords)

plotTangentSpace(Y.gpa$coords, groups = paste(plethodon$species, plethodon$site))
}
\author{
Dean Adams
}
\keyword{plotTangentSpace}

