\name{plotRefToTarget}
\alias{plotRefToTarget}
\title{Plot shape differences between a reference and target specimen}
\usage{
  plotRefToTarget(M1, M2,
    method = c("TPS", "vector", "points"), mag = 1,
    links = NULL, ...)
}
\arguments{
  \item{M1}{Matrix of landmark coordinates for the first
  (reference) specimen}

  \item{M2}{Matrix of landmark coordinates for the second
  (target) specimen}

  \item{method}{Method used to visualize shape difference;
  see below for details}

  \item{mag}{The desired magnification to be used when
  visualizing the shape difference (e.g., mag=2)}

  \item{links}{An optional matrix defining for links
  between landmarks}

  \item{...}{Additional parameters to be passed to 'plot'
  or 'plot3d'}
}
\description{
  Function plots shape differences between a reference and
  target specimen
}
\details{
  The function generates a plot of the shape differences of
  a target specimen relative to a reference specimen. The
  option "mag" allows the user to indicates the degree of
  magnification to be used when displaying the shape
  difference. If "method=TPS" a thin-plate spline
  deformation grid is generated.  If "method=vector" a plot
  showing the vector displacements between corresponding
  landmarks in the reference and target specimen is shown.
  If "method=points" a plot is displayed with the landmarks
  in the target (black) overlaying those of the reference
  (gray). Additionally, if a matrix of links is provided,
  the landmarks of the mean shape will be connected by
  lines.  The link matrix is an M x 2 matrix, where M is
  the desired number of links. Each row of the link matrix
  designates the two landmarks to be connected by that
  link. The function will plot either two- or
  three-dimensional data (note: for 3D, "method=TPS" will
  generate thin-plate spline deformations in the x-y and
  x-z planes).  This function combines numerous plotting
  functions found in Claude (2008).
}
\examples{
data(plethodon)
Y.gpa<-gpagen(plethodon$land)    #GPA-alignment
ref<-mshape(Y.gpa$coords)

# Differnt plotting options
plotRefToTarget(ref,Y.gpa$coords[,,39])

plotRefToTarget(ref,Y.gpa$coords[,,39],mag=3)   #magnify difference by 3X

plotRefToTarget(ref,Y.gpa$coords[,,39],method="vector")

plotRefToTarget(ref,Y.gpa$coords[,,39],method="points")

# Three dimensional data

data(scallops)
Y.gpa<-gpagen(A=scallops$coorddata, curves=scallops$curvslide, surfaces=scallops$surfslide)
ref<-mshape(Y.gpa$coords)
plotRefToTarget(ref,Y.gpa$coords[,,1],method="points")
}
\author{
  Dean Adams
}
\references{
  Claude, J. 2008. Morphometrics with R. Springer, New
  York.
}
\keyword{plotRefToTarget}

