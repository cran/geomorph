\name{read.ply}
\alias{read.ply}
\title{Read mesh data (vertices and faces) from ply files}
\usage{
  read.ply(file, ShowSpecimen = TRUE)
}
\arguments{
  \item{file}{An ASCII ply file}

  \item{ShowSpecimen}{logical: A logical value indicating
  whether or not the ply file should be displayed}
}
\value{
  Function returns the following components:
  \item{mesh3d}{list of class mesh3d- see rgl for details}
}
\description{
  Read ply files to obtain landmark coordinates
}
\details{
  This function reads three-dimensional surface data in the
  form of a single ply file (Polygon File Format; ASCII
  format only, from 3D scanners such as NextEngine and
  David scanners). Vertices of the surface may then be used
  to digitize three-dimensional points, and semilandmarks
  on curves and surfaces. The function opens the ply file
  and plots the mesh, with faces rendered if file contains
  faces, and colored if the file contains vertex color.
}
\author{
  Dean Adams & Emma Sherratt
}
\keyword{read.ply}

