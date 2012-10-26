\name{read.vrml}
\alias{read.vrml}
\title{Read landmark data from .vrml files}
\usage{
  read.vrml(file, plotspec = TRUE, plottri = TRUE,
    write.nts = FALSE)
}
\arguments{
  \item{file}{A .vrml file with coordinates in a
  "centralized" format, or in "stitched" format. File names
  can be written in manually, including paths, or obtained
  using directory/file manipulation functions e.g.,
  \code{\link{list.files}}}

  \item{plotspec}{Logical should object be plotted.
  Defaults to TRUE.}

  \item{plottri}{Logical should triangles be plotted.
  Defaults to TRUE.}

  \item{write.nts}{Logical should .nts file be created.
  Defaults to FALSE.}
}
\value{
  Function returns a list with the following components:
  \item{coords}{The x,y,z coordinates of the .vrml surface}
  \item{triangles}{Triangle facet connections between
  coordinate points (if available)}
}
\description{
  Read vrml files (Virtual Reality Modeling Language) to
  obtain landmark coordinates and triangulations
}
\details{
  This function reads three-dimensional surface data in the
  form of a single vrml file (Virtual Reality Modeling
  Language). The landmarks of this surface may then be used
  to digitize three-dimensional points, and semilandmarks
  on curves and surfaces. .vrml files are stored either as
  centralized data, where 3D scanned object information
  (i.e., coordinates, triangles, color, and surface) are
  stored within single data blocks for the complete object.
  This version of read.vrml will import .wrl files written
  in "utf8" and "ascii" format. Mesh triangle facets will
  be imported if present as coordinate connections.
  Argument plotspec allows users to plot file. This is
  helpful to help inspect object for potential errors.
  Argument plottri plots triangle facets if present.
  write.nts provides users with the option of writing a
  .nts coordinate file.
}
\author{
  Erik Otarola-Castillo and Dean Adams
}
\seealso{
  \code{\link{list.files}}
}
\keyword{read.vrml}

