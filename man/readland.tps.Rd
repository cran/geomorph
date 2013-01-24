\name{readland.tps}
\alias{readland.tps}
\title{Read landmark data from tps file}
\usage{
  readland.tps(file)
}
\arguments{
  \item{file}{A *.tps file containing two- or
  three-dimensional landmark data}
}
\value{
  Function returns a (p x k x n) array, where p is the
  number of landmark points, k is the number of landmark
  dimensions (2 or 3), and n is the number of specimens.
  The third dimension of this array contains names for each
  specimen, which are obtained from the image names in the
  *.tps file.
}
\description{
  Read *.tps file to obtain landmark coordinates
}
\details{
  This function reads a *.tps file containing two- or
  three-dimensional landmark coordinates.  Tps files are
  text files in one of the standard formats for geometric
  morphometrics (see Rohlf 2010).  Two-dimensional
  landmarks coordinates are designated by the identifier
  "LM=", while three-dimensional data are designated by
  "LM3=". Landmark coordinates are multiplied by their
  scale factor if this is provided for all specimens. If
  one or more specimens are missing the scale factor,
  landmarks are treated in their original units.  NOTE: all
  other information present in tps files (curves, comments,
  variables, radii, etc.)  is ignored.
}
\author{
  Dean Adams
}
\references{
  Rohlf, F. J. 2010. tpsRelw: Relative warps analysis.
  Version 1.49. Department of Ecology and Evolution, State
  University of New York at Stony Brook, Stony Brook, NY.
}
\keyword{readland.tps}

