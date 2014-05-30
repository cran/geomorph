\name{readland.tps}
\alias{readland.tps}
\title{Read landmark data from tps file}
\usage{
readland.tps(file, specID = c("None", "ID", "imageID"), warnmsg = T)
}
\arguments{
  \item{file}{A *.tps file containing two- or
  three-dimensional landmark data}

  \item{specID}{a character specifying whether to extract
  the specimen ID names from the ID or IMAGE lines (default
  is "None").}

  \item{warnmsg}{A logical value stating whether warnings
  should be printed}
}
\value{
Function returns a (p x k x n) array, where p is the number
of landmark points, k is the number of landmark dimensions
(2 or 3), and n is the number of specimens. The third
dimension of this array contains names for each specimen,
which are obtained from the image names in the *.tps file.
}
\description{
Read *.tps file to obtain landmark coordinates
}
\details{
This function reads a *.tps file containing two- or
three-dimensional landmark coordinates.  Tps files are text
files in one of the standard formats for geometric
morphometrics (see Rohlf 2010).  Two-dimensional landmarks
coordinates are designated by the identifier "LM=", while
three-dimensional data are designated by "LM3=". Landmark
coordinates are multiplied by their scale factor if this is
provided for all specimens. If one or more specimens are
missing the scale factor, landmarks are treated in their
original units.

Missing data may be present in the file. In this case, they
must be designated by 'NA'. The positions of missing
landmarks may then be estimated using estimate.missing.

The user may specify whether specimen names are to be
extracted from the 'ID=' field or 'IMAGE=' field and
included in the resulting 3D array. e.g., for 'ID=' use
(file, specID = "ID") and for 'IMAGE=' use (file, specID =
"imageID"). The default is specID="None".

NOTE: At present, all other information that can be
contained in tps files (curves, comments, variables, radii,
etc.)  is ignored.
}
\author{
Dean Adams & Emma Sherratt
}
\references{
Rohlf, F. J. 2010. tpsRelw: Relative warps analysis.
Version 1.49. Department of Ecology and Evolution, State
University of New York at Stony Brook, Stony Brook, NY.
}
\keyword{IO}

