\name{readmulti.nts}
\alias{readmulti.nts}
\title{Read landmark data from multiple nts files}
\usage{
  readmulti.nts(filelist)
}
\arguments{
  \item{filelist}{A list of names for the *.nts files to be
  read by the function. The names in the list require
  quotes (").}
}
\value{
  Function returns a (p x k x n) array, where p is the
  number of landmark points, k is the number of landmark
  dimensions (2 or 3), and n is the number of specimens.
  The third dimension of this array contains names for each
  specimen, which are obtained from the original file
  names.
}
\description{
  Read a list of names for several *.nts files to obtain
  landmark coordinates for a set of specimens
}
\details{
  This function reads a list containing the names of
  multiple *.nts files, where each contains the landmark
  coordinates for a single specimen. For these files, the
  number of variables (columns) of the data matrix will
  equal the number of dimensions of the landmark data (k=2
  or 3). When the function is called a dialog box is
  opened, from which the user may select multiple *.nts
  files. These are then read and concatenated into a single
  matrix for all specimens.
}
\author{
  Dean Adams
}
\keyword{readmulti.nts}

