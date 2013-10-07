\name{readland.nts}
\alias{readland.nts}
\title{Read landmark data matrix from nts file}
\usage{
readland.nts(file)
}
\arguments{
  \item{file}{A *.nts file containing two- or
  three-dimensional landmark data}
}
\value{
Function returns a (p x k x n) array, where p is the number
of landmark points, k is the number of landmark dimensions
(2 or 3), and n is the number of specimens. The third
dimension of this array contains names for each specimen,
which are obtained from the image names in the *.nts file.
}
\description{
Read *.nts file to obtain landmark coordinates for a set of
specimens
}
\details{
Function reads a *.nts file containing a matrix of two- or
three-dimensional landmark coordinates.  NTS files are text
files in one of the standard formats for geometric
morphometrics (see Rohlf 2012).  The parameter line
contains 5 or 6 elements, and must begin with a "1" to
designate a rectangular matrix. The second and third values
designate how many specimens (n) and how many total
variables (p x k) are in the data matrix. The fourth value
is a "0" if the data matrix is complete and a "1" if there
are missing values. If missing values are present, the '1'
is followed by the arbitrary numeric code used to represent
missing values (e.g., -999). These values will be replaced
with "-999" in the output array. The final value of the
parameter line denotes the dimensionality of the landmarks
(2,3) and begins with "DIM=". If specimen and variable
labels are included, these are designated placing an "L"
immediately following the specimen or variable values in
the parameter file. The labels then precede the data
matrix.

Function is for *.nts file containing landmark coordinates
for multiple specimens. Function can read *.dta files in
the nts format written by Landmark Editor
\url{http://graphics.idav.ucdavis.edu/research/projects/EvoMorph}.
}
\author{
Dean Adams & Emma Sherratt
}
\references{
Rohlf, F. J. 2012 NTSYSpc: Numerical taxonomy and
multivariate analysis system. Version 2.2. Exeter Software,
New York.
}
\keyword{readland.nts}

