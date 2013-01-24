\name{read.morphologika}
\alias{read.morphologika}
\title{Read 3D landmark data from Morphologika files}
\usage{
  read.morphologika(file, plot = FALSE)
}
\arguments{
  \item{file}{A morphologika text file. File name can be
  written in manually, including path, or obtained using
  directory/file manipulation functions e.g.,
  \code{\link{list.files}}}

  \item{plot}{Logical should (unaligned) specimens be
  plotted? Defaults to FALSE.}
}
\value{
  Function returns a list with the following components:
  \item{coords}{If Morphologika headers "[labels]" and
  "[labelvalues]" are not present in the file, function
  only returns a (p x k x n) array, where p is the number
  of landmark points, k is the number of landmark
  dimensions, and n is the number of specimens}
  \item{dataframe}{If Morphologika headers "[labels]" and
  "[labelvalues]" are present read.morphologika returns the
  above p x k x n array and a dataframe containing specimen
  specific information stored in file}
}
\description{
  Read Morphologika files to obtain 3D landmark coordinates
  and specimen information
}
\details{
  This function reads the commonly used Morphologika file
  format for 2D and 3D landmarks. Landmark coordinates and
  specimen information may then be used to conduct GPA
  using \code{\link{gpagen}}, and select three-dimensional
  semilandmarks on curves (if present) using
  \code{\link{digit.curves}}. If the Morphologika headers
  "[labels]" and "[labelvalues]" are present, then a data
  matrix containing all individual specimen information is
  returned. If "[labels]" and "[labelvalues]" are not
  present in the file, then only landmark coordinates are
  returned.
}
\author{
  Erik Otarola-Castillo and Dean Adams
}
\seealso{
  \code{\link{list.files}}
}
\keyword{read.morphologika}

