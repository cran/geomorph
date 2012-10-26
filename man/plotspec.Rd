\name{plotspec}
\alias{plotspec}
\title{Plot 3D specimen, fixed landmarks and surface semilandmarks.}
\usage{
  plotspec(specimen, digitspec, fixed)
}
\arguments{
  \item{specimen}{Name of data matrix containing "raw"
  three-dimensional landmark coordinates.}

  \item{digitspec}{Name of data matrix containing
  three-dimensional fixed and/or surface sliding
  coordinates.}

  \item{fixed}{Numeric: the number of fixed template
  landmarks}
}
\description{
  A function to plot three dimensional specimen along with
  its fixed landmarks and surface semilandmarks.
}
\details{
  Function to plot 3D specimens along with their digitized
  "fixed" and surface sliding semilandmarks.
}
\examples{
data(Specimen4Raw)
rawdat<-as.matrix(Specimen4Raw)
data(scallops)
digitdat<-scallops$coorddata[,,4]
plotspec(specimen=rawdat,digitspec=digitdat,fixed=16)
}
\author{
  \href{http://www.people.fas.harvard.edu/~eotarolacastillo}{Erik
  Otarola-Castillo} and
  \href{http://www.public.iastate.edu/~dcadams}{Dean Adams}
}
\keyword{plotspec}

