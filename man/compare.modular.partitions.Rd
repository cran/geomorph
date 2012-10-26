\name{compare.modular.partitions}
\alias{compare.modular.partitions}
\title{Compare modular signal to alternative landmark subsets}
\usage{
  compare.modular.partitions(A, land1, land2, iter = 999)
}
\arguments{
  \item{A}{An array (p x k x n) containing GPA-aligned
  coordinates for all specimens}

  \item{land1}{A list of landmark numbers for the first
  module}

  \item{land2}{A list of landmark numbers for the second
  module}

  \item{iter}{Number of iterations for significance
  testing}
}
\value{
  Function returns a list with the following components:
  \item{RV}{The estimate of morphological integration}
  \item{pvalue}{The significance level of the observed
  signal}
}
\description{
  Function quantifies the degree of morphological
  integration between two modules of Procrustes-aligned
  landmark coordinates and compares this to patterns found
  by randomly assigning landmarks into subsets
}
\details{
  The function quantifies the degree of morphological
  integration between two modules of shape data as defined
  by landmark coordinates, and compares this to modular
  signals found by randomly assigning landmarks to the two
  subsets. It is assumed that the landmarks have previously
  been aligned using Generalized Procrustes Analysis (GPA)
  [e.g., with \code{\link{gpagen}}]. The degree of
  morphological integration is quantified using the RV
  coefficient (Klingenberg 2009). The RV coefficient for
  the observed modular hypothesis is then compared to a
  distribution of values obtained by randomly assigning
  landmarks into subsets, with the restriction that the
  number of landmarks in each subset is identical to that
  observed in each of the original partitions. A
  significant modular signal is found when the observed RV
  coefficient is small relative to this distribution (see
  Klingenberg 2009). A histogram of coefficients obtained
  via resampling is presented, with the observed value
  designated by an arrow in the plot.

  The function currently evaluates morphological
  integration between two modules only.
}
\examples{
data(plethodon)
Y.gpa<-gpagen(plethodon$land)    #GPA-alignment

compare.modular.partitions(Y.gpa$coords,seq(1:5),seq(6:12),iter=99)
#Result implies that the skull and mandible are not independent modules
}
\author{
  Dean Adams
}
\references{
  Klingenberg, C. P. 2009. Morphometric integration and
  modularity in configurations of landmarks: tools for
  evaluating a priori hypotheses. Evol. Develop.
  11:405-421.
}
\keyword{compare.modular.partitions}

