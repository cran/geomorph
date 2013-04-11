\name{physignal}
\alias{physignal}
\title{Assessing phylogenetic signal in morphometric data}
\usage{
  physignal(phy, A, iter = 999)
}
\arguments{
  \item{phy}{A phylogenetic tree of type 'phylo'}

  \item{A}{An array (p x k x n) containing GPA-aligned
  coordinates for a set of specimens}

  \item{iter}{Number of iterations for significance
  testing}
}
\value{
  Function returns a list with the following components:
  \item{phy.signal}{The estimate of phylogenetic signal}
  \item{pvalue}{The significance level of the observed
  signal}
}
\description{
  Function calculates the degree of phylogenetic signal
  from a set of Procrustes-aligned specimens
}
\details{
  The function estimates the degree of phylogenetic signal
  present in shape data for a given phylogeny based on a
  Brownian motion model of evolution. It is assumed that
  the landmarks have previously been aligned using
  Generalized Procrustes Analysis (GPA) [e.g., with
  \code{\link{gpagen}}]. To assess phylogenetic signal
  under alternative models, first consider branch-length
  transformations of the phylogeny.  Phylogenetic signal is
  estimated as the sum of squared changes (SSC) in shape
  along all branches of the phylogeny (Klingenberg and
  Gidasqewski 2010). Significance testing is found by
  permuting the shape data among the tips of the phylogeny.
  A plot of the specimens in tangent space with the
  phylogeny superimposed is included, and ancestral states
  are returned. Note that the method can be slow as
  ancestral states must be estimated for every iteration.
}
\examples{
data(plethspecies)
Y.gpa<-gpagen(plethspecies$land)    #GPA-alignment

physignal(plethspecies$phy,Y.gpa$coords,iter=5)
}
\author{
  Dean Adams
}
\references{
  Klingenberg, C. P., and N. A. Gidaszewski. 2010. Testing
  and quantifying phylogenetic signals and homoplasy in
  morphometric data. Syst. Biol. 59:245-261.
}
\keyword{physignal}

