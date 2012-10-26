\name{plotGMPhyloMorphoSpace}
\alias{plotGMPhyloMorphoSpace}
\title{Plot phylogenetic tree and specimens in tangent space}
\usage{
  plotGMPhyloMorphoSpace(phy, A, labels = TRUE,
    ancStates = T)
}
\arguments{
  \item{phy}{A phylogenetic tree of type 'phylo'}

  \item{A}{An array (p x k x n) containing landmark
  coordinates for a set of specimens}

  \item{labels}{A logical value indicating whether taxa
  labels should be included}

  \item{ancStates}{A logical value indicating whether
  ancestral state values should be returned}
}
\description{
  Function plots a phylogenetic tree and a set of
  Procrustes-aligned specimens in tangent space
}
\details{
  The function creates a plot of the first two dimensions
  of tangent space for a set of Procrustes-aligned
  specimens. The phylogenetic tree for these specimens is
  superimposed in this plot revealing how shape evolves
  (e.g., Rohlf 2002; Klingenberg and Gidaszewski 2010). The
  plot also displays the ancestral states for each node of
  the phylogenetic tree, whose state values can optionally
  be returned.
}
\examples{
data(plethspecies)
Y.gpa<-gpagen(plethspecies$land)    #GPA-alignment

plotGMPhyloMorphoSpace(plethspecies$phy,Y.gpa$coords)
}
\author{
  Dean Adams
}
\references{
  Klingenberg, C. P., and N. A. Gidaszewski. 2010. Testing
  and quantifying phylogenetic signals and homoplasy in
  morphometric data. Syst. Biol. 59:245-261.

  Rohlf, F. J. 2002. Geometric morphometrics and phylogeny.
  Pp. 175'193 in N. Macleod, and P. Forey, eds. Morphology,
  shape, and phylogeny. Taylor & Francis, London.
}
\keyword{plotGMPhyloMorphoSpace}

