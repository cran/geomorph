\name{compare.evol.rates}
\alias{compare.evol.rates}
\title{Comparing rates of shape evolution on phylogenies}
\usage{
compare.evol.rates(phy, A, gp, iter = 999)
}
\arguments{
  \item{phy}{A phylogenetic tree of type 'phylo'}

  \item{A}{An array (p x k x n) containing GPA-aligned
  coordinates for a set of specimens}

  \item{gp}{A factor array designating group membership}

  \item{iter}{Number of iterations for significance
  testing}
}
\value{
Function returns a list with the following components:
\item{sigmad.all}{The phylogenetic evolutionary rate for
all species on the phylogeny} \item{sigmad.ratio}{The ratio
of maximum to minimum evolutionary rates} \item{pvalue}{The
significance level of the observed ratio}
}
\description{
Function calculates rates of shape evolution for two or
more groups of species on a phylogeny from a set of
Procrustes-aligned specimens
}
\details{
The function compares rates of morphological evolution for
two or more groups of species on a phylogeny, under a
Brownian motion model of evolution. It is assumed that the
landmarks have previously been aligned using Generalized
Procrustes Analysis (GPA) [e.g., with
\code{\link{gpagen}}]. The approach is based on the
distances between species in morphospace after phylogenetic
transformation (Adams 2014). From the data the rate of
shape evolution for each group is calculated, and a ratio
of rates is obtained. If three or more groups of species
are used, the ratio of the maximum to minimum rate is used
as a test statistic (see Adams 2014). Significance testing
is accomplished by phylogenetic simulation in which tips
data are obtained under Brownian motion using a single
evolutionary rate for all species on the phylogeny. A
histogram of evolutionary rate ratios obtained via
phylogenetic simulation is presented, with the observed
value designated by an arrow in the plot.
}
\examples{
data(plethspecies)
Y.gpa<-gpagen(plethspecies$land)    #GPA-alignment

 gp.end<-factor(c(0,0,1,0,0,1,1,0,0))  #endangered species vs. rest
 names(gp.end)<-plethspecies$phy$tip

compare.evol.rates(plethspecies$phy,Y.gpa$coords,gp=gp.end,iter=49)
}
\author{
Dean Adams
}
\references{
Adams, D.C. 2014. Quantifying and comparing phylogenetic
evolutionary rates for shape and other high-dimensional
phenotypic data. Syst. Biol. (Accepted).
}
\keyword{compare.phylo.shape.rate}

