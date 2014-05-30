\name{procD.pgls}
\alias{procD.pgls}
\title{Phylogenetic ANOVA/regression for shape data}
\usage{
procD.pgls(f1, phy, iter = 999)
}
\arguments{
  \item{f1}{A formula for the linear model (e.g., y~x1+x2)}

  \item{phy}{A phylogenetic tree of {class phylo} - see
  \code{\link[ape]{read.tree}} in library ape}

  \item{iter}{Number of iterations for significance
  testing}
}
\value{
Function returns an ANOVA table of statistical results for
all factors: df (for each factor), SS, MS, F ratio, Prand,
and Rsquare.
}
\description{
Function performs Procrustes ANOVA in a phylogenetic
framework and uses permutation procedures to assess
statistical hypotheses describing patterns of shape
variation and covariation for a set of Procrustes-aligned
coordinates
}
\details{
The function performs ANOVA and regression models in a
phylogenetic context under a Brownian motion model of
evolution, in a manner that can accommodate
high-dimensional datasets. The approach is derived from the
statistical equivalency between parametric methods
utilizing covariance matrices and methods based on distance
matrices (Adams 2014). Data input is specified by a formula
(e.g., y~X), where 'y' specifies the response variables
(shape data), and 'X' contains one or more independent
variables (discrete or continuous). The response matrix 'y'
must be in the form of a two-dimensional data matrix of
dimension (n x [p x k]), rather than a 3D array.  It is
assumed that the landmarks have previously been aligned
using Generalized Procrustes Analysis (GPA) [e.g., with
\code{\link{gpagen}}].  The user must also specify a
phylogeny describing the evolutionary relationships among
species (of class phylo).  Note that the specimen labels
for both x and y must match the labels on the tips of the
phylogeny.

From the phylogeny, a phylogenetic transformation matrix is
obtained under a Brownian motion model, and used to
transform the x and y variables. Next, the Gower-centered
distance matrix is obtained from predicted values from the
model (y~x), from which sums-of-squares, F-ratios, and R^2
are estimated for each factor in the model (see Adams,
2014).  Data are then permuted across the tips of the
phylogeny, and estimates of statistical values are obtained
for the permuted data, which are compared to the observed
value to assess significance.
}
\examples{
### Example of D-PGLS for high-dimensional data
data(plethspecies)
Y.gpa<-gpagen(plethspecies$land)    #GPA-alignment

procD.pgls(two.d.array(Y.gpa$coords)~Y.gpa$Csize,plethspecies$phy,iter=49)
}
\author{
Dean Adams
}
\references{
Adams, D.C. 2014. A method for assessing phylogenetic least
squares models for shape and other high-dimensional
multivariate data. Evolution. (In Press).
}
\keyword{analysis}

