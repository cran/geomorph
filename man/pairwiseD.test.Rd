\name{pairwiseD.test}
\alias{pairwiseD.test}
\title{Pairwise Group Comparisons}
\usage{
pairwiseD.test(y, x, iter = 999)
}
\arguments{
  \item{y}{A two-dimensional array of shape data}

  \item{x}{A factor defining groups}

  \item{iter}{Number of iterations for permutation test}
}
\value{
Function returns a list with the following components:
\item{Dist.obs}{A matrix of Euclidean distances among group
means} \item{Prob.Dist}{A matrix of pairwise significance
levels based on permutation}
}
\description{
Function performs pairwise comparisons among groups using
the Euclidean distances among group means.
}
\details{
The function performs pairwise comparisons to identify
shape among groups. The function is designed as a post-hoc
test to Procrustes ANOVA, where the latter has identified
significant shape variation explained by a grouping factor.
The function takes as input the shape data (y), and a
grouping factor (x). It then estimates the Euclidean
distances among group means, which are used as test values.
These are then statistically evaluated through permutation,
where the rows of the shape matrix are randomized relative
to the grouping variable.

The input for the shape data (y) must be in the form of a
two-dimensional data matrix of dimension (n x [p x k]),
rather than a 3D array.  It is assumed that the landmarks
have previously been aligned using Generalized Procrustes
Analysis (GPA) [e.g., with \code{\link{gpagen}}]. The
function \code{\link{two.d.array}} can be used to obtain a
two-dimensional data matrix from a 3D array of landmark
coordinates.
}
\examples{
data(plethodon)
Y.gpa<-gpagen(plethodon$land)    #GPA-alignment
y<-two.d.array(Y.gpa$coords)
### Procrustes ANOVA
procD.lm(y~plethodon$species,iter=99)

### Pairwise comparisons
pairwiseD.test(y,plethodon$species,iter=99)
}
\author{
Dean Adams
}
\keyword{PairwiseD.test}

