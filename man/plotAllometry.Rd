\name{plotAllometry}
\alias{plotAllometry}
\title{Plot allometric patterns in landmark data}
\usage{
  plotAllometry(A, sz, groups = NULL,
    method = c("CAC", "RegScore", "PredLine"),
    warpgrids = TRUE, label = FALSE)
}
\arguments{
  \item{A}{An array (p x k x n) containing landmark
  coordinates for a set of specimens}

  \item{sz}{A vector of centroid size measures for all
  specimens}

  \item{groups}{An optional vector containing group labels
  for each specimen if available}

  \item{method}{Method for estimating allometric shape
  components; see below for details}

  \item{warpgrids}{A logical value indicating whether
  deformation grids for small and large shapes should be
  displayed}

  \item{label}{A logical value indicating whether labels
  for each specimen should be displayed}
}
\description{
  Function plots allometry curves for a set of specimens
}
\details{
  The function generates a plot that describes the
  multivariate relationship between size and shape derived
  from landmark data (ie. allometry). It is assumed that
  the landmarks have previously been aligned using
  Generalized Procrustes Analysis (GPA) [e.g., with
  \code{\link{gpagen}}]. The abscissa of the plot is
  log(centroid size) while the ordinate represents shape.
  Three complementary approaches can be implemented to
  visualize allometry. If "method=CAC" (the default) the
  function calculates the common allometric component of
  the shape data, which is an estimate of the average
  allometric trend within groups (Mitteroecker et al.
  2004). If "method=RegScore" the function calculates shape
  scores from the regression of shape on size, and plots
  these versus size (Drake and Klingenberg 2008).  For a
  single group, these shape scores are mathematically
  identical to the CAC (Adams et al. 2012).  If
  "method=PredLine" the function calculates predicted
  values from a regression of shape on size, and plots the
  first principal component of the predicted values versus
  size as a stylized graphic of the allometric trend (Adams
  and Nistri 2010). For all methods, both centroid size and
  allometry scores are returned. Optionally, deformation
  grids can be requested, which display the shape of the
  smallest and largest specimens relative to the average
  specimen (using 'warpgrid=T' or 'warpgrid=F'). Finally,
  if groups are provided, the above approaches are
  implemented while accounting for within-group patterns of
  covariation (see references for explanation).
}
\examples{
data(ratland)
Y.gpa<-gpagen(ratland)    #GPA-alignment

#Using CAC for plot
plotAllometry(Y.gpa$coords,Y.gpa$Csize,method="CAC")

#Using Regression Scores for plot
plotAllometry(Y.gpa$coords,Y.gpa$Csize,method="RegScore")

#Using predicted allometry curve for plot
plotAllometry(Y.gpa$coords,Y.gpa$Csize,method="PredLine")
}
\author{
  Dean Adams
}
\references{
  Adams, D.C., F.J. Rohlf, and D.E. Slice. 2013. A field
  comes of age: geometric morphometrics in the 21st
  century. Hystrix. 24:7-14.

  Adams, D. C., and A. Nistri. 2010. Ontogenetic
  convergence and evolution of foot morphology in European
  cave salamanders (Family: Plethodontidae). BMC Evol.
  Biol. 10:1-10.

  Drake, A. G., and C. P. Klingenberg. 2008. The pace of
  morphological change: Historical transformation of skull
  shape in St Bernard dogs. Proceedings of the Royal
  Society B, Biological Sciences 275:71'76.

  Mitteroecker, P., P. Gunz, M. Bernhard, K. Schaefer, and
  F. L. Bookstein. 2004.  Comparison of cranial ontogenetic
  trajectories among great apes and humans. J. Hum. Evol.
  46:679-698.
}
\keyword{plotAllometry}

