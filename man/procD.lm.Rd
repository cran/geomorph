\name{procD.lm}
\alias{procD.lm}
\title{Procrustes ANOVA/regression for shape data}
\usage{
  procD.lm(f1, data = NULL, iter = 999)
}
\arguments{
  \item{f1}{A formula for the linear model (e.g., y~x1+x2)}

  \item{data}{An optional value specifying a data frame
  containing all data (not required)}

  \item{iter}{Number of iterations for significance
  testing}
}
\value{
  Function returns an ANOVA table of statistical results
  for all factors: df (for each factor), SS, MS, Prand.
}
\description{
  Function performs Procrustes ANOVA with permutation
  procedures to assess statistical hypotheses describing
  patterns of shape variation and covariation for a set of
  Procrustes-aligned coordinates
}
\details{
  The function quantifies the relative amount of shape
  variation attributable to one or more factors in a linear
  model and assesses this variation via permutation. Data
  input is specified by a formula (e.g., y~X), where 'y'
  specifies the response variables (shape data), and 'X'
  contains one or more independent variables (discrete or
  continuous). The response matrix 'y' must be in the form
  of a two-dimensional data matrix of dimension (n x
  [p*k]), rather than a 3D array.  It is assumed that the
  landmarks have previously been aligned using Generalized
  Procrustes Analysis (GPA) [e.g., with
  \code{\link{gpagen}}]. The function
  \code{\link{two.d.array}} can be used to obtain a
  two-dimensional data matrix from a 3D array oflandmark
  coordinates. The names specified for the independent (x)
  variables in the formula represent one or more vectors
  containing continuous data or factors. It is assumed that
  the order of the specimens in the shape matrix matches
  the order of values in the independent variables.

  The function performs statistical assessment of the terms
  in the model using Procrustes distances among specimens,
  rather than explained covariance matrices among
  variables. With this approach, the sum-of-squared
  Procrustes distances are used as a measure of SS (see
  Goodall 1991). The observed SS are evaluated through
  permutation, where the rows of the shape matrix are
  randomized relative to the design matrix. Procedurally,
  Procrustes ANOVA is identical to permutational-MANOVA as
  used in other fields (Anderson 2001). For several
  reasons, Procrustes ANOVA is particularly useful for
  shape data. First, covariance matrices from GPA-aligned
  Procrustes coordinates are singular, and thus standard
  approaches such as MANOVA cannot be accomplished unless
  generalized inverses are utilized. This problem is
  accentuated when using sliding semilandmarks.
  Additionally, GM datasets often have more variables than
  specimens (the 'small N large P' problem). In these
  cases, distance-based procedures can still be utilized to
  assess statistical hypotheses, whereas standard linear
  models cannot.
}
\examples{
### MANOVA example for Goodall's F test (multivariate shape vs. factors)
data(plethodon)
Y.gpa<-gpagen(plethodon$land)    #GPA-alignment
y<-two.d.array(Y.gpa$coords)

procD.lm(y~plethodon$species*plethodon$site,iter=99)

### Regression example
data(rats)
rat.gpa<-gpagen(ratland)         #GPA-alignment

procD.lm(two.d.array(rat.gpa$coords)~rat.gpa$Csize,iter=99)
}
\author{
  Dean Adams
}
\references{
  Anderson MJ. 2001. A new method for non-parametric
  multivariate analysis of variance.  Austral Ecology 26:
  32'46.

  Goodall, C. R. 1991. Procrustes methods in the
  statistical analysis of shape. Journal of the Royal
  Statistical Society B 53:285-339.
}
\keyword{procD.lm}

