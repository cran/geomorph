\name{morphol.integr}
\alias{morphol.integr}
\title{Quantify morphological integration between two modules}
\usage{
  morphol.integr(A, WithinConfig = FALSE, A2 = NULL,
    landgroups = NULL, method = c("PLS", "RV"), iter = 999)
}
\arguments{
  \item{A}{An array (p x k x n) containing landmark
  coordinates for the first module}

  \item{WithinConfig}{A logical value indicating whether
  morphological integration is to be assessed within or
  across structures}

  \item{A2}{An optional array (p x k x n) containing
  landmark coordinates for the second module
  (WithinConfig=FALSE)}

  \item{landgroups}{A list of which landmarks belong in
  which partition (WithinConfig=TRUE)}

  \item{method}{Method to estimate morphological
  integration; see below for details}

  \item{iter}{Number of iterations for significance
  testing}
}
\value{
  Function returns a list with the following components:
  \item{value}{The estimate of morphological integration:
  PLS.corr or RV} \item{pvalue}{The significance level of
  the observed signal} \item{Xscores}{PLS scores for the
  first block of landmarks (PLS method only)}
  \item{Yscores}{PLS scores for the second block of
  landmarks (PLS method only)}
}
\description{
  Function quantifies the degree of morphological
  integration between two modules of Procrustes-aligned
  coordinates
}
\details{
  The function quantifies the degree of morphological
  integration between two modules of shape data as defined
  by landmark coordinates. It is assumed that the landmarks
  have previously been aligned using Generalized Procrustes
  Analysis (GPA) [e.g., with \code{\link{gpagen}}]. The
  function may be used to assess the degree of
  morphological integration between two separate structures
  (WithinConfig=FALSE) or between two modules defined
  within the same landmark configuration
  (WithinConfig=TRUE).

  For both options, two analytical approaches are currently
  implemented to assess the degree of morphological
  integration. If "method=PLS" (the default) the function
  estimates the degree of morphological integration using
  two-block partial least squares, or PLS. When used with
  landmark data, this analysis is referred to as singular
  warps analysis (Bookstein et al. 2003). When
  "method=PLS", the scores along the X & Y PLS axes are
  also returned.  Alternatively, if "method=RV" the
  function estimates the degree of morphological
  integration using the RV coefficient (Klingenberg 2009).
  Significance testing for both approaches is found by
  permuting the objects in one data matrix relative to
  those in the other. A histogram of coefficients obtained
  via resampling is presented, with the observed value
  designated by an arrow in the plot.

  If the degree of morphological integration between more
  than two sets of landmarks is of interest, one may use
  the average RV coefficent as implemented in the function
  \code{\link{compare.modular.partitions}}.
}
\examples{
data(plethodon)
Y.gpa<-gpagen(plethodon$land)    #GPA-alignment
land.gps<-c("A","A","A","A","A","B","B","B","B","B","B","B")

#Morphological integration using PLS
morphol.integr(Y.gpa$coords,landgroups=land.gps,WithinConfig=TRUE,method="PLS",iter=99)

#Morphological integration using RV
morphol.integr(Y.gpa$coords,landgroups=land.gps,WithinConfig=TRUE,method="RV",iter=99)
}
\author{
  Dean Adams
}
\references{
  Bookstein, F. L., P. Gunz, P. Mitteroecker, H.
  Prossinger, K. Schaefer, and H. Seidler.  2003. Cranial
  integration in Homo: singular warps analysis of the
  midsagittal plane in ontogeny and evolution. J. Hum.
  Evol. 44:167-187.

  Klingenberg, C. P. 2009. Morphometric integration and
  modularity in configurations of landmarks: tools for
  evaluating a priori hypotheses. Evol. Develop.
  11:405-421.
}
\keyword{morphol.integr}

