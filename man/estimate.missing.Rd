\name{estimate.missing}
\alias{estimate.missing}
\title{Estimate locations of missing landmarks using the thin-plate spline}
\usage{
  estimate.missing(ref, target)
}
\arguments{
  \item{ref}{Landmark coordinates of a reference specimen}

  \item{target}{Landmark coordinates of a target specimen
  (contains -999 for missing landmarks)}
}
\value{
  Function returns a n * p matrix of coordinates for the
  target specimen that includes the original landmarks plus
  the estimated coordinates for the missing landmarks.
}
\description{
  A function for estimating the locations of missing
  landmarks
}
\details{
  The function estimates the locations of missing landmarks
  on a target specimen based on the locations of the
  corresponding landmarks on a complete, reference
  specimen. Missing landmarks in a target specimen are
  designated by '-999' in place of the x,y,z coordinates.
  First, a complete reference specimen is aligned to a
  target specimen (missing one or more landmarks), using
  the set of landmarks common to both specimens. Next, the
  thin-plate spline interpolating function is used to
  estimate the locations of the missing landmarks in the
  target specimen.
}
\examples{
data(plethodon)
ref<-plethodon$land[,,1]
target<-plethodon$land[,,11]
target[2,]<-target[6,]<- -999    #create some missing landmarks
estimate.missing(ref,target)
}
\author{
  Dean Adams
}
\references{
  Bookstein, F. L., K. Schafer, H. Prossinger, H. Seidler,
  M. Fieder, G. Stringer, G. W. Weber, J.-L. Arsuaga, D. E.
  Slice, F. J. Rohlf, W. Recheis, A. J. Mariam, and L. F.
  Marcus. 1999. Comparing frontal cranial profiles in
  archaic and modern Homo by morphometric analysis. Anat.
  Rec. (New Anat.) 257:217-224.

  Gunz, P., P. Mitteroecker, S. Neubauer, G. W. Weber, and
  F. L. Bookstein. 2009. Principles for the virtual
  reconstruction of hominin crania. J. Hum. Evol. 57:48-62.
}

