\name{trajectory.analysis}
\alias{trajectory.analysis}
\title{Quantify and compare shape change trajectories}
\usage{
  trajectory.analysis(A, traj, lvls, estimate.traj = TRUE,
    iter = 99)
}
\arguments{
  \item{A}{An array (p x k x n) containing GPA-aligned
  coordinates for a set of specimens}

  \item{traj}{A vector of labels assigning each shape to a
  trajectory or group}

  \item{lvls}{A vector of labels assigning each shape to a
  sub-level within each trajectory}

  \item{estimate.traj}{A logical value indicating whether
  trajectories are estimated from original data; described
  below}

  \item{iter}{Number of iterations for significance
  testing}
}
\value{
  Function returns a list with the following components:
  \item{procDist.lm}{Procrustes ANOVA table
  (if"estimate.traj=TRUE")} \item{traj.size}{A matrix of
  pairwise differences in trajectory size} \item{p.size}{A
  matrix of pairwise significance levels for trajectory
  size} \item{traj.orient}{A matrix of pairwise differences
  in trajectory orientation} \item{p.orient}{A matrix of
  pairwise significance levels for trajectory orientation}
  \item{traj.shape}{A matrix of pairwise differences in
  trajectory shape (if applicable)} \item{p.shape}{A matrix
  of pairwise significance levels for trajectory shape}
}
\description{
  Function estimates attributes of shape change
  trajectories for a set of Procrustes-aligned specimens
  and compares them statistically
}
\details{
  The function quantifies phenotypic shape change
  trajectories from a set of specimens, and compares them
  statistically. It is assumed that the landmarks have
  previously been aligned using Generalized Procrustes
  Analysis (GPA) [e.g., with \code{\link{gpagen}}]. A shape
  change trajectory is defined by a sequence of shapes in
  tangent space. These trajectories can be quantified
  various attributes (their size, orientation, and shape),
  and comparisons of these attribute enables the
  statistical comparison of shape change trajectories (see
  Collyer and Adams 2007; Adams and Collyer 2007; Adams and
  Collyer 2009).

  If "estimate.traj=TRUE" the function partitions the set
  of shapes into their trajectory-level groups, and
  estimates the mean shape for each group (e.g.,
  Traj1-Level1; Traj1-Level2, etc.). These are then
  concatenated to form shape change trajectories. A
  Procrustes ANOVA is performed, and differences in
  trajectory attributes (size, orientation, and shape) are
  statistically assessed via residual randomization. This
  approach is useful when summarizing shape changes across
  levels for experimental factors in which the shape of
  multiple individuals have been measured (e.g., comparing
  the shape change from allopatry-to-sympatry across
  species: see Adams 2010).

  If "estimate.traj=FALSE" the trajectories are assembled
  directly from the set of shapes provided in the three
  dimensional array (A). This approach is useful when the
  set of shapes forming each trajectory have been
  quantified directly (e.g., when motion paths are
  compared: see Adams and Cerney 2007). For both methods, a
  plot of specimens and trajectories in tangent space is
  included.
}
\examples{
data(plethodon)
Y.gpa<-gpagen(plethodon$land)    #GPA-alignment

trajectory.analysis(Y.gpa$coords,plethodon$species,plethodon$site,iter=99)
}
\author{
  Dean Adams
}
\references{
  Adams, D. C. 2010. Parallel evolution of character
  displacement driven by competitive selection in
  terrestrial salamanders. BMC Evol. Biol. 10:1-10.

  Adams, D. C., and M. M. Cerney. 2007. Quantifying
  biomechanical motion using Procrustes motion analysis. J.
  Biomech. 40:437-444.

  Adams, D. C., and M. L. Collyer. 2007. The analysis of
  character divergence along environmental gradients and
  other covariates. Evolution 61:510-515.

  Adams, D. C., and M. L. Collyer. 2009. A general
  framework for the analysis of phenotypic trajectories in
  evolutionary studies. Evolution 63:1143-1154.

  Collyer, M. L., and D. C. Adams. 2007. Analysis of
  two-state multivariate phenotypic change in ecological
  studies. Ecology 88:683-692.
}
\keyword{trajectory.analysis}

