% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/phylo.integration.r
\name{phylo.integration}
\alias{phylo.integration}
\title{Quantify phylogenetic morphological integration between two or more sets of variables under Brownian motion}
\usage{
phylo.integration(
  A,
  A2 = NULL,
  phy,
  partition.gp = NULL,
  iter = 999,
  seed = NULL,
  print.progress = TRUE
)
}
\arguments{
\item{A}{A 2D array (n x [p1 x k1]) or 3D array (p1 x k1 x n) containing Procrustes shape variables for the first block}

\item{A2}{An optional 2D array (n x [p2 x k2]) or 3D array (p2 x k2 x n) containing Procrustes shape variables for the second block}

\item{phy}{A phylogenetic tree of class = "phylo" - see \code{\link[ape]{read.tree}} in library ape}

\item{partition.gp}{A list of which landmarks (or variables) belong in which partition: 
(e.g. A, A, A, B, B, B, C, C, C). This is required when only 1 dataset provided.}

\item{iter}{Number of iterations for significance testing}

\item{seed}{An optional argument for setting the seed for random permutations of the resampling procedure.
If left NULL (the default), the exact same P-values will be found for repeated runs of the analysis (with the same number of iterations).
If seed = "random", a random seed will be used, and P-values will vary.  One can also specify an integer for specific seed values,
which might be of interest for advanced users.}

\item{print.progress}{A logical value to indicate whether a progress bar should be printed to the screen.
This is helpful for long-running analyses.}
}
\value{
Objects of class "pls" from integration.test return a list of the following:
 \item{r.pls}{The estimate of morphological integration: PLS.corr. The mean of pairwise
 PLS correlations between partitions is used when there are more than two partitions.}
   \item{r.pls.mat}{The pairwise r.pls, if the number of partitions is greater than 2.}
   \item{P.value}{The empirically calculated P-value from the resampling procedure.}
  \item{Effect.Size}{The multivariate effect size associated with sigma.d.ratio.}
   \item{left.pls.vectors}{The singular vectors of the left (x) block (for 2 modules only).}
   \item{right.pls.vectors}{The singular vectors of the right (y) block (for 2 modules only).}
   \item{random.r}{The correlation coefficients found in each random permutation of the 
  resampling procedure.}
   \item{XScores}{Values of left (x) block projected onto singular vectors 
  (for 2 modules only).}
   \item{YScores}{Values of right (y) block projected onto singular vectors
  (for 2 modules only).}
   \item{svd}{The singular value decomposition of the cross-covariances (for 2 modules only).}
   \item{A1}{Input values for the left block (for 2 modules only).}
   \item{A2}{Input values for the right block (for 2 modules only).}
   \item{A1.matrix}{Left block (matrix) found from A1 (for 2 modules only).}
   \item{A2.matrix}{Right block (matrix) found from A2 (for 2 modules only).}
   \item{Pcov}{The phylogenetic transformation matrix, needed for certain other analyses.}
   \item{permutations}{The number of random permutations used in the resampling procedure.}
   \item{call}{The match call.}
}
\description{
Function quantifies the degree of phylogenetic morphological covariation between two or more sets of
Procrustes shape variables using partial least squares.
}
\details{
The function quantifies the degree of phylogenetic morphological integration between two or more sets of Procrustes shape variables. 
The approach is based on a Brownian motion model of evolution. It is 
  assumed that the landmarks have previously been aligned using 
  Generalized Procrustes Analysis (GPA) [e.g., with \code{\link{gpagen}}].   

The function estimates the degree of  morphological covariation between two or sets of variables 
while accounting for phylogeny using partial least squares (Adams and Felice 2014), and under a Brownian
motion model of evolution. If more than two partitions are defined, the average pairwise PLS correlation is 
utilized as the test statistic. The observed value is statistically assessed using permutation, where data for 
one partition are permuted relative to the other partitions. In addition, a multivariate effect size 
describing the strength of the effect is estimated from the empirically-generated sampling distribution 
(see details in Adams and Collyer 2019). Note that this permutation is performed on phylogenetically-
transformed data, so that the probability of phylogenetic association of A vs. B is similar to that of B vs. A: 
i.e., prob(A,B|phy)~prob(B,A|phy); thus, shuffling the correct exchangeable units under the null 
hypothesis of no integration (Adams and Collyer 2018). 

 Input for the analysis can take one of two forms. First, one can input a single dataset (as a matrix or 3D array, along with 
 a vector describing which variables correspond to which partitions (for the case of a 3D array, which landmarks belong to which 
 partitions is specified). Alternatively, when evaluating the integration between two structures or partitions, two datasets may be provided.

 The generic functions, \code{\link{print}}, \code{\link{summary}}, and \code{\link{plot}} all work with \code{\link{phylo.integration}}.  
 The generic function, \code{\link{plot}}, produces a two-block.pls plot.  This function calls \code{\link{plot.pls}}, which produces an ordination plot.  
 An additional argument allows one to include a vector to label points.  Starting with version 3.1.0, warpgrids are no longer available with \code{\link{plot.pls}}
 but after making a plot, the function returns values that can be used with \code{\link{picknplot.shape}} or a combination of 
\code{\link{shape.predictor}} and \code{\link{plotRefToTarget}} to visualize shape changes in the plot (via warpgrids).
 
\subsection{Similarity to \code{\link{two.b.pls}} and \code{\link{compare.pls}} }{ 
Note that \code{phylo.integration} performed on two matrices or arrays returns the same results as a phylogenetic variation of
 \code{\link{two.b.pls}}.  It might be of interest with 3+ modules to perform separate phylogenetic integration tests
between all pairwise comparisons of modules.  This can be done, test by test, and the levels of integration can be compared with
\code{\link{compare.pls}}.  Such results are different than using the average amount of integration when more than two modules 
are input, as found with \code{phylo.integration}.
}
 
 \subsection{Notes for geomorph 3.0.4 and subsequent versions}{ 
 Compared to previous versions of geomorph, users might notice differences in effect sizes.  Previous versions used z-scores calculated with 
 expected values of statistics from null hypotheses (sensu Collyer et al. 2015); however Adams and Collyer (2016) showed that expected values 
 for some statistics can vary with sample size and variable number, and recommended finding the expected value, empirically, as the mean from the set 
 of random outcomes.  Geomorph 3.0.4 and subsequent versions now center z-scores on their empirically estimated expected values and where appropriate, 
 log-transform values to assure statistics are normally distributed.  This can result in negative effect sizes, when statistics are smaller than 
 expected compared to the average random outcome.  For ANOVA-based functions, the option to choose among different statistics to measure effect size 
 is now a function argument.
}
}
\examples{
\dontrun{
data(plethspecies) 
Y.gpa <- gpagen(plethspecies$land)    #GPA-alignment
land.gps <- c("A","A","A","A","A","B","B","B","B","B","B") 

IT <- phylo.integration(Y.gpa$coords, partition.gp = land.gps,
  phy = plethspecies$phy)
summary(IT) # Test summary
P <- plot(IT) # PLS plot

 # Block 1 
 minx <- min(P$plot_args$x)
 maxx <- max(P$plot_args$x)
 preds <- shape.predictor(P$A1, 
 x = P$plot.args$x,
 min = minx, max = maxx)
 plotRefToTarget(mshape(P$A1), preds$min)
 plotRefToTarget(mshape(P$A1), preds$max)

 # Block 2 
 miny <- min(P$plot_args$y)
 maxy <- max(P$plot_args$y)
 preds <- shape.predictor(P$A2, 
 x = P$plot.args$y,
 min = miny, max = maxy)
 plotRefToTarget(mshape(P$A2), preds$min)
 plotRefToTarget(mshape(P$A2), preds$max)
 
 ### Visualize shape variation using picknplot.shape Because picknplot  
 ### requires user decisions, the following example
 ### is not run.
 ### For detailed options, see the picknplot help file
 # picknplot.shape(P)
}
}
\references{
Adams, D.C. and R. Felice. 2014. Assessing phylogenetic morphological 
integration and trait covariation in morphometric data using evolutionary covariance 
matrices. PLOS ONE. 9(4):e94335.

Collyer, M.L., D.J. Sekora, and D.C. Adams. 2015. A method for analysis of phenotypic change for phenotypes described 
by high-dimensional data. Heredity. 115:357-365.

Adams, D.C. and M.L. Collyer. 2016.  On the comparison of the strength of morphological integration across morphometric 
datasets. Evolution. 70:2623-2631.

Adams, D.C. and M.L. Collyer. 2018. Multivariate comparative methods: evaluations, comparisons, and
recommendations. Systematic Biology. 67:14-31.

Adams, D.C. and M.L. Collyer. 2019. Comparing the strength of modular signal, and evaluating 
alternative modular hypotheses, using covariance ratio effect sizes with morphometric data. 
Evolution. 73:2352-2367.
}
\seealso{
\code{\link{integration.test}}, \code{\link{modularity.test}}, and 
\code{\link{two.b.pls}}
}
\author{
Dean Adams
}
\keyword{analysis}
