% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/compare.evol.rates.r
\name{compare.evol.rates}
\alias{compare.evol.rates}
\title{Comparing net rates of shape evolution on phylogenies}
\usage{
compare.evol.rates(
  A,
  phy,
  gp,
  iter = 999,
  seed = NULL,
  method = c("permutation", "simulation"),
  print.progress = TRUE
)
}
\arguments{
\item{A}{A 3D array (p x k x n) containing GPA-aligned coordinates for all specimens, or a matrix (n x variables)}

\item{phy}{A phylogenetic tree of class = "phylo" - see \code{\link[ape]{read.tree}} in library ape}

\item{gp}{A factor array designating group membership for individuals}

\item{iter}{Number of iterations for significance testing}

\item{seed}{An optional argument for setting the seed for random permutations of the resampling procedure.
If left NULL (the default), the exact same P-values will be found for repeated runs of the analysis (with the same number of iterations).
If seed = "random", a random seed will be used, and P-values will vary.  One can also specify an integer for specific seed values,
which might be of interest for advanced users.}

\item{method}{One of "simulation" or "permutation", to choose which approach should be used to assess significance.}

\item{print.progress}{A logical value to indicate whether a progress bar should be printed to the screen.
This is helpful for long-running analyses.}
}
\value{
An object of class "evolrate" returns a list with the following components: 
  \item{sigma.d.ratio}{The ratio of maximum to minimum net evolutionary rates.}
  \item{P.value}{The significance level of the observed ratio.}
  \item{Effect.Size}{The multivariate effect size associated with sigma.d.ratio.}
  \item{sigma.d.gp}{The phylogenetic net evolutionary rate for each group of species on the phylogeny.}
  \item{random.sigma}{The sigma values found in random permutations of the resampling procedure.}
  \item{permutations}{The number of random permutations used.}
}
\description{
Function calculates net rates of shape evolution for two or more groups of species on a phylogeny from a set of Procrustes-aligned specimens
}
\details{
The function compares net rates of morphological evolution for two or more groups of species on a phylogeny, under a 
 Brownian motion model of evolution. It is assumed that the landmarks have previously been aligned 
 using Generalized Procrustes Analysis (GPA) [e.g., with \code{\link{gpagen}}]. The approach is based on the outer-product matrix of  
 between species differences in morphospace after phylogenetic transformation (Adams 2014). From the data the net rate of shape evolution
 for each group in the multi-dimensional space is calculated, and a ratio of rates is obtained. If three or more groups of species are used, the ratio of 
 the maximum to minimum rate is used as a test statistic (see Adams 2014). The function can be used with univariate data (i.e. 
 centroid size) if imported as matrix with rownames giving the taxa names.
 
 The generic functions, \code{\link{print}}, \code{\link{summary}}, and \code{\link{plot}} all work with \code{\link{compare.evol.rates}}.
 The generic function, \code{\link{plot}}, produces a histogram of random rate-ratios associated with
 the resampling procedure.

\subsection{Notes for geomorph 3.0.4 and subsequent versions}{ 
Significance testing is now accomplished in one of two ways. First, phylogenetic simulation may be used, in which tips data are 
obtained under Brownian motion using a common evolutionary rate pattern for all species on the phylogeny. Specifically, the 
common evolutionary rate matrix for all species is used, with the multi-dimensional rate used along the diagonal elements (see 
Denton and Adams 2015). This procedure is more general than the original simulation procedure, and retains the desirable 
statistical properties of earlier methods, and under a wider array of data types.  Second, significance may be accomplished via 
permutation, where data values at the tips are permuted relative to the (see Adams and Collyer 2018). This procedure is shown to 
retain all appropriate statistical properties, including rotation-invariance of significance levels (see results of Adams and Collyer 2018).
In addition, a multivariate effect size describing the strength of the effect is estimated from the 
empirically-generated sampling distribution (see details in Adams and Collyer 2019). Values from these 
distributions are log-transformed prior to effect size estimation, to assure normally distributed data. 
}
}
\examples{
\dontrun{

data(plethspecies) 
Y.gpa <- gpagen(plethspecies$land)    #GPA-alignment    
 gp.end <- factor(c(0,0,1,0,0,1,1,0,0))  #endangered species vs. rest
 names(gp.end) <- plethspecies$phy$tip

ER<-compare.evol.rates(A = Y.gpa$coords, phy = plethspecies$phy,
  method = "simulation", gp = gp.end)
summary(ER)
plot(ER)
}
}
\references{
Adams, D.C. 2014. Quantifying and comparing phylogenetic evolutionary rates for 
 shape and other high-dimensional phenotypic data. Syst. Biol. 63:166-177.

Denton, J.S.S., and D.C. Adams. 2015. A new phylogenetic test for comparing 
multiple high-dimensional evolutionary rates suggests interplay of evolutionary rates and 
modularity in lanternfishes (Myctophiformes; Myctophidae). Evolution. 69:2425-2440.

Adams, D.C. and M.L. Collyer. 2018. Multivariate comparative methods: evaluations, comparisons, and
recommendations. Systematic Biology. 67:14-31.

Adams, D.C. and M.L. Collyer. 2019. Comparing the strength of modular signal, and evaluating 
alternative modular hypotheses, using covariance ratio effect sizes with morphometric data. 
Evolution. 73:2352-2367.
}
\author{
Dean Adams & Emma Sherratt
}
\keyword{analysis}
