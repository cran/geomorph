#' Two-block partial least squares analysis for Procrustes shape variables
#'
#' Function performs two-block partial least squares analysis to assess the degree of association between 
#' to blocks of Procrustes shape variables (or other variables)
#'
#' The function quantifies the degree of association between two blocks of shape data as 
#'   defined by Procrustes shape variables using partial least squares (see Rohlf and Corti 2000). If geometric morphometric data are 
#'   used, it is assumed 
#'   that the landmarks have previously been aligned using 
#'   Generalized Procrustes Analysis (GPA) [e.g., with \code{\link{gpagen}}]. If other variables are used, they must be input as a 
#'   2-Dimensional matrix (rows = specimens, columns = variables).  It is also assumed that the separate inputs
#'   have specimens (observations) in the same order.  Additionally, if names for the objects are specified, these must be the same for both datasets.
#'  The observed test value is then compared to a distribution of values obtained by randomly permuting 
#'   the individuals (rows) in one partition relative to those in the other. A significant result is found when the 
#'   observed PLS correlation is large relative to this distribution. In addition, a multivariate effect size describing the strength of the effect is 
#'   estimated from the empirically-generated sampling distribution (see details in Adams and Collyer 2016; 
#'   Adams and Collyer 2019).   
#'   
#'  The generic function, \code{\link{plot}}, produces a two-block.pls plot.  This function calls \code{\link{plot.pls}}, which produces an ordination plot.  
#'  An additional argument allows one to include a vector to label points.  Starting with version 3.1.0, warpgrids are no longer available with \code{\link{plot.pls}}
#'  but after making a plot, the function returns values that can be used with \code{\link{picknplot.shape}} or a combination of 
#' \code{\link{shape.predictor}} and \code{\link{plotRefToTarget}} to visualize shape changes in the plot (via warpgrids).
#'  
#'  \subsection{For more than two blocks}{ 
#' If one wishes to consider 3+ arrays or matrices, there are multiple options.  First, one could perform multiple two.b.pls analyses and use
#' \code{\link{compare.pls}} to ascertain which blocks are more "integrated".  Second, one could use \code{\link{integration.test}} and perform a test that
#' averages the amount of integration (correlations) across multiple pairwise blocks.  Note that performing \code{\link{integration.test}} performed on two matrices or
#' arrays returns the same results as \code{\link{two.b.pls}}.  Thus, \code{\link{integration.test}} is more flexible and thorough.
#' }
#' 
#'  \subsection{Using phylogenies and PGLS}{ 
#' If one wishes to incorporate a phylogeny, \code{\link{phylo.integration}} is the function to use.  This function is exactly the same as \code{\link{integration.test}}
#' but allows PGLS estimation of PLS vectors.  Because \code{\link{integration.test}} can be used on two blocks, \code{\link{phylo.integration}} likewise allows one to
#' perform a phylogenetic two-block PLS analysis.
#' }
#'  
#'  \subsection{Notes for geomorph 3.0}{ 
#' There is a slight change in two.b.pls plots with geomorph 3.0.  Rather than use the shapes of specimens that matched minimum and maximum PLS
#' scores, major-axis regression is used and the extreme fitted values are used to generate deformation grids.  This ensures that shape deformations
#' are exactly along the major axis of shape covariation.  This axis is also shown as a best-fit line in the plot.
#' }
#' 
#' 
#' @param A1 A 3D array (p x k x n) containing Procrustes shape variables for the first block, or a matrix (n x variables)
#' @param A2 A 3D array (p x k x n) containing Procrustes shape variables for the second block, or a matrix (n x variables)
#' @param iter Number of iterations for significance testing
#' @param seed An optional argument for setting the seed for random permutations of the resampling procedure.  
#' If left NULL (the default), the exact same P-values will be found for repeated runs of the analysis (with the same number of iterations).
#' If seed = "random", a random seed will be used, and P-values will vary.  One can also specify an integer for specific seed values,
#' which might be of interest for advanced users.
#' @param print.progress A logical value to indicate whether a progress bar should be printed to the screen.  
#' This is helpful for long-running analyses.
#' @export
#' @keywords analysis
#' @author Dean Adams and Michael Collyer
#' @return Object of class "pls" that returns a list of the following:
#'   \item{r.pls}{The correlation coefficient between scores of projected values on the first
#'   singular vectors of left (x) and right (y) blocks of landmarks (or other variables).  This value can only be negative
#'   if single variables are input, as it reduces to the Pearson correlation coefficient.}
#'   \item{P.value}{The empirically calculated P-value from the resampling procedure.}
#'   \item{Effect.Size}{The multivariate effect size associated with sigma.d.ratio.}
#'   \item{left.pls.vectors}{The singular vectors of the left (x) block}
#'   \item{right.pls.vectors}{The singular vectors of the right (y) block}
#'   \item{random.r}{The correlation coefficients found in each random permutation of the 
#'   resampling procedure.}
#'   \item{XScores}{Values of left (x) block projected onto singular vectors.}
#'   \item{YScores}{Values of right (y) block projected onto singular vectors.}
#'   \item{svd}{The singular value decomposition of the cross-covariances.  See \code{\link{svd}} for further details.}
#'   \item{A1}{Input values for the left block.}
#'   \item{A2}{Input values for the right block.}
#'   \item{A1.matrix}{Left block (matrix) found from A1.}
#'   \item{A2.matrix}{Right block (matrix) found from A2.}
#'   \item{permutations}{The number of random permutations used in the resampling procedure.}
#'   \item{call}{The match call.}
#'   
#' @seealso \code{\link{integration.test}}, \code{\link{modularity.test}}, 
#' \code{\link{phylo.integration}}, and \code{\link{compare.pls}}
#' @references  Rohlf, F.J., and M. Corti. 2000. The use of partial least-squares to study covariation in shape. 
#' Systematic Biology 49: 740-753.
#' @references Adams, D.C. and M.L. Collyer. 2016.  On the comparison of the strength of morphological integration across morphometric 
#' datasets. Evolution. 70:2623-2631.
#' @references Adams, D.C. and M.L. Collyer. 2019. Comparing the strength of modular signal, and evaluating 
#' alternative modular hypotheses, using covariance ratio effect sizes with morphometric data. 
#' Evolution. 73:2352-2367.
#' @examples
#' data(plethShapeFood) 
#' Y.gpa<-gpagen(plethShapeFood$land)    #GPA-alignment    
#'
#' #2B-PLS between head shape and food use data
#' PLS <-two.b.pls(Y.gpa$coords,plethShapeFood$food,iter=999)
#' summary(PLS)
#' plot(PLS)
#'  
#'  ### Visualize shape variation using picknplot.shape Because picknplot  
#'  ### requires user decisions, the following example
#'  ### is not run (but can be with removal of #).
#'  ### For detailed options, see the picknplot help file
#'  # picknplot.shape(plot(PLS))
#'  
#' 

two.b.pls <- function (A1, A2,  iter = 999, seed = NULL, print.progress=TRUE){
    if (any(is.na(A1))) 
      stop("\nData matrix 1 contains missing values. Estimate these first (see 'estimate.missing').",
           call. = FALSE)
    if (any(is.na(A2))) 
      stop("\nData matrix 2 contains missing values. Estimate these first (see 'estimate.missing').",
           call. = FALSE)
  
  x <- try(two.d.array(A1), silent = TRUE)
  if(inherits(x, "try-error")) x <- try(as.matrix(A1), silent = TRUE)
  if(inherits(x, "try-error"))
    stop("\nA is not a suitable data array for analysis. ", call. = FALSE)
  
  y <- try(two.d.array(A2), silent = TRUE)
  if(inherits(y, "try-error")) y <- try(as.matrix(A2), silent = TRUE)
  if(inherits(y, "try-error"))
    stop("\nA is not a suitable data array for analysis. ", call. = FALSE)
  
  x <- as.matrix(x)
  y <- as.matrix(y)
  n <- nrow(x)
  namesX <- rownames(x)
  namesY <- rownames(y)
  cnamesY <- colnames(y)
  if(is.null(namesX) || is.null(namesY))
    cat("Data in either A1 or A2 do not have names.  It is assumed data in both A1 and A2 are ordered the same.\n")
  
  if (is.null(namesX)) namesX <- 1:NROW(x)
  if (is.null(namesY)) {namesY <- namesX
     rownames(y) <- namesY}
  
  if (length(namesX) != length(namesY)) stop("\nData matrices have different numbers of specimens.",
                               call. = FALSE)
  if(length(unique(c(namesX, namesY))) != n) 
    stop("\nMismatched specimen names for A1 and A2.\n", call. = FALSE)
  
  y <- as.matrix(y[match(namesX, namesY), ]); colnames(y) <- cnamesY

  pls.obs <- pls(x, y, RV=FALSE, verbose=TRUE)
  
  if(NCOL(x) > n){
    pcax <- prcomp(x)
    d <- which(zapsmall(pcax$sdev) > 0)
    x <- pcax$x[,d]
  }
  if(NCOL(y) > n){
    pcay <- prcomp(y)
    d <- which(zapsmall(pcay$sdev) > 0)
    y <- pcay$x[,d]
  }
  
  if(!is.null(seed) && seed == "random") seed = sample(1:iter, 1)
  ind <- perm.index(n, iter, seed = seed)
  perms <- length(ind)
  
  if(print.progress){
    cat(paste("\nRandom PLS calculations:", perms, "permutations.\n"))
    pb <- txtProgressBar(min = 0, max = perms+1, initial = 0, style=3)
  }
  
  xc <- as.matrix(center(x))
  yc <- as.matrix(center(y))
  pls.rand <- sapply(1:perms, function(j) {
    step <- j
    if(print.progress) setTxtProgressBar(pb,step)
    s <- ind[[j]]
    xs <- as.matrix(xc[s,])
    quick.pls(xs, yc)
  })

  p.val <- pval(abs(pls.rand))
  Z <- effect.size(pls.rand, center=TRUE) 
  XScores <- pls.obs$XScores
  YScores <- pls.obs$YScores
  out <- list(r.pls = pls.rand[1], P.value = p.val, Z = Z,
              left.pls.vectors = pls.obs$left.vectors,
              right.pls.vectors = pls.obs$right.vectors,
              random.r = pls.rand, 
              XScores = pls.obs$XScores,
              YScores = pls.obs$YScores,
              svd = pls.obs$pls.svd,
              A1 = A1, A2 = A2,
              A1.matrix = x, A2.matrix =y,
              permutations = iter+1, call=match.call(),
              method="PLS")
  class(out) <- "pls"
  out
  }
