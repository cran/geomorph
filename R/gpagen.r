#' Generalized Procrustes analysis of points, curves, and surfaces
#'
#' A general function to perform Procrustes analysis of two- or three-dimensional landmark data that 
#'  can include both fixed landmarks and sliding semilandmarks
#'
#' The function performs a Generalized Procrustes Analysis (GPA) on two-dimensional or three-dimensional
#'  landmark coordinates. The analysis can be performed on fixed landmark points, semilandmarks on 
#'  curves, semilandmarks on surfaces, or any combination. If data are provided in the form of a 3D array, all
#'  landmarks and semilandmarks are contained in this object. If this is the only component provided, the function
#'  will treat all points as if they were fixed landmarks. To designate some points as semilandmarks, one uses 
#'  the "curves =" or "surfaces =" options (or both). To include semilandmarks on curves, a matrix defining 
#'  which landmarks are to be treated as semilandmarks is provided using the "curves =" option. This matrix contains
#'  three columns that specify the semilandmarks and two neighboring landmarks which are used to specify the tangent 
#'  direction for sliding. The matrix may be generated using the function \code{\link{define.sliders}}). Likewise, 
#'  to include semilandmarks 
#'  on surfaces, one must specify a vector listing which landmarks are to be treated as surface semilandmarks 
#'  using the "surfaces=" option. The "ProcD = FALSE" option (the default) will slide the semilandmarks 
#'  based on minimizing bending energy, while "ProcD = TRUE" will slide the semilandmarks along their tangent 
#'  directions using the Procrustes distance criterion. The Procrustes aligned specimens may be projected into tangent
#'  space using the "Proj = TRUE" option. 
#'  The function also outputs a matrix of pairwise Procrustes Distances, which correspond to Euclidean distances between specimens in tangent space if "Proj=TRUE", or to the geodesic distances in shape space if "Proj=FALSE".   
#'  NOTE: Large datasets may exceed the memory limitations of R. 
#'
#'  Generalized Procrustes Analysis (GPA: Gower 1975, Rohlf and Slice 1990) is the primary means by which 
#'   shape variables are obtained from landmark data (for a general overview of geometric morphometrics see 
#'   Bookstein 1991, Rohlf and Marcus 1993, Adams et al. 2004, Zelditch et al. 2012, Mitteroecker and 
#'   Gunz 2009, Adams et al. 2013). GPA translates all specimens to the origin, scales them to unit-centroid 
#'   size, and optimally rotates them (using a least-squares criterion) until the coordinates of corresponding
#'   points align as closely as possible. The resulting aligned Procrustes coordinates represent the shape 
#'   of each specimen, and are found in a curved space related to Kendall's shape space (Kendall 1984). 
#'   Typically, these are projected into a linear tangent space yielding Kendall's tangent space coordinates 
#'   (i.e., Procrustes shape variables), which are used for subsequent multivariate analyses (Dryden and Mardia 1993, Rohlf 1999). 
#'   Additionally, any semilandmarks on curves and surfaces are slid along their tangent directions or tangent planes during the 
#'   superimposition (see Bookstein 1997; Gunz et al. 2005). Presently, two implementations are possible: 
#'   1) the locations of semilandmarks can be optimized by minimizing the bending energy between the 
#'   reference and target specimen (Bookstein 1997), or by minimizing the Procrustes distance between the two 
#'   (Rohlf 2010).  Note that specimens are NOT automatically reflected to improve the GPA-alignment.
#'   
#'   The generic functions, \code{\link{print}}, \code{\link{summary}}, and \code{\link{plot}} all work with \code{\link{gpagen}}.
#'   The generic function, \code{\link{plot}}, calls \code{\link{plotAllSpecimens}}.
#'
#'  \subsection{Notes for geomorph 4.0}{ 
#'  Starting with geomorph version 4.0, it is possible to use approximated thin-plate spline (TPS) mapping
#'  to estimate bending energy (when bending energy is used for sliding semi-landmarks).  Approximated TPS
#'  has some computational benefit for large data sets (more GPA iterations are possible in a shorter time), but
#'  one should make sure first that results are reasonable.  Approximated TPS uses a subset of points (semilandmarks) 
#'  to estimate bending energy locally, where points are free to slide.  Some subsets might be misrepresentative of
#'  the global (full configuration) bending energy, and strange results are possible.  A comparison between TPS
#'  and approximated TPS outcomes could be tried with a few specimens to verify consistency of results before using
#'  approximated TPS on a large data set.
#'  
#'  Choosing to use approximated TPS with only a few sliders is dangerous, as the subset of points might be too small
#'  to be reliable.  Approximated TPS will work best for data sets with a sufficient number of surface landmarks.  If a
#'  data set is nearly 100 \% semilandmarks, TPS and approximated TPS should converge, but no time savings is likely.
#' }
#'  \subsection{Notes for geomorph 3.0}{ 
#' Compared to older versions of geomorph, users might notice subtle differences in Procrustes shape variables when using
#' semilandmarks (curves or surfaces).  This difference is a result of using recursive updates of the 
#' consensus configuration with the sliding algorithms (minimized bending energy or Procrustes distances).  
#' (Previous versions used a single consensus through the sliding algorithms.)  Shape differences using the recursive 
#' updates of the consensus configuration should be highly correlated with shape differences using a single consensus 
#' during the sliding algorithm, but rotational "flutter" can be expected.  This should have no qualitative effect on 
#' inferential analyses using Procrustes residuals. 
#' }

#' @param A Either an object of class geomorphShapes or a 3D array (p x k x n) containing landmark coordinates 
#' for a set of specimens.  If A is a geomorphShapes object, the curves argument is not needed.
#' @param PrinAxes A logical value indicating whether or not to align the shape data by principal axes 
#' @param curves An optional matrix defining which landmarks should be treated as semilandmarks on boundary 
#'   curves, and which landmarks specify the tangent directions for their sliding.  This matrix is generated 
#'   with \code{\link{readland.shapes}} following digitizing of curves in StereoMorph, or may be generated
#'   using the function \code{\link{define.sliders}}.
#' @param surfaces An optional vector defining which landmarks should be treated as semilandmarks on surfaces
#' @param rot.pts An optional vector defining which landmarks are used for rotation.  If NULL or "all", all fixed landmarks
#' and semilandmarks will be used.  If "fixed", only fixed landmarks will be used.  If "sliders", only sliding
#' semilandmarks will be used. Any other vector of values will
#' direct which landmarks or semilandmarks will be used for rotation.  Illogical vectors (like choosing points that 
#' do not exist or only choosing one point) will default to using all landmarks and semilandmarks.
#' @param ProcD A logical value indicating whether or not Procrustes distance should be used as the criterion
#'   for optimizing the positions of semilandmarks (if not, bending energy is used)
#' @param approxBE A logical value for whether bending energy should be estimated via approximate 
#' thin-plate spline (TPS) mapping for sliding semilandmarks.  Approximate TPS mapping is much faster and allows 
#' for more iterations, which might return more reliable results than few iterations with full TPS.  
#' If using full TPS, one should probably change max.iter to be few for large data sets.
#' @param sen A numeric value between 0.1 and 1 to adjust the sensitivity of true bending 
#' energy to use in an approximated thin plate mapping of bending energy. This value is the proportion of landmarks
#' (excluding semilandmarks) to seek for estimating bending energy.  Sliding semilandmarks are always included in
#' the final set of landmarks used to make estimates, so the actual sensitivity is higher than the chosen value.
#' @param max.iter The maximum number of GPA iterations to perform before superimposition is halted.  The final
#' number of iterations will be larger than max.iter, if curves or surface semilandmarks are involved, as max.iter
#' will pertain only to iterations with sliding of landmarks.
#' @param tol A numeric value that can be adjusted for evaluation of the convergence criterion.  If the criterion
#' is lower than the tolerance, iterations will be stopped.
#' @param Proj A logical value indicating whether or not the Procrustes aligned specimens should be projected 
#'   into tangent space 
#' @param verbose A logical value for whether to include statistics that take a long time to 
#' calculate with large data sets, including a Procrustes distance matrix among specimens, 
#' a variance-covariance matrix among 
#' Procrustes coordinates (shape variables), and variances of landmark points.  Formatting a data frame
#' for coordinates and centroid size is also embedded within this option.
#' This argument should be FALSE unless 
#' explicitly needed.  For large data sets, it will slow down the analysis, extensively.  
#' @param print.progress A logical value to indicate whether a progress bar should be printed to the screen.  
#' @param Parallel For sliding semi-landmarks only, either a logical value to indicate whether 
#' parallel processing should be used or a numeric value to indicate the number of cores to use in 
#' parallel processing via the \code{parallel} library. 
#' If TRUE, this argument invokes forking of all processor cores, except one.  If
#' FALSE, only one core is used. A numeric value directs the number of cores to use,
#' but one core will always be spared.  Parallel processing is probably only valuable with
#' large data sets.
#' @keywords analysis
#' @export
#' @author Dean Adams and Michael Collyer
#' @return An object of class gpagen returns a list with the following components:
#'  \item{coords}{A (p x k x n) array of Procrustes shape variables, where p is the number of landmark 
#'     points, k is the number of landmark dimensions (2 or 3), and n is the number of specimens. The third 
#'     dimension of this array contains names for each specimen if specified in the original input array.}
#'  \item{Csize}{A vector of centroid sizes for each specimen, containing the names for each specimen if 
#'     specified in the original input array.}
#'  \item{iter}{The number of GPA iterations until convergence was found (or GPA halted).}
#'  \item{points.VCV}{Variance-covariance matrix among Procrustes shape variables.}
#'  \item{points.var}{Variances of landmark points.}
#'  \item{consensus}{The consensus (mean) configuration.}
#'  \item{procD}{Procrustes distance matrix for all specimens (see details). Note that for large data
#'  sets, R might return a memory allocation error, in which case the error will be suppressed and this component will be NULL.
#'  For such cases, users can augment memory allocation and create distances with the dist function, independent from gpagen,
#'  using the coords or data output.}
#'  \item{p}{Number of landmarks.}
#'  \item{k}{Number of landmark dimensions.}
#'  \item{nsliders}{Number of semilandmarks along curves.}
#'  \item{nsurf}{Number of semilandmarks as surface points.}
#'  \item{data}{Data frame with an n x (pk) matrix of Procrustes shape variables and centroid size.}
#'  \item{Q}{Final convergence criterion value.}
#'  \item{slide.method}{Method used to slide semilandmarks.}
#'  \item{call}{The match call.}
#' @references  Adams, D. C., F. J. Rohlf, and D. E. Slice. 2004. Geometric morphometrics: ten years of 
#'    progress following the 'revolution'. It. J. Zool. 71:5-16.
#' @references Adams, D. C., F. J. Rohlf, and D. E. Slice. 2013. A field comes of age: Geometric 
#'   morphometrics in the 21st century. Hystrix.24:7-14.
#' @references Bookstein, F. L. 1991. Morphometric tools for landmark data: Geometry and Biology. 
#'  Cambridge Univ. Press, New York.
#' @references Bookstein, F. L. 1997. Landmark methods for forms without landmarks: morphometrics of 
#'   group differences in outline shape.  1:225-243.
#' @references Dryden, I. L., and K. V. Mardia. 1993. Multivariate shape analysis. Sankhya 55:460-480.
#' @references Gower, J. C. 1975. Generalized Procrustes analysis. Psychometrika 40:33-51.
#' @references Gunz, P., P. Mitteroecker, and F. L. Bookstein. 2005. semilandmarks in three dimensions. 
#'   Pp. 73-98 in D. E. Slice, ed. Modern morphometrics in physical anthropology. Klewer Academic/Plenum, New York.
#' @references Kendall, D. G. 1984. Shape-manifolds, Procrustean metrics and complex projective spaces. 
#'   Bulletin of the London Mathematical Society 16:81-121.
#' @references Mitteroecker, P., and P. Gunz. 2009. Advances in geometric morphometrics. Evol. Biol. 36:235-247.
#' @references Rohlf, F. J., and D. E. Slice. 1990. Extensions of the Procrustes method for the optimal 
#'   superimposition of landmarks. Syst. Zool. 39:40-59.
#' @references Rohlf, F. J., and L. F. Marcus. 1993. A revolution in morphometrics. Trends Ecol. Evol. 8:129-132.
#' @references Rohlf, F. J. 1999. Shape statistics: Procrustes superimpositions and tangent spaces. 
#'   Journal of Classification 16:197-223.
#' @references Rohlf, F. J. 2010. tpsRelw: Relative warps analysis. Version 1.49. Department of Ecology and 
#'   Evolution, State University of New York at Stony Brook, Stony Brook, NY.
#' @references Zelditch, M. L., D. L. Swiderski, H. D. Sheets, and W. L. Fink. 2012. Geometric morphometrics 
#'   for biologists: a primer. 2nd edition. Elsevier/Academic Press, Amsterdam.
#' @examples
#' \dontrun{
#' # Example 1: fixed points only
#' data(plethodon) 
#' Y.gpa <- gpagen(plethodon$land, PrinAxes = FALSE)
#' summary(Y.gpa)
#' plot(Y.gpa)
#' 
#' # Example 2: points and semilandmarks on curves
#' data(hummingbirds)
#' 
#' ###Slider matrix
#' hummingbirds$curvepts
#'
#' # Using bending energy for sliding
#' Y.gpa <- gpagen(hummingbirds$land, curves = hummingbirds$curvepts,
#' ProcD = FALSE)   
#' summary(Y.gpa)
#' plot(Y.gpa)
#' 
#' 
#' # Using Procrustes Distance for sliding
#' Y.gpa <- gpagen(hummingbirds$land, curves = hummingbirds$curvepts,
#' ProcD = TRUE)   
#' summary(Y.gpa)
#' plot(Y.gpa)
#' 
#' 
#' # Example 3: points, curves and surfaces
#' data(scallops)
#' 
#' # Using Procrustes Distance for sliding
#' Y.gpa <- gpagen(A = scallops$coorddata, curves = scallops$curvslide, 
#' surfaces = scallops$surfslide)
#' # NOTE can summarize as: summary(Y.gpa)
#' # NOTE can plot as: plot(Y.gpa) 
#' }
gpagen = function(A, curves=NULL, surfaces=NULL, 
                  rot.pts = NULL,
                  PrinAxes = TRUE, 
                  max.iter = NULL, tol = 1e-4, 
                  ProcD=FALSE, approxBE = FALSE, 
                  sen = 0.5,
                  Proj = TRUE, verbose = FALSE,
                  print.progress = TRUE, Parallel = FALSE){
  
  if(Parallel != FALSE) {
    if(is.numeric(Parallel) && Parallel <= 1) Parallel <- FALSE
  }
  
  if(inherits(A, "geomorphShapes")) {
    Y <- A$landmarks
    if(any(unlist(lapply(Y, is.na)))) stop("Data matrix contains missing values. Estimate these first (see 'estimate.missing').")
    curves <- A$curves
    n <- A$n
    p <- A$p
    k <- A$k
    
    spec.names <- names(Y)
    p.names <- dimnames(Y[[1]])[[1]]
    k.names <- c("X", "Y", "Z")[1:k] 
    
  } else {
    
    if(!is.array(A)) stop("Coordinates must be a 3D array")
    if(length(dim(A)) != 3) stop("Coordinates array does not have proper dimensions")
    if(any(is.na(A))) stop("Data matrix contains missing values. Estimate these first (see 'estimate.missing').")
    n <- dim(A)[[3]]; p <- dim(A)[[1]]; k <- dim(A)[[2]]
    
    spec.names <- 1:n
    p.names <- 1:p
    k.names <- c("X", "Y", "Z")[1:k] 
    
    dim.names <- dimnames(A)
    if(length(dim.names) != 0) {
      dim.name.check <- sapply(1:length(dim.names), is.null)
      if(!dim.name.check[[1]]) spec.names <- dim.names[[1]]
      if(!dim.name.check[[2]]) spec.names <- dim.names[[2]]
      if(!dim.name.check[[3]]) spec.names <- dim.names[[3]]
    }
    
    Y <- lapply(1:n, function(j) A[,,j])
  }
  
  if(!is.logical(ProcD)) prD <- TRUE else prD <- ProcD
  if(is.null(max.iter)) max.it <- 10 else max.it <- as.numeric(max.iter)
  if(is.numeric(max.it) & max.it > 100) {
    warning("GPA might be halted ahead of maximum iterations, 
            as the number chosen is exceedingly large", immediate. = TRUE,
            call. = FALSE, noBreaks. = TRUE)
    max.it = 100
  }
  if(is.na(max.it)) max.it <- 10
  if(max.it < 0) max.it <- 10
  if(!is.null(curves)) {
    curves <- as.matrix(curves) 
    if(ncol(curves) != 3) stop("curves must be a matrix of three columns")
  } else curves <- NULL
  if(!is.null(surfaces)) surf <- as.vector(surfaces) else surf <- NULL
  
  fixed <- 1:p
  if(!is.null(curves)) fixed <- fixed[which(!fixed %in% curves[,2])]
  if(!is.null(surf)) fixed <- fixed[which(!fixed %in% surf)]
  if(length(fixed) == 0) fixed <- 1:p
  if(is.null(rot.pts)) rot.pts <- 1:p
  if(length(rot.pts) == 1 && rot.pts  == "all") rot.pts <- 1:p
  if(length(rot.pts) == 1 && rot.pts  == "fixed") rot.pts <- fixed
  if(length(rot.pts) == 1 && rot.pts  == "sliders") rot.pts <- setdiff(1:p, fixed)
  if(length(rot.pts) > p) rot.pts <- 1:p
  if(length(rot.pts) < k) rot.pts <- 1:p
  
  if(print.progress && Parallel != FALSE) {
        print.progress = FALSE
        message("\nResults status turned off for parallel processing.\n")
      }  
  
  
  if(print.progress == TRUE){
    cat("\nPerforming GPA\n")
    if(!is.null(curves) || !is.null(surf)) 
      gpa <- pGpa.wSliders(Y, curves = curves, surf=surf, rot.pts = rot.pts,
                           PrinAxes = PrinAxes, max.iter=max.it, 
                           tol = tol, appBE = approxBE, sen = sen, ProcD=prD) else
                             gpa <- pGpa(Y, rot.pts = rot.pts, PrinAxes = PrinAxes, 
                                         max.iter=max.it, tol = tol)
                                                                
  } else {
    if(!is.null(curves) || !is.null(surf)) 
      gpa <- .pGpa.wSliders(Y, curves = curves, surf=surf, rot.pts = rot.pts,
                            PrinAxes = PrinAxes, max.iter=max.it, 
                            tol = tol, appBE= approxBE, sen = sen, ProcD=prD, 
                            Parallel = Parallel) else
                              gpa <- .pGpa(Y, rot.pts = rot.pts, PrinAxes = PrinAxes, 
                                           max.iter=max.it, tol = tol, 
                                           Parallel = Parallel)
  }
  
  coords <- gpa$coords
  M <- gpa$consensus
  dimnames(M) <- list(p.names, k.names)
  if(is.null(colnames(M))) colnames(M) <- c("X", "Y", "Z")[1:k] 
  
  if (Proj == TRUE) {
    if(print.progress) cat("\nMaking projections... ")
    coords <- orp(coords)
    M <- Reduce("+",coords)/n
    dimnames(M) <- list(p.names, k.names)
    if(print.progress) cat("Finished!\n")
  }
  Csize <- gpa$CS
  names(Csize) <- spec.names
  iter <- gpa$iter
  coordsL <- coords
  coords <- simplify2array(coords)
  dimnames(coords) <- list(p.names, k.names, spec.names)
  
  
  if(!is.null(curves) || !is.null(surf)) {
    nsliders <- nrow(curves)
    nsurf <- length(surf)
    if(ProcD == TRUE) smeth <- "ProcD" else smeth <- "BE"
  } else {
    nsliders <- 0
    nsurf <- 0
    smeth <- NULL
  }
  if(is.null(nsliders)) nsliders <- 0; if(is.null(nsurf)) nsurf <- 0
  
  if(verbose) {
    if(print.progress) cat("\nCompiling additional results ... ")
    
    two.d.coords <- two.d.array(coords)
    rownames(two.d.coords) <- spec.names
    dat <- data.frame(coords = two.d.coords, Csize = Csize)
    pt.var <- Reduce("+", Map(function(y) y^2/n, coordsL))
    rownames(pt.var) <- p.names
    colnames(pt.var) <- c("Var.X", "Var.Y", "Var.Z")[1:k]
    pt.VCV <- var(two.d.coords)
    procD <- try(dist(two.d.coords), silent = TRUE)
    if(inherits(procD, "try-error")) procD <- NULL
    if(print.progress) cat("Finished!\n")
  } else pt.var <- pt.VCV <- procD <- two.d.coords <- dat <- NULL
  

  out <- list(coords=coords, Csize=Csize, 
              iter=iter, rot.pts = rot.pts,
              consensus = M, procD = procD, 
              p = p,k = k, nsliders=nsliders, nsurf = nsurf,
              points.VCV = pt.VCV, points.var = pt.var, 
              data = dat,
              Q = gpa$Q, slide.method = smeth, call= match.call())
  class(out) <- "gpagen"
  

  out
}
