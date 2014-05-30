#' Procrustes ANOVA/regression for shape data
#'
#' Function performs Procrustes ANOVA with permutation procedures to assess statistical hypotheses describing 
#'   patterns of shape variation and covariation for a set of Procrustes-aligned coordinates
#'
#' The function quantifies the relative amount of shape variation attributable to one or more factors in a 
#'   linear model and assesses this variation via permutation. Data input is specified by a formula (e.g., 
#'   y~X), where 'y' specifies the response variables (shape data), and 'X' contains one or more independent 
#'   variables (discrete or continuous). The response matrix 'y' must be in the form of a two-dimensional data 
#'   matrix of dimension (n x [p x k]), rather than a 3D array.  It is assumed that the landmarks have previously 
#'   been aligned using Generalized Procrustes Analysis (GPA) [e.g., with \code{\link{gpagen}}]. The function
#'   \code{\link{two.d.array}} can be used to obtain a two-dimensional data matrix from a 3D array of landmark
#'   coordinates. The names specified for the independent (x) variables in the formula represent one or more 
#'   vectors containing continuous data or factors. It is assumed that the order of the specimens in the 
#'   shape matrix matches the order of values in the independent variables.
#'
#'   The function performs statistical assessment of the terms in the model using Procrustes distances among 
#'   specimens, rather than explained covariance matrices among variables. With this approach, the sum-of-squared 
#'   Procrustes distances are used as a measure of SS (see Goodall 1991). The observed SS are evaluated through 
#'   permutation. In morphometrics this approach is known as a Procrustes ANOVA (Goodall 1991), which is equivalent
#'   to distance-based anova designs (Anderson 2001). Two possible resampling procedures are provided. First, if RRPP=FALSE, 
#'   the rows of the matrix of shape variables 
#'   are randomized relative to the design matrix. This is analogous to a 'full' randomization. Second, if RRPP=TRUE,
#'   a residual randomization permutation procedure is utilized (Collyer et al. 2014). Here, residual shape values from a reduced model are
#'   obtained, and are randomized with respect to the linear model under consideration. These are then added to 
#'   predicted values from the remaining effects to obtain pseudo-values from which SS are calculated. NOTE: for
#'   single-factor designs, the two approaches are identical.  However, when evaluating factorial models it has been
#'   shown that RRPP attains higher statistical power and thus has greater ability to identify patterns in data should
#'   they be present (see Anderson and terBraak 2003). 
#'
#' @param f1 A formula for the linear model (e.g., y~x1+x2)
#' @param iter Number of iterations for significance testing
#' @param RRPP a logical value indicating whether residual randomization should be used for significance testing
#' @param verbose A logical value specifying whether additional output should be displayed
#' @keywords analysis
#' @export
#' @author Dean Adams and Mike Collyer
#' @return Function returns an ANOVA table of statistical results for all factors: df (for each factor), SS, MS,
#' F ratio, Prand, and Rsquare.
#' @references Anderson MJ. 2001. A new method for non-parametric multivariate analysis of variance. 
#'    Austral Ecology 26: 32-46.
#' @references Anderson MJ. and C.J.F. terBraak. 2003. Permutation tests for multi-factorial analysis of variance.
#'    Journal of Statistical Copmutation and Simulation 73: 85-113.
#' @references Collyer, M.L., D.J. Sekora, and D.C. Adams. 2014. A method for analysis of phenotypic change for phenotypes described 
#' by high-dimensional data. Heredity. (In Press).
#' @references Goodall, C. R. 1991. Procrustes methods in the statistical analysis of shape. Journal of the 
#'    Royal Statistical Society B 53:285-339.
#' @examples
#' ### MANOVA example for Goodall's F test (multivariate shape vs. factors)
#' data(plethodon) 
#' Y.gpa<-gpagen(plethodon$land)    #GPA-alignment    
#' y<-two.d.array(Y.gpa$coords)
#'
#' procD.lm(y~plethodon$species*plethodon$site,iter=99)
#'
#' ### Regression example
#' data(ratland)
#' rat.gpa<-gpagen(ratland)         #GPA-alignment
#'
#' procD.lm(two.d.array(rat.gpa$coords)~rat.gpa$Csize,iter=99)
#' 
#' ## using RRPP
#'  procD.lm(two.d.array(rat.gpa$coords)~rat.gpa$Csize,iter=99,RRPP=TRUE)
procD.lm <- function(f1, iter = 999, RRPP = FALSE, verbose=FALSE){ 
  data=NULL
	form.in <- formula(f1)
    	Terms <- terms(form.in, keep.order = TRUE)
    	Y <- eval(form.in[[2]], parent.frame())
    	if (length(dim(Y)) != 2) {
        stop("Response matrix (shape) not a 2D array. Use 'two.d.array' first.")
    	}
    	if (any(is.na(Y)) == T) {
        stop("Response data matrix (shape) contains missing values. Estimate these first (see 'estimate.missing').")
    }
    	if (is.null(dimnames(Y)[[1]])) {
        print("No specimen names in response matrix. Assuming specimens in same order.")
    }
    	df <- df.tmp <- SS.tmp <- SS.obs <- F <- Rsq <- array()
    	dat <- model.frame(form.in, data=NULL)
	Xdims <- dim(as.matrix(model.matrix(Terms)))
	j <- ncol(attr(Terms, "factors"))
	
	if(j==1){
		SS.null <- SSE(lm(Y~1))
		mod.mat <- model.matrix(Terms[1], data = dat)
		SS.tmp <- SSE(lm(Y~mod.mat-1))
		SS.obs <- SS.null-SS.tmp
		SS.tot <- SS.null
		SS.res <- SS.tot-SS.obs
		df <- ncol(mod.mat)-1
		MS <- SS.obs/df
		df.tot <- nrow(Y) - 1
        df.res <- nrow(Y) - df -1
        MS.tot <- SS.tot/df.tot
        MS.res <- SS.res/df.res
        Rsq <- SS.obs/SS.tot
        F <- MS/MS.res
        P <- array(c(SS.obs,rep(0,iter)))
    	for(i in 1:iter){
			Yr <- Y[sample(nrow(Y)),]
    		SS.tmp <- SSE(lm(Yr~mod.mat-1))
    		P[i+1] <- SS.null-SS.tmp		
		}
    	P.val <- pval(P)
	}
	
	if(j>1){
	Xs <- array(0, c(Xdims, (j+1)))
	Xs[,1,1] <- 1
	for(i in 1:j){
		x <- as.matrix(model.matrix(Terms[1:i]))
		Xs[,1:ncol(x),(i+1)] <- x
        mod.mat <- model.matrix(Terms[1:i], data = dat)
        SS.tmp[i] <- SSE(lm(Y ~ mod.mat-1))
        df.tmp[i] <- ifelse(ncol(mod.mat) == 1, 1, (ncol(mod.mat) - 1))
        ifelse(i == 1, df[i] <- df.tmp[i], df[i] <- (df.tmp[i] - df.tmp[i - 1]))
    	}
    SS.null <- (c(SSE(lm(Y~1)),SS.tmp))[1:j]
    SS.obs <- SS.null - SS.tmp
    MS <- SS.obs/df
    SS.tot <- SSE(lm(Y~1))
    SS.res <- SS.tot - sum(SS.obs)
    df.tot <- nrow(Y) - 1
    df.res <- nrow(Y) - 1 - sum(df)
    MS.tot <- SS.tot/df.tot
    MS.res <- SS.res/df.res
    Rsq <- SS.obs/SS.tot
    F <- MS/MS.res
    P <- array(0,c(dim(matrix(SS.obs)),iter+1))
    P[,,1]=SS.obs
    
    if(RRPP==TRUE){
    
    		for(i in 1:iter){
			Yr <- RRP.submodels(Xs,Y)
    			for (ii in 1:j) {
        			SS.tmp[ii] <- SSE(lm(Yr[,,ii] ~ Xs[,,ii+1] -1))				
			}
    	  	P[,,i+1] <- SS.null-SS.tmp		
    		}
		P.val <- Pval.matrix(P)
		}
    
        if(RRPP==FALSE){
    
    		for(i in 1:iter){
			Yr <- Y[sample(nrow(Y)),]
    			for (ii in 1:j) {
        			SS.tmp[ii] <- SSE(lm(Yr ~ Xs[,,ii+1] -1))				
			}
    	  	P[,,i+1] <- SS.null-SS.tmp		
    		}
		P.val <- Pval.matrix(P)
		}
    }
    
    anova.tab <- data.frame(df = c(df,df.res,df.tot), 
    	SS = c(SS.obs, SS.res, SS.tot), 
    	MS = c(MS, MS.res, MS.tot),
    	Rsq = c(Rsq, NA, NA),
    	F = c(F, NA, NA),
    	P.val = c(P.val, NA, NA))
    	rownames(anova.tab) <- c(attr(Terms, "term.labels"), "Residuals","Total")
    	if(RRPP == TRUE) anova.title = "\nRandomized Residual Permutation Procedure used\n"
    	if(RRPP == FALSE) anova.title = "\nRandomization of Raw Values used\n"
    attr(anova.tab, "heading") <- paste("\nType I (Sequential) Sums of Squares and Cross-products\n",anova.title)
    class(anova.tab) <- c("anova", class(anova.tab))
    if(verbose==TRUE)  {
    	out <- list(call=match.call(), anova.tab = anova.tab, SS.rand = P, model.matrix = mod.mat, terms = Terms)
    	class(out) <- "prodcD.lm"
    	out} else {
    		return(anova.tab)
    	}
}  