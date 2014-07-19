#' Pairwise Group Comparisons
#'
#' Function performs pairwise comparisons among groups using the Euclidean distances among group means.
#'
#' The function performs pairwise comparisons to identify shape differences among groups. The function is designed as a post-hoc
#'  test to Procrustes ANOVA, where the latter has identified significant shape variation explained by a grouping factor. 
#'  
#'  As input the user provides a formula describing the linear model of how shape (y) varies as a function of groups (x). 
#'  Multiple factors may be examined. The shape data (y) must be in the form of a two-dimensional data matrix of dimension 
#'  (n x [p x k]), rather than a 3D array.  It is assumed that the landmarks have previously 
#'   been aligned using Generalized Procrustes Analysis (GPA) [e.g., with \code{\link{gpagen}}]. The function
#'   \code{\link{two.d.array}} can be used to obtain a two-dimensional data matrix from a 3D array of landmark
#'   coordinates. From the data, the Euclidean distances among group means are estimated, and used as test values.
#'   
#'   To evaluate significance of group differences, two possible resampling procedures are provided. First, if 
#'   RRPP=FALSE, the rows of the matrix of shape variables are randomized relative to the design matrix. This is 
#'   analogous to a 'full' randomization. Second, if RRPP=TRUE, a residual randomization permutation procedure 
#'   is utilized (Collyer et al. 2014). Here, residual shape values from a reduced model are
#'   obtained, and are randomized with respect to the linear model under consideration. These are then added to 
#'   predicted values from the remaining effects to obtain pseudo-values from which SS are calculated. NOTE: for
#'   single-factor designs, the two approaches are identical.  However, when evaluating factorial models it has been
#'   shown that RRPP attains higher statistical power and thus has greater ability to identify patterns in data should
#'   they be present (see Anderson and terBraak 2003). 
#'
#' @param f1 A formula for the linear model from which groups are to be compared (e.g., y~x1+x2)
#' @param iter Number of iterations for permutation test
#' @param RRPP a logical value indicating whether residual randomization should be used for significance testing
#' @keywords analysis
#' @export
#' @author Mike Collyer and Dean Adams
#' @references Anderson MJ. and C.J.F. terBraak. 2003. Permutation tests for multi-factorial analysis of variance.
#'    Journal of Statistical Copmutation and Simulation 73: 85-113.
#' @references Collyer, M.L., D.J. Sekora, and D.C. Adams. 2014. A method for analysis of phenotypic change for phenotypes described 
#' by high-dimensional data. Heredity. (In Press).
#' @return Function returns a list with the following components: 
#'   \item{Dist.obs}{A matrix of Euclidean distances among group means}
#'   \item{Prob.Dist}{A matrix of pairwise significance levels based on permutation}
#' @examples
#' data(plethodon) 
#' Y.gpa<-gpagen(plethodon$land)    #GPA-alignment    
#' y<-two.d.array(Y.gpa$coords)
#' ### Procrustes ANOVA
#' procD.lm(y~plethodon$species,iter=99)
#'
#' ### Pairwise comparisons: full randomization
#' pairwiseD.test(y~plethodon$species*plethodon$site,iter=99)
#' 
#' ## Pairwise comparisons: residual randomization
#' #' pairwiseD.test(y~plethodon$species*plethodon$site,iter=99,RRPP=TRUE)
pairwiseD.test <- function (f1, iter = 999, RRPP = FALSE) {
	data=NULL
	form.in <- formula(f1)
    Terms <- terms(form.in, keep.order = T)
    Y <- eval(form.in[[2]], parent.frame())
    X <- model.matrix(Terms)
    j <- ncol(attr(Terms, "factors"))
    dat <- model.frame(form.in, data=NULL)
    Xdims <- dim(X)
    df <- df.tmp <- SS.tmp <- SS.obs <- F <- Rsq <- array()
    	
    if(any(X!=1 & X!=0)) 
    stop("pairwise D test only meaningful for factors\nConsider using pairwise.slope.test for covariates")
    
    if (length(dim(Y)) != 2) {
        stop("Response matrix (shape) not a 2D array. Use 'two.d.array' first.")
    }
    if (any(is.na(Y)) == T) {
        stop("Response data matrix (shape) contains missing values. Estimate these first (see 'estimate.missing').")
    }
 	
 	if(j==1){
		SS.null <- SSE(lm(Y~1))
		fit <- lm(Y~X-1)
		SS.tmp <- SSE(fit)
		SS.obs <- SS.null-SS.tmp
		SS.tot <- SS.null
		SS.res <- SS.tot-SS.obs
		df <- ncol(X)-1
		MS <- SS.obs/df
		df.tot <- nrow(Y) - 1
        df.res <- nrow(Y) - df -1
        MS.tot <- SS.tot/df.tot
        MS.res <- SS.res/df.res
        Rsq <- SS.obs/SS.tot
        F <- MS/MS.res
        Xm <-diag(1,Xdims[2])
        Xm[,1] <- 1
        m <- Xm%*%coef(fit)
        rownames(m) <-levels(model.frame(Terms)[,2])
        P <- array(c(SS.obs,rep(0,iter)))
        D <- array(0,c(nrow(m),nrow(m),iter+1))
        D[,,1] <- as.matrix(dist(m))
    	for(i in 1:iter){
			Yr <- Y[sample(nrow(Y)),]
			fit.r <- lm(Yr~X-1)
    		SS.tmp <- SSE(fit.r)
    		P[i+1] <- SS.null-SS.tmp
    		mr <- Xm%*%coef(fit.r)
    		D[,,i+1] <-as.matrix(dist(mr))
		}
		dimnames(D) <- dimnames(as.matrix(dist(m)))
    		P.val <- pval(P)
    		D.p.val <- Pval.matrix(D)
 
	}
 
 	if(j>1){
 
		Xs <- array(0, c(Xdims, (j+1)))
		Xs[,1,1] <- 1
		for(i in 1:j){
			x <- as.matrix(model.matrix(Terms[1:i], data=dat))
			Xs[,1:ncol(x),(i+1)] <- x
        	SS.tmp[i] <- SSE(lm(Y ~ x -1))
        	df.tmp[i] <- ifelse(ncol(x) == 1, 1, (ncol(x) - 1))
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
    	Xm <- unique(X[do.call(order, lapply(1:ncol(X), function(i) X[, i])), ])
    	m <- Xm%*%coef(lm(Y~X-1))
    	nm=(model.frame(Terms))[,2]
        for(i in 3:j){nm <- factor(paste(nm,(model.frame(Terms))[,j]))}
        rownames(m)<-levels(nm)
    	P <- array(0,c(dim(matrix(SS.obs)),iter+1))
    	P[,,1]=SS.obs
    	D <- array(0,c(nrow(m),nrow(m),iter+1))
    	D[,,1] <- as.matrix(dist(m))
    	dimnames(D)[1:2] <- dimnames(as.matrix(dist(m)))
    
    if(RRPP==TRUE){
    		for(i in 1:iter){
			Yr <- RRP.submodels(Xs,Y)
    			for (ii in 1:j) {
        			SS.tmp[ii] <- SSE(lm(Yr[,,ii] ~ Xs[,,ii+1] -1))				
				}
    	  	P[,,i+1] <- SS.null-SS.tmp
    	  	D[,,i+1] <- as.matrix(dist(Xm%*%coef(lm(Yr[,,j]~ X -1))))
    		}
		P.val <- Pval.matrix(P)
		D.p.val <- Pval.matrix(D)
		}
    
	if(RRPP==FALSE){
    		for(i in 1:iter){
			Yr <- Y[sample(nrow(Y)),]
    			for (ii in 1:j) {
        			SS.tmp[ii] <- SSE(lm(Yr ~ Xs[,,ii+1] -1))	
        			SS.null <- (c(SSE(lm(Y~1)),SS.tmp))[1:j]			
				}
    	  	P[,,i+1] <- SS.null-SS.tmp
    	  	D[,,i+1] <- as.matrix(dist(Xm%*%coef(lm(Yr~ X -1))))		
    		}
		P.val <- Pval.matrix(P)
		D.p.val <- Pval.matrix(D)
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
 	
 	dm <- data.frame(as.matrix(D[,,1]))
 	
    list(anova.tab=anova.tab, Obs.dist = dm, Prob.dist = D.p.val)
}
