#' Pairwise Comparisons of Slopes
#'
#' Function performs pairwise comparisons among slopes for groups as specified by a linear model.
#'
#' The function performs pairwise comparisons to identify differences in slopes between groups. The function is 
#' designed as a post-hoc test to MANCOVA, where the latter has identified significant shape variation explained by a 
#' covariate*group interaction term. 
#' 
#'  As input the user provides a formula describing the linear model of how shape varies as a function of several explanatory 
#'  variables. This MUST be in the form of: [y~covariate + group], and the shape data (y) must be in the form of a 
#'  two-dimensional data matrix of dimension 
#'  (n x [p x k]), rather than a 3D array.  It is assumed that the landmarks have previously 
#'   been aligned using Generalized Procrustes Analysis (GPA) [e.g., with \code{\link{gpagen}}]. The function
#'   \code{\link{two.d.array}} can be used to obtain a two-dimensional data matrix from a 3D array of landmark
#'   coordinates. From the data, the slopes for each group are estimated, and pairwise differences in slopes determined.
#'   
#'   For the model, one can specify whether slopes or intercepts are to be evaluated. Slopes are compared if (heterogenous slopes) 
#'   het.slopes=TRUE.To evaluate significance of the pairwise differences, two possible resampling procedures are provided. First, if 
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
#' @param het.slopes A logical value indicting whether slopes are to be compared
#' @param angle.type A value specifying whether differences between slopes should be represented by vector
#' correlations (r), radians (rad) or degrees (deg)
#' @param RRPP a logical value indicating whether residual randomization should be used for significance testing
#' @keywords analysis
#' @export
#' @author Mike Collyer
#' @references Anderson MJ. and C.J.F. terBraak. 2003. Permutation tests for multi-factorial analysis of variance.
#'    Journal of Statistical Copmutation and Simulation 73: 85-113.
#' @references Collyer, M.L., D.J. Sekora, and D.C. Adams. 2014. A method for analysis of phenotypic change for phenotypes described 
#' by high-dimensional data. Heredity. (In Press).
#' @return Function returns a list with the following components: 
#'   \item{ANOVA.table}{An ANOVA table assessing the linear model}
#'   \item{Obs.LS.Dist}{A matrix of pairwise differences between intercepts (least squares means) if het.slopes = FALSE }
#'   \item{Slope.Dist}{A matrix of pairwise differences between slopes represented as correlations, radians, or degrees, if het.slopes = TRUE}
#'   \item{Prob.Dist}{A matrix of pairwise significance levels based on permutation for either intercepts or slopes}
#'   \item{Magnitude.Diff}{A matrix of pairwise differences in magnitude between regression lines (from smallest to largest specimen, if het.slopes = TRUE)}
#'   \item{Prob.Mag}{A matrix of pairwise significance levels based on permutation for magnitude differences}
#'   @examples
#' ### MANCOVA example for Goodall's F test (multivariate shape vs. factors)
#' data(plethodon) 
#' Y.gpa<-gpagen(plethodon$land)    #GPA-alignment    
#' y<-two.d.array(Y.gpa$coords)
#' 
#' ## Pairwise slope test
#' # Assuming heterogenous slopes
#' pairwise.slope.test(y~Y.gpa$Csize+plethodon$site,iter=49,angle.type="rad")
#' 
#' # Assuming parallel slopes
#' pairwise.slope.test(y~Y.gpa$Csize+plethodon$site,het.slopes=FALSE, iter=49, angle.type="rad") 
#' 
#' ## Using RRPP
#' # Assuming heterogenous slopes
#' pairwise.slope.test(y~Y.gpa$Csize+plethodon$site,iter=49, angle.type="rad", RRPP=TRUE)
#' # Assuming parallel slopes
#' pairwise.slope.test(y~Y.gpa$Csize+plethodon$site, het.slopes=FALSE, 
#'       iter=49, angle.type="rad", RRPP=TRUE)
pairwise.slope.test <- function (f1, iter = 999, het.slopes = T, angle.type = "r", RRPP = FALSE) {
  data = NULL
  form.in <- formula(f1)
  Terms <- terms(form.in, keep.order = T)
  Term.labels = c(attr(Terms, "term.labels"), paste(attr(Terms, "term.labels")[1], attr(Terms, "term.labels")[2], sep=":"))
  Y <- eval(form.in[[2]], parent.frame())
  X <- model.matrix(Terms)
  j <- ncol(attr(Terms, "factors"))
  dat <- model.frame(form.in, data=NULL)
  Xdims <- dim(X)
  
  if (length(dim(Y)) != 2) {
    stop("\nResponse matrix (shape) not a 2D array. Use 'two.d.array' first.")
  }
  if (any(is.na(Y)) == T) {
    stop("\nResponse data matrix (shape) contains missing values. Estimate these first (see 'estimate.missing').")
  } 
  if(any(angle.type == c("r", "deg","rad")) == FALSE){
  	print("angle.type not one of r, deg, or rad; assuming angle.type = r")
  	angle.type = "r"
  }
  
  if (j < 2) stop ("\nFormula must contain at least one covariate and one factor")
  
  if (j > 2) stop ("\n Formula can only contain one covariate and one factor, in that order")
  
  if (class(dat[,2]) != "numeric") stop("\nFirst variable in formula must be a covariate")
  if (class(dat[,3]) != "factor") stop("\nSecond variable in formula must be a factor")
  
  if(het.slopes == FALSE){
    g <- Xdims[2]-2
    SS.tmp <- numeric(j)
    Xs <- array(0, c(Xdims, 3))
    Xs[,1,1] <- 1
    for(i in 1:2){
      x <- as.matrix(model.matrix(Terms[1:i], data = dat))
      Xs[,1:ncol(x),(i+1)] <- x
      SS.tmp[i] <- SSE(lm(Y ~ x -1))
    }
    df <- c(1,g)
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
    
    Xm <- rbind(0,diag(1,g))
    Xm <- cbind(1, mean(X[,2]),Xm)
    m <- Xm%*%coef(lm(Y~X-1))
    rownames(m) <- levels(dat[,3])
    
    P <- array(0,c(dim(matrix(SS.obs)),iter+1))
    P[,,1]=SS.obs
    D <- array(0,c(nrow(m),nrow(m),iter+1))
    D[,,1] <- as.matrix(dist(m))
    dimnames(D)[1:2] <- dimnames(as.matrix(dist(m)))
    
    if(RRPP==TRUE){
      for(i in 1:iter){
        Yr <- RRP.submodels(Xs,Y)
        for (ii in 1:2) {
          SS.tmp[ii] <- SSE(lm(Yr[,,ii] ~ Xs[,,ii+1] -1))				
        }
        P[,,i+1] <- SS.null-SS.tmp
        D[,,i+1]	 <- as.matrix(dist(Xm%*%coef(lm(Yr[,,2]~ Xs[,,3] -1))))
      }
      P.val <- Pval.matrix(P)
      D.p.val <- Pval.matrix(D)
    }
    
    
    if(RRPP==FALSE){
      for(i in 1:iter){
        Yr <- Y[sample(nrow(Y)),]
        for (ii in 1:2) {
          SS.tmp[ii] <- SSE(lm(Yr ~ Xs[,,ii+1] -1))
          SS.null <- c(SSE(lm(Yr ~ 1)), SS.tmp)[1:2]			
        }
        P[,,i+1] <- SS.null-SS.tmp
        D[,,i+1] <- as.matrix(dist(Xm%*%coef(lm(Yr ~ Xs[,,3] -1))))		
      }
      P.val <- Pval.matrix(P)
      D.p.val <- Pval.matrix(D)
    }	
    
    anova.tab <- data.frame(df = c(df,df.res,df.tot), 
                            SS = c(SS.obs, SS.res, SS.tot), 
                            MS = c(MS, MS.res, MS.tot),
                            Rsq = c(Rsq, NA, NA),
                            F = c(F, NA, NA),
                            P.val = c(P.val, NA, NA))
    rownames(anova.tab) <- c(attr(Terms, "term.labels")[1:2], "Residuals","Total")
    if(RRPP == TRUE) anova.title = "\nRandomized Residual Permutation Procedure used\n"
    if(RRPP == FALSE) anova.title = "\nRandomization of Raw Values used\n"
    attr(anova.tab, "heading") <- paste("\nType I (Sequential) Sums of Squares and Cross-products\n",anova.title)
    class(anova.tab) <- c("anova", class(anova.tab))
    dm <- data.frame(as.matrix(D[,,1]))
    result <- list(anova.tab=anova.tab, Obs.LS.dist = dm, Prob.dist = D.p.val)
  }
  
  if(het.slopes == T){
    g <- ncol(model.matrix(Terms[1:2]))
    Xdims[2] = Xdims[2] + Xdims[2] - 2
    SS.tmp <- numeric(3)
    Xs <- array(0, c(Xdims, 4))
    Xs[,1,1] <- 1
    for(i in 1:2){
      x <- as.matrix(model.matrix(Terms[1:i]))
      Xs[,1:ncol(x),(i+1)] <- x
      mod.mat <- model.matrix(Terms[1:i], data = dat)
      SS.tmp[i] <- SSE(lm(Y ~ mod.mat-1))
    }
    Xs[,,4] <- cbind(X, X[,2]*X[,-(1:2)])
    SS.tmp[3] <- SSE(lm(Y ~ Xs[,,4]-1))
    df <- c(1,g-2,g-2)
    SS.null <- (c(SSE(lm(Y~1)),SS.tmp))[1:3]
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
    Xm <- rbind(0,diag(1,g-2))
    Xm <- cbind(1, mean(X[,2]),Xm, mean(X[,2])*Xm)
    B = coef(lm(Y~Xs[,,4]-1))
    m = Xm%*%B
    rownames(m) <- levels(dat[,3])
    Bslopes = rbind(B[2,], B[2,]+B[-(1:g),])
    P <- array(0,c(dim(matrix(SS.obs)),iter+1))
    P[,,1]=SS.obs
    D <- V <- array(0,c(nrow(Bslopes),nrow(Bslopes),iter+1))
    D[,,1] <- as.matrix(dist(Bslopes))
    dimnames(D)[1:2] <- dimnames(as.matrix(dist(m)))
    if(angle.type == "r") V[,,1] = vec.cor.matrix(Bslopes) 
    if(angle.type == "rad") V[,,1] = vec.ang.matrix(Bslopes)  
    if(angle.type == "deg") V[,,1] = vec.ang.matrix(Bslopes)*180/pi   
    
    if(RRPP==TRUE){ 
      for(i in 1:iter){
        Yr <- RRP.submodels(Xs,Y)
        for (ii in 1:3) {
          SS.tmp[ii] <- SSE(lm(Yr[,,ii] ~ Xs[,,ii+1] -1))				
        }
        P[,,i+1] <- SS.null-SS.tmp
        Br = coef(lm(Yr[,,3]~Xs[,,4]-1))
        Brslopes = rbind(Br[2,], Br[2,]+Br[-(1:g),])
        D[,,i+1] <- as.matrix(dist(Brslopes))
        if(angle.type == "rad") {
          V[,,i+1] <- vec.ang.matrix(Brslopes)
        } else { 
          if(angle.type == "deg"){
            V[,,i+1] <- vec.ang.matrix(Brslopes)*180/pi
          } else {V[,,i+1] <- vec.cor.matrix(Brslopes)	}
        }}
      
      P.val <- Pval.matrix(P)
      D.p.val <- Pval.matrix(D)
      V.p.val <- Pval.matrix(abs(V))
    }
    
    
    if(RRPP==FALSE){
      for(i in 1:iter){
        Yr <- Y[sample(nrow(Y)),]
        for (ii in 1:3) {
          SS.tmp[ii] <- SSE(lm(Yr ~ Xs[,,ii+1] -1))
          SS.null <- (c(SSE(lm(Y~1)),SS.tmp))[1:3]				
        }
        P[,,i+1] <- SS.null-SS.tmp
        Br = coef(lm(Yr~Xs[,,4]-1))
        Brslopes = rbind(Br[2,], Br[2,]+Br[-(1:g),])
        D[,,i+1] <- as.matrix(dist(Brslopes))
        if(angle.type == "rad") {
          V[,,i+1] <- vec.ang.matrix(Brslopes)
        } else { 
          if(angle.type == "deg"){
            V[,,i+1] <- vec.ang.matrix(Brslopes)*180/pi
          } else {V[,,i+1] <- vec.cor.matrix(Brslopes)	}
        }}
      P.val <- Pval.matrix(P)
      D.p.val <- Pval.matrix(D)
      V.p.val <- Pval.matrix(abs(V))
    }	
    
    if(angle.type == "r")	{
      anova.tab <- data.frame(df = c(df,df.res,df.tot), 
                              SS = c(SS.obs, SS.res, SS.tot), 
                              MS = c(MS, MS.res, MS.tot),
                              Rsq = c(Rsq, NA, NA),
                              F = c(F, NA, NA),
                              P.val = c(P.val, NA, NA))
      rownames(anova.tab) <- c(Term.labels, "Residuals","Total")
      if(RRPP == TRUE) anova.title = "\nRandomized Residual Permutation Procedure used\n"
      if(RRPP == FALSE) anova.title = "\nRandomization of Raw Values used\n"
      attr(anova.tab, "heading") <- paste("\nType I (Sequential) Sums of Squares and Cross-products\n",anova.title)
      class(anova.tab) <- c("anova", class(anova.tab))
      dm <- data.frame(as.matrix(D[,,1]))
      result <- list(anova.tab=anova.tab, Vector.magnitude.difference = dm, VM.prob.dist = D.p.val, Vector.correlation = V[,,1], VC.prob.dist = V.p.val)
    } 
    else {
      anova.tab <- data.frame(df = c(df,df.res,df.tot), 
                              SS = c(SS.obs, SS.res, SS.tot), 
                              MS = c(MS, MS.res, MS.tot),
                              Rsq = c(Rsq, NA, NA),
                              F = c(F, NA, NA),
                              P.val = c(P.val, NA, NA))
      rownames(anova.tab) <- c(Term.labels, "Residuals","Total")
      if(RRPP == TRUE) anova.title = "\nRandomized Residual Permutation Procedure used\n"
      if(RRPP == FALSE) anova.title = "\nRandomization of Raw Values used\n"
      attr(anova.tab, "heading") <- paste("\nType I (Sequential) Sums of Squares and Cross-products\n",anova.title)
      class(anova.tab) <- c("anova", class(anova.tab))
      dm <- data.frame(as.matrix(D[,,1]))
      result <- list(anova.tab=anova.tab, Vector.magnitude.difference = dm, VM.prob.dist = D.p.val, Angle = V[,,1], Angle.prob.dist = V.p.val)
    }
    
  }
  result
}