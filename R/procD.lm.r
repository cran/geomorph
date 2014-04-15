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
#'   permutation, where the rows of the shape matrix are randomized relative to the design matrix. Procedurally, 
#'   Procrustes ANOVA is identical to permutational-MANOVA as used in other fields (Anderson 2001). For several 
#'   reasons, Procrustes ANOVA is particularly useful for shape data. First, covariance matrices from GPA-aligned 
#'   Procrustes coordinates are singular, and thus standard approaches such as MANOVA cannot be accomplished 
#'   unless generalized inverses are utilized. This problem is accentuated when using sliding semilandmarks. 
#'   Additionally, GM datasets often have more variables than specimens (the 'small N large P' problem). In 
#'   these cases, distance-based procedures can still be utilized to assess statistical hypotheses, whereas 
#'   standard linear models cannot. 
#'
#' @param f1 A formula for the linear model (e.g., y~x1+x2)
#' @param data An optional value specifying a data frame containing all data (not required)
#' @param iter Number of iterations for significance testing
#' @keywords analysis
#' @export
#' @author Dean Adams
#' @return Function returns an ANOVA table of statistical results for all factors: df (for each factor), SS, MS,
#' F ratio, Prand, and Rsquare.
#' @references Anderson MJ. 2001. A new method for non-parametric multivariate analysis of variance. 
#'    Austral Ecology 26: 32-46.
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
procD.lm<-function(f1,data=NULL,iter=999){
  form.in<-formula(f1)
  Terms<-terms(form.in,keep.order=TRUE)
  Y<-eval(form.in[[2]],parent.frame())
  if (length(dim(Y))!=2){
    stop("Response matrix (shape) not a 2D array. Use 'two.d.array' first.")  }
  if(any(is.na(Y))==T){
    stop("Response data matrix (shape) contains missing values. Estimate these first (see 'estimate.missing').")  }
  if(is.null(dimnames(Y)[[1]])){
    print("No specimen names in response matrix. Assuming specimens in same order.")  }
  df<-df.tmp<-SS.tmp<-SS.obs<-F<-Rsq<-array()
  dat<-model.frame(form.in,data)
  for (i in 1:ncol(attr(Terms, "factors"))){
    mod.mat<-model.matrix(Terms[1:i],data=dat)
    SS.tmp[i]<-sum(dist(predict(lm(Y~mod.mat)))^2)/(nrow(Y))
    ifelse(i==1, SS.obs[i]<-SS.tmp[i], SS.obs[i]<-(SS.tmp[i]-SS.tmp[i-1]))
    df.tmp[i]<-ifelse(ncol(mod.mat)==1,1,(ncol(mod.mat)-1))
    ifelse(i==1, df[i]<-df.tmp[i], df[i]<-(df.tmp[i]-df.tmp[i-1]))
  }
  MS<-SS.obs/df
  mod.mat<-model.matrix(Terms[1])
  SS.tot <- sum(dist(predict(lm(Y ~ mod.mat)))^2)/(nrow(Y)) + 
    sum(dist(resid(lm(Y ~ mod.mat)))^2)/(nrow(Y))
  SS.res<-SS.tot-sum(SS.obs)
  df.tot<-nrow(Y)-1
  df.res<-nrow(Y)-1-sum(df)
  MS.tot<-SS.tot/df.tot
  MS.res<-SS.res/df.res
  Rsq<-SS.obs/SS.tot
  F<-MS/MS.res
  P.val<-array(1,dim=length(SS.obs))
  for(i in 1:iter){
    SS.tmp<-SS.r<-array()
    Y.r<-Y[sample(nrow(Y)),]  
    for (ii in 1:ncol(attr(Terms, "factors"))){
      mod.mat<-model.matrix(Terms[1:ii])
      SS.tmp[ii]<-sum(dist(predict(lm(Y.r~mod.mat)))^2)/(nrow(Y))
      ifelse(ii==1, SS.r[ii]<-SS.tmp[ii], SS.r[ii]<-(SS.tmp[ii]-SS.tmp[ii-1]))
    }
    P.val<-ifelse(SS.r>=SS.obs, P.val+1,P.val) 
  }
  P.val<-P.val/(iter+1)
  anova.tab<-cbind(df,SS.obs,MS,F,P.val,Rsq)
  anova.tab<-rbind(anova.tab,c(df.tot,SS.tot,MS.tot,NA,NA,NA))
  rownames(anova.tab)<-c(colnames(attr(Terms, "factors")), "Total")
  return(anova.tab)
}