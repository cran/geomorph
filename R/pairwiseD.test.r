#' Pairwise Group Comparisons
#'
#' Function performs pairwise comparisons among groups using the Euclidean distances among group means.
#'
#' The function performs pairwise comparisons to identify shape among groups. The function is designed as a post-hoc
#'  test to Procrustes ANOVA, where the latter has identified significant shape variation explained by a grouping factor. 
#'  The function takes as input the shape data (y), and a grouping factor (x). It then estimates the Euclidean distances 
#'  among group means, which are used as test values. These are then statistically evaluated through permutation, 
#'  where the rows of the shape matrix are randomized relative to the grouping variable. 
#'
#' The input for the shape data (y) must be in the form of a two-dimensional data matrix of dimension (n x [p x k]), 
#'   rather than a 3D array.  It is assumed that the landmarks have previously 
#'   been aligned using Generalized Procrustes Analysis (GPA) [e.g., with \code{\link{gpagen}}]. The function
#'   \code{\link{two.d.array}} can be used to obtain a two-dimensional data matrix from a 3D array of landmark
#'   coordinates. 
#'
#' @param y A two-dimensional array of shape data
#' @param x A factor defining groups
#' @param iter Number of iterations for permutation test
#' @keywords analysis
#' @export
#' @author Dean Adams
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
#' ### Pairwise comparisons
#' pairwiseD.test(y,plethodon$species,iter=99)
pairwiseD.test<-function(y,x,iter=999){
  if (length(dim(y))!=2){
    stop("Response matrix (shape) not a 2D array. Use 'two.d.array' first.")  } 
  if(any(is.na(y))==T){
    stop("Response data matrix (shape) contains missing values. Estimate these first (see 'estimate.missing').")  }
  if(!is.factor(x)){
    stop("X-variable is not a factor.")}
  yhat.obs<-predict(lm(y~x))
  lsmeans.obs <- rowsum(yhat.obs, x)/as.vector(table(x))     
  D.obs<-as.matrix(dist(lsmeans.obs))   
  PDist<-array(1, dim=c(dim(lsmeans.obs)[1],dim(lsmeans.obs)[1]))
  for(i in 1:iter){
    y.r<-y[sample(nrow(y)),] 
    yhat.rand<-predict(lm(y.r~x))
    lsmeans.rand <- rowsum(yhat.rand, x)/as.vector(table(x))     
    D.rand<-as.matrix(dist(lsmeans.rand))   
    PDist<-ifelse(D.rand>=D.obs, PDist+1, PDist)
  }
  PDist<-PDist/(iter+1)
  return(list(Dist.obs=D.obs,Prob.Dist=PDist))
}