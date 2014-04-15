#' Two-block partial least squares analysis for shape data
#'
#' Function performs two-block partial least squares analysis to assess the degree of association between 
#' to blocks of Procrustes-aligned coordinates (or other variables)
#'
#' The function quantifies the degree of association between two blocks of shape data as 
#'   defined by landmark coordinates using partial least squares (see Rohlf and Corti 2000). It is assumed 
#'   that the landmarks have previously been aligned using 
#'   Generalized Procrustes Analysis (GPA) [e.g., with \code{\link{gpagen}}]. 
#'   
#'   A plot of PLS scores from Block1 versus Block2 is provided for the first set of PLS axes. Thin-plate spline 
#'   deformation grids along these axes are also shown (if data were input as a 3D array).
#'   
#' @param A1 A matrix (n x [p x k]) or 3D array (p x k x n) containing GPA-aligned coordinates for the first block
#' @param A2 A matrix (n x [p x k]) or 3D array (p x k x n) containing GPA-aligned coordinates for the second block 
#' @param iter Number of iterations for significance testing
#' @param warpgrids A logical value indicating whether deformation grids for shapes along PC1 should be displayed
#'  (only relevant if data for A1 or A2 [or both] were input as 3D array)
#' @param verbose A logical value indicating whether the output is basic or verbose (see Value below)
#' @export
#' @keywords analysis
#' @author Dean Adams
#' @return Function returns a list with the following components: 
#'   \item{value}{The estimate of association between block}
#'   \item{pvalue}{The significance level of the observed association}
#'   \item{Xscores}{PLS scores for the first block of landmarks (when {verbose=TRUE})}
#'   \item{Yscores}{PLS scores for the second block of landmarks (when {verbose=TRUE})}
#' @references  Rohlf, F.J., and M. Corti. 2000. The use of partial least-squares to study covariation in shape. 
#' Systematic Biology 49: 740-753.
#' @examples
#' data(plethShapeFood) 
#' Y.gpa<-gpagen(plethShapeFood$land)    #GPA-alignment    
#'
#' #2B-PLS between head shape and food use data
#' two.b.pls(Y.gpa$coords,plethShapeFood$food,iter=99)
#' 
two.b.pls<-function(A1,A2,warpgrids = TRUE,iter=999, verbose=FALSE){
  if(any(is.na(A1))==T){
    stop("Data matrix 1 contains missing values. Estimate these first (see 'estimate.missing').")  }
  if(any(is.na(A2))==T){
    stop("Data matrix 2 contains missing values. Estimate these first (see 'estimate.missing').")  }
  if (length(dim(A1))==3){ 
    if(is.null(dimnames(A1)[[3]])){
      print("No specimen names in data matrix 1. Assuming specimens in same order.")  }
    x<-two.d.array(A1)}
  if (length(dim(A1))==2){ 
    if(is.null(dimnames(A1)[[1]])){
      print("No specimen names in data matrix 1. Assuming specimens in same order.")  }
    x<-A1 }
  if (length(dim(A2))==3){ 
    if(is.null(dimnames(A2)[[3]])){
      print("No specimen names in data matrix 2. Assuming specimens in same order.")  }
    y<-two.d.array(A2)}
  if (length(dim(A2))==2){ 
    if(is.null(dimnames(A2)[[1]])){
      print("No specimen names in data matrix 2. Assuming specimens in same order.")  }
    y<-A2 }
  if(nrow(x)!=nrow(y)){
    stop("Data matrices have different numbers of specimens.")  }
  if(is.null(rownames(x))==FALSE && is.null(rownames(y))==FALSE){
    mtch<-x[is.na( match(rownames(x),rownames(y)))]
    if (length(mtch)>0){stop("Specimen names in data sets are not the same.")  }
  }
  if(is.null(rownames(x))==FALSE && is.null(rownames(y))==FALSE){
    y<-y[rownames(x),]
  }
  XY.vcv<-cov(cbind(x,y))
  S12<-XY.vcv[1:dim(x)[2],(dim(x)[2]+1):(dim(x)[2]+dim(y)[2])]; S21<-t(S12)
  S11<-XY.vcv[1:dim(x)[2],1:dim(x)[2]]
  S22<-XY.vcv[(dim(x)[2]+1):(dim(x)[2]+dim(y)[2]),(dim(x)[2]+1):(dim(x)[2]+dim(y)[2])]
  pls<-svd(S12)
  U<-pls$u; V<-pls$v
  XScores<-x%*%U; YScores<-y%*%V
  PLS.obs<-cor(XScores[,1],YScores[,1])
  P.val<-1
  pls.val<-rep(0,iter)
  for(i in 1:iter){
    y.r<-y[sample(nrow(y)),]  
    XY.vcv.r<-cov(cbind(x,y.r))
    S12.r<-XY.vcv.r[1:dim(x)[2],(dim(x)[2]+1):(dim(x)[2]+dim(y.r)[2])]; S21.r<-t(S12.r)
    S11.r<-XY.vcv.r[1:dim(x)[2],1:dim(x)[2]]
    S22.r<-XY.vcv.r[(dim(x)[2]+1):(dim(x)[2]+dim(y.r)[2]),(dim(x)[2]+1):(dim(x)[2]+dim(y.r)[2])]
    pls.r<-svd(S12.r)
    U.r<-pls.r$u; V.r<-pls.r$v
    XScores.r<-x%*%U.r[,1]; YScores.r<-y.r%*%V.r[,1]
    PLS.r<-cor(XScores.r,YScores.r)
    pls.val[i]<-PLS.r
    P.val<-ifelse(PLS.r>=PLS.obs, P.val+1,P.val) 
  }  
  pls.val[iter+1]=PLS.obs
  P.val<-P.val/(iter+1)
  if (length(dim(A1))==2 && length(dim(A2))==2){
    plot(XScores[,1],YScores[,1],pch=21,bg="black",main="PLS Plot",xlab = "PLS1 Block 1",ylab = "PLS1 Block 2")
  }
  if (length(dim(A1))==3){A1.ref<-mshape(A1);
                          pls1.min<-A1[,,which.min(XScores[,1])];pls1.max<-A1[,,which.max(XScores[,1])]}
  if (length(dim(A2))==3){A2.ref<-mshape(A2);
                          pls2.min<-A2[,,which.min(XScores[,1])];pls2.max<-A2[,,which.max(XScores[,1])]}
  if (dim(A1)[2] == 2 ||dim(A2)[2] == 2) {par(mar=c(1,1,1,1)+0.1)
    split.screen(matrix(c(0.22,1,0.22,1,.19,.39,0,.19,.8,1,0,.19,0,.19,.19,.39,0,.19,.8,1), byrow=T, ncol=4))
    screen(1)
    plot(XScores[,1],YScores[,1],pch=21,bg="black",main="PLS1 Plot: Block 1 (X) vs. Block 2 (Y) ",xlab = "PLS1 Block 1",ylab = "PLS1 Block 2")
    if(warpgrids==TRUE){
      if (length(dim(A1))==3  && dim(A1)[2]==2){
        screen(2);       tps(A1.ref, pls1.min, 20,sz=.7)
        screen(3);       tps(A1.ref, pls1.max, 20,sz=.7)     
      }
      if (length(dim(A2))==3  && dim(A2)[2]==2){
        screen(4);       tps(A2.ref, pls2.min, 20,sz=.7)
        screen(5);       tps(A2.ref, pls2.max, 20,sz=.7)
      }
    }
    close.screen(all.screens=TRUE)
    par(mar=c(5.1, 4.1, 4.1, 2.1))
  }
  if (length(dim(A1))==3  && dim(A1)[2] == 3) {
    plot(XScores[,1],YScores[,1],pch=21,bg="black",main="PLS Plot",xlab = "PLS1 Block 1",ylab = "PLS1 Block 2")
    open3d()
    plot3d(pls1.min, type = "s", col = "gray", main = paste("PLS Block1 negative"),size = 1.25, aspect = FALSE)
    open3d()
    plot3d(pls1.max, type = "s", col = "gray", main = paste("PLS Block1 positive"),size = 1.25, aspect = FALSE)
  }  
  if (length(dim(A2))==3  && dim(A2)[2] == 3){
    open3d()
    plot3d(pls2.min, type = "s", col = "gray", main = paste("PLS Block2 negative"),size = 1.25, aspect = FALSE)
    open3d()
    plot3d(pls2.max, type = "s", col = "gray", main = paste("PLS Block2 positive"),size = 1.25, aspect = FALSE)
  }
  if(verbose==TRUE){ return(list(x.scores=XScores[,1],y.scores=YScores[,1],PLS.corr=PLS.obs,pvalue=P.val)) }
  if(verbose==FALSE){ return(list(PLS.corr=PLS.obs,pvalue=P.val))}
}