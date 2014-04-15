#' Plot allometric patterns in landmark data
#'
#' Function plots allometry curves for a set of specimens
#'
#' The function performs a regression of shape on size, and generates a plot that describes the 
#' multivariate relationship between size and shape 
#'   derived from landmark data (i.e., allometry). It is assumed that the landmarks have previously been 
#'   aligned using Generalized Procrustes Analysis (GPA) [e.g., with \code{\link{gpagen}}]. The abscissa 
#'   of the plot is log(centroid size) while the ordinate represents shape. Three complementary approaches 
#'   can be implemented to visualize allometry: 
#'  \enumerate{
#'   \item {If "method=CAC" (the default) the function calculates the 
#'   common allometric component of the shape data, which is an estimate of the average allometric trend 
#'   within groups (Mitteroecker et al. 2004). The function also calculates the residual shape component (RSC) for 
#'   the data.}
#'   \item {If "method=RegScore" the function calculates shape scores 
#'   from the regression of shape on size, and plots these versus size (Drake and Klingenberg 2008). 
#'   For a single group, these shape scores are mathematically identical to the CAC (Adams et al. 2013).}
#'   \item {If "method=PredLine" the function calculates predicted values from a regression of shape on size, and 
#'   plots the first principal component of the predicted values versus size as a stylized graphic of the 
#'   allometric trend (Adams and Nistri 2010). }
#'   }
#'   For all methods, both centroid size and allometry scores are returned. Optionally, deformation grids can be 
#'   requested, which display the shape of the smallest and largest specimens relative to the average specimen (using 
#'   'warpgrids=T' or 'warpgrids=F'). 
#'   Finally, if groups are provided, the above approaches are implemented while 
#'   accounting for within-group patterns of covariation (see references for explanation). In this case,
#'   the regression is of the form: shape~size+groups (Note: to examine the interaction term use \code{\link{procD.lm}}).
#'   Specimens from each group are plotted using distinct colors based on the order in which the groups are
#'   found in the dataset, and using R's standard color palette: black, red, green, blue, cyan, magenta,
#'   yellow, and gray. 
#'
#' @param A An array (p x k x n) containing landmark coordinates for a set of specimens 
#' @param sz A vector of centroid size measures for all specimens 
#' @param groups An optional vector containing group labels for each specimen if available 
#' @param method Method for estimating allometric shape components; see below for details
#' @param warpgrids A logical value indicating whether deformation grids for small and large shapes 
#'  should be displayed
#'  @param mesh A mesh3d object to be warped to represent shape deformation of the directional and fluctuating components
#' of asymmetry if {warpgrids= TRUE} (see \code{\link{warpRefMesh}}).
#' @param iter Number of iterations for significance testing
#' @param label A logical value indicating whether labels for each specimen should be displayed
#' @param verbose A logical value indicating whether the output is basic or verbose (see Value below)
#' @keywords analysis
#' @keywords visualization
#' @export
#' @return Function returns an ANOVA table of statistical results for log centroid size: df, SS, MS, Prand.
#' If verbose=TRUE, function returns a list with the following components:
#'  \item{ProcDist.lm}{An ANOVA table as above}
#'  \item{allom.score}{ A matrix of the allometry shape scores}
#'  \item{logCsize}{ A matrix of log centroid size}
#'  \item{pred.shape}{A matrix containing the predicted shapes from the regression}
#'  \item{resid.shape}{ The residual shape component (RSC) of the data ("method=CAC" only)}
#' @author Dean Adams
#' @references Adams, D.C., F.J. Rohlf, and D.E. Slice. 2013. A field comes of age: geometric morphometrics 
#'   in the 21st century. Hystrix. 24:7-14. 
#' @references Adams, D. C., and A. Nistri. 2010. Ontogenetic convergence and evolution of foot morphology 
#'   in European cave salamanders (Family: Plethodontidae). BMC Evol. Biol. 10:1-10.
#' @references Drake, A. G., and C. P. Klingenberg. 2008. The pace of morphological change: Historical 
#'   transformation of skull shape in St Bernard dogs. Proceedings of the Royal Society B, Biological Sciences 275:71'76.
#' @references Mitteroecker, P., P. Gunz, M. Bernhard, K. Schaefer, and F. L. Bookstein. 2004. 
#'   Comparison of cranial ontogenetic trajectories among great apes and humans. J. Hum. Evol. 46:679-698.
#' @examples
#' data(ratland) 
#' Y.gpa<-gpagen(ratland)    #GPA-alignment
#' 
#' #Using CAC for plot
#' plotAllometry(Y.gpa$coords,Y.gpa$Csize,method="CAC", iter=5)
#'
#' #Using Regression Scores for plot
#' plotAllometry(Y.gpa$coords,Y.gpa$Csize,method="RegScore", iter=5)
#'
#' #Using predicted allometry curve for plot
#' plotAllometry(Y.gpa$coords,Y.gpa$Csize,method="PredLine", iter=5)
plotAllometry<-function(A,sz,groups=NULL,method=c("CAC","RegScore","PredLine"),warpgrids=TRUE,
                        iter=99,label=FALSE, mesh=NULL, verbose=FALSE){
  method <- match.arg(method)
  if (length(dim(A))!=3){
    stop("Data matrix 1 not a 3D array (see 'arrayspecs').")  }
  if(any(is.na(A))==T){
    stop("Data matrix contains missing values. Estimate these first (see 'estimate.missing').")  }
  if(is.null(dimnames(A)[[3]])){
    print("No specimen names in data matrix. Assuming specimens in same order.")  }
  csz<-as.matrix(log(sz))
  n<-nrow(csz)
  if(is.null(rownames(csz))){
    print("No specimen names in size vector. Assuming specimens in same order.")  }
  y<-two.d.array(A)
  if(nrow(y)!=nrow(csz)){
    stop("Number of specimens differs from number of values in size vector.")  }
  if(is.null(rownames(y))==FALSE && is.null(rownames(csz))==FALSE){
    mtch<-y[is.na( match(rownames(y),rownames(csz)))]
    if (length(mtch)>0){stop("Specimen names in data set don't match those in size vector.")  }
  }
  if(is.null(rownames(y))==FALSE && is.null(rownames(csz))==FALSE){
    csz<-csz[rownames(y),]
  }
  if(!is.null(groups)){
    groups<-as.factor(groups)    
    if(is.null(names(groups))){
      print("No specimen names in grouping variable. Assuming specimens in same order.")  }
  }
  if(is.null(rownames(y))==FALSE && is.null(names(groups))==FALSE){
    mtch<-y[is.na( match(rownames(y),names(groups)))]
    if (length(mtch)>0){stop("Specimen names in data set don't match those in grouping variable.")  }
  }
  if(is.null(rownames(y))==FALSE && is.null(names(groups))==FALSE){
    groups<-groups[rownames(y)]
  }
  if(is.null(groups)){lm.res<-procD.lm(y~csz,iter=iter)}  
  if(!is.null(groups)){lm.res<-procD.lm(y~csz+groups,iter=iter)}
  if(is.null(groups)){
    y.mn<-predict(lm(y~1))
    B<-coef(lm(y~csz))
    yhat<-predict(lm(y~csz))
  }
  if(!is.null(groups)){
    y.mn<-predict(lm(y~groups))
    B<-coef(lm(y~csz+groups))
    yhat<-predict(lm(y~csz*groups))
  }
  y.cent<-y-y.mn
  a<-(t(y.cent)%*%csz)%*%(1/(t(csz)%*%csz)); a<-a%*%(1/sqrt(t(a)%*%a))
  CAC<-y.cent%*%a  
    resid<-y.cent%*%(diag(dim(y.cent)[2])-a%*%t(a))
  RSC<-prcomp(resid)$x
  Reg.proj<-y%*%B[2,]%*%sqrt(solve(t(B[2,])%*%B[2,])) 
  pred.val<-prcomp(yhat)$x[,1] 
  Ahat<-arrayspecs(yhat,dim(A)[1],dim(A)[2])
  ref<-mshape(A)
  if(method!="CAC"){
    layout(matrix(c(2,1,1,1,1,1,1,1,3),3,3))   
    if(method=="RegScore"){
      plot(csz,Reg.proj,xlab="log(CSize)", ylab="Shape (Regression Score)",pch=21,bg="black",cex=1.25)
      if(!is.null(groups)){points(csz,Reg.proj,pch=21,bg=groups,cex=1.25)}
      if(label ==T){text(csz,Reg.proj,seq(1,n),adj=c(-.7,-.7))}
      if(warpgrids==T && dim(A)[2]==2){
        arrows(min(csz), (0.7*max(Reg.proj)), min(csz), 0, length = 0.1,lwd = 2)
        arrows(max(csz), (0.7 * min(Reg.proj)), max(csz), 0, length = 0.1,lwd = 2)
      }
    } 
    if(method=="PredLine"){
      plot(csz,pred.val,xlab="log(CSize)", ylab="Shape (Predicted)",pch=21,bg="black",cex=1.25)
      if(!is.null(groups)){points(csz,pred.val,pch=21,bg=groups,cex=1.25)}
      if(label ==T){text(csz,pred.val,seq(1,n),adj=c(-.7,-.7))}
      if(warpgrids==T && dim(A)[2]==2){
        arrows(min(csz), (0.7*max(pred.val)), min(csz), 0, length = 0.1,lwd = 2)
        arrows(max(csz), (0.7 * min(pred.val)), max(csz), 0, length = 0.1,lwd = 2)
      }
    }
    if(warpgrids==T && dim(A)[2]==2){
      tps(ref,Ahat[,,which.min(csz)],20)
      tps(ref,Ahat[,,which.max(csz)],20)
    }
    layout(1)    
  }
  if(method=="CAC"){
    layout(matrix(c(3,1,1,1,1,1,1,1,4,2,2,2,2,2,2,2,2,2),3,6))   
    plot(csz,CAC,xlab="log(CSize)", ylab="CAC",pch=21,bg="black",cex=1.25)
    if(warpgrids==T && dim(A)[2]==2){
      arrows(min(csz), (0.7*max(CAC)), min(csz), 0, length = 0.1,lwd = 2)
      arrows(max(csz), (0.7 * min(CAC)), max(csz), 0, length = 0.1,lwd = 2)
    }
    if(!is.null(groups)){points(csz,CAC,pch=21,bg=groups,cex=1.25)}
    if(label ==T){text(csz,CAC,seq(1,n),adj=c(-.7,-.7))}
    plot(CAC,RSC[,1], xlab="CAC",ylab="RSC 1", pch=21,bg="black",cex=1.25)
    if(!is.null(groups)){points(CAC,RSC,pch=21,bg=groups,cex=1.25)}
    if(label ==T){text(CAC,RSC,seq(1,n),adj=c(-.7,-.7))}
    if(warpgrids==T && dim(A)[2]==2){
      tps(ref,Ahat[,,which.min(csz)],20)
      tps(ref,Ahat[,,which.max(csz)],20)
    }
    layout(1)
  }
  if(warpgrids==T && dim(A)[2]==3){
    if (is.null(mesh)==TRUE){
      open3d()
      plot3d(Ahat[,,which.min(csz)],type="s",col="gray",main="Shape at minimum size",size=1.25,aspect=FALSE)
      open3d()
      plot3d(Ahat[,,which.max(csz)],type="s",col="gray",main="Shape at maximum size",size=1.25,aspect=FALSE)
    }
    if(is.null(mesh)==FALSE){
      plotRefToTarget(ref, Ahat[,,which.min(csz)], mesh, method = "surface")
      title3d(main="Shape at minimum size")
      plotRefToTarget(ref, Ahat[,,which.max(csz)], mesh, method = "surface")
      title3d(main="Shape at maximum size")
    }
  }
  if(verbose==TRUE){ 
    if(method=="CAC"){return(list(allom.score=CAC,resid.shape=RSC,logCsize=csz,ProcDist.lm=lm.res,pred.shape=Ahat))}
    if(method=="RegScore"){return(list(allom.score=Reg.proj,logCsize=csz,ProcDist.lm=lm.res,pred.shape=Ahat))}
    if(method=="PredLine"){return(list(allom.score=pred.val,logCsize=csz,ProcDist.lm=lm.res,pred.shape=Ahat))}
  }
  if(verbose==FALSE){ return(list(ProcDist.lm=lm.res))}
}