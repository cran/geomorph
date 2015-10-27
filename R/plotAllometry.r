#' Plot allometric patterns in landmark data
#'
#' Function plots allometry curves for a set of specimens
#'
#' The function performs a regression of shape on size, and generates a plot that describes the 
#' multivariate relationship between size and shape 
#'   derived from landmark data (i.e., allometry). It is assumed that the landmarks have previously been 
#'   aligned using Generalized Procrustes Analysis (GPA) [e.g., with \code{\link{gpagen}}]. The abscissa 
#'   of the plot is size (or log of size) while the ordinate represents shape [NOTE: the function takes the 
#'   input size and performed log-transformation automatically by default (logsz = TRUE), as log(Size) should be used]. 
#'   Three complementary approaches can be implemented to visualize allometry: 
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
#'   'warpgrids=TRUE' or 'warpgrids=FALSE'). 
#'   Finally, if groups are provided, the above approaches are implemented while 
#'   accounting for within-group patterns of covariation (see references for explanation). In this case,
#'   the regression is of the form: shape~size+groups (Note: to examine the interaction term use \code{\link{procD.lm}}).
#'   Specimens from each group are plotted using distinct colors based on the order in which the groups are
#'   found in the dataset, and using R's standard color palette: black, red, green, blue, cyan, magenta,
#'   yellow, and gray. 
#'
#' @param f1 A formula (of the form shape ~ size); shape can be a matrix (n x [p1 x k]) or 3D array (p1 x k x n) containing 
#' GPA-aligned coordinates for the specimens
#' @param f2 An optional right hand formula for groups (e.g., ~ groups); must be a single factor
#' @param method Method for estimating allometric shape components; see below for details
#' @param warpgrids A logical value indicating whether deformation grids for small and large shapes 
#'  should be displayed (note: if groups are provided no TPS grids are shown)
#' @param iter Number of iterations for significance testing
#' @param label An optional vector indicating labels for each specimen that are to be displayed
#' @param mesh A mesh3d object to be warped to represent shape deformation of the minimum and maximum size if {warpgrids=TRUE} (see \code{\link{warpRefMesh}}).
#' @param logsz A logical value indicating whether log(size) is used 
#' @param RRPP A logical value to indicate if a randomized residual permutation procedure (RRPP) should be used for statistical tests
#' @param verbose A logical value indicating whether the output is basic or verbose (see Value below)
#' @keywords analysis
#' @keywords visualization
#' @export
#' @return Function returns an ANOVA table of statistical results for size: df, SS, MS, F, Z,Prand.
#' If verbose=TRUE, function returns a list with the following components:
#'  \item{ProcDist.lm}{An ANOVA table as above}
#'  \item{allom.score}{ A matrix of the allometry shape scores}
#'  \item{logSize}{ A matrix of log size}
#'  \item{pred.shape}{A matrix containing the predicted shapes from the regression}
#'  \item{resid.shape}{ The residual shape component (RSC) of the data ("method=CAC" only)}
#' @author Dean Adams and Michael Collyer
#' @references Adams, D.C., F.J. Rohlf, and D.E. Slice. 2013. A field comes of age: geometric morphometrics 
#'   in the 21st century. Hystrix. 24:7-14. 
#' @references Adams, D. C., and A. Nistri. 2010. Ontogenetic convergence and evolution of foot morphology 
#'   in European cave salamanders (Family: Plethodontidae). BMC Evol. Biol. 10:1-10.
#' @references Drake, A. G., and C. P. Klingenberg. 2008. The pace of morphological change: Historical 
#'   transformation of skull shape in St Bernard dogs. Proc. R. Soc. B. 275:71-76.
#' @references Mitteroecker, P., P. Gunz, M. Bernhard, K. Schaefer, and F. L. Bookstein. 2004. 
#'   Comparison of cranial ontogenetic trajectories among great apes and humans. J. Hum. Evol. 46:679-698.
#' @examples
#' data(ratland) 
#' Y.gpa<-gpagen(ratland,PrinAxes=FALSE)    #GPA-alignment
#' 
#' #Using CAC for plot
#' plotAllometry(Y.gpa$coords ~ Y.gpa$Csize, method="CAC", iter=9)
#'
#' #Using Regression Scores for plot
#' plotAllometry(Y.gpa$coords ~ Y.gpa$Csize, method="RegScore", iter=9)
#'
#' #Using predicted allometry curve for plot
#' plotAllometry(Y.gpa$coords ~ Y.gpa$Csize, method="PredLine", iter=9)
plotAllometry<-function(f1, f2 = NULL, method=c("CAC","RegScore","PredLine"),warpgrids=TRUE,
                        iter=999,label=NULL, mesh=NULL, logsz = TRUE, RRPP=FALSE, verbose=FALSE){
  form.in <- as.formula(f1)
  A <- eval(form.in[[2]], parent.frame())
  size.df <- allometry.data.frame(form.in)
  Y <- as.matrix(size.df$Y)
  Size <- size.df$Size
  n<-nrow(Y)
  if(logsz == TRUE) {
    xlab <- "log(Size)"
    fnew <- as.formula("Y~log(Size)")
    print(noquote("Natural log of size is used."))
    } else 
      { 
        xlab <- "Size"
        fnew <- as.formula("Y~Size")
        print(noquote("Size has not been log transformed."))
      }                                                                                          
  pf1<- procD.fit(fnew, data=size.df, keep.order=FALSE)
  method <- match.arg(method)

  if(!is.null(f2)) {
    fac.mf <- model.frame(f2)
    if(dim(fac.mf)[[2]] > 1) stop("Only a single grouping variable can be used")
    fTerms <- terms(f2)
    cTerms <- terms(fnew, data = size.df)
    all.terms <- c(attr(cTerms, "term.labels"), attr(fTerms, "term.labels"))
    f3 <- as.formula(paste(" Y ~",paste(all.terms, collapse="+")))
    f4 <- as.formula(paste(" Y ~",paste(all.terms, collapse="*")))
    pf4 <-procD.fit(f4, data=size.df)
    pf3 <-procD.fit(f3, data=size.df)
    Xs2 <- pf4$Xs
    Xs1 <- pf3$Xs
    k <- length(pf4$Terms) 
    X2 <- pf4$Xs[[k+1]]
    X1 <- pf3$Xs[[k]]
    anova.parts.obs2 <- anova.parts(pf4, X=Xs2)
    anova.parts.obs1 <- anova.parts(pf3, X=Xs1)
    anova.tab2 <-anova.parts.obs2$table 
    anova.tab <-anova.parts.obs1$table 
    P <- array(0, c(k, 1, iter+1))
    if(RRPP == TRUE) P <- SS.random(pf4,Yalt="RRPP", iter=iter) else P <- SS.random(pf4, Yalt="resample", iter=iter)
    P.val <- Pval.matrix(P)
    Z <- Effect.size.matrix(P)
    anova.tab2 <- data.frame(anova.tab2, Z = c(Z, NA, NA), P.value = c(P.val, NA, NA))
    anova.tab <- data.frame(anova.tab, Z = c(Z[1:(k-1)], NA, NA), P.value = c(P.val[1:(k-1)], NA, NA))
      
  } else {
    fac.mf <- NULL 
    Xs <- pf1$Xs
    k <- length(pf1$Terms) 
    X <- as.matrix(Xs[[k+1]])
    anova.parts.obs <- anova.parts(pf1, X=Xs, Yalt="observed")
    anova.tab <-anova.parts.obs$table
    if(RRPP == TRUE) P <- SS.random(pf1,Yalt="RRPP", iter=iter) else P <- SS.random(pf1, Yalt="resample", iter=iter)
    P.val <- Pval.matrix(P)
    Z <- Effect.size.matrix(P)
    anova.tab <- data.frame(anova.tab, Z = c(Z, NA, NA), P.value = c(P.val, NA, NA))
  }
  
  if(is.null(f2)){
    y.mn<-predict(lm(Y~1))
    B<-as.matrix(pf1$fit)
    yhat<-X%*%B
  }
  if(!is.null(f2)){
    y.mn<-predict(lm(Y~model.matrix(fac.mf) - 1), data = data.frame(model.frame(f2)))
    B<-coef(lm(Y~X1 -1))
    yhat<-predict(lm(Y~X2 - 1))
    if(anova.tab2[3,7]>0.05){
      yhat<-predict(lm(Y~X1 - 1))      
    }
  } 
  asp = NULL
  if(anova.tab[1,7]>0.05){ asp <- 1}
  y.cent<-Y-y.mn
  size<-pf1$mf[[2]]
  a<-(t(y.cent)%*%size)%*%(1/(t(size)%*%size)); a<-a%*%(1/sqrt(t(a)%*%a))
  CAC<-y.cent%*%a  
  resid<-y.cent%*%(diag(dim(y.cent)[2])-a%*%t(a))
  RSC<-prcomp(resid)$x
  Reg.proj<-Y%*%B[2,]%*%sqrt(solve(t(B[2,])%*%B[2,])) 
  pred.val<-prcomp(yhat)$x[,1] 
  Ahat<-arrayspecs(yhat,dim(A)[1],dim(A)[2])
  ref<-mshape(A)
  pt.cols <- as.numeric(fac.mf[[1]])
  if(method!="CAC"){
    layout(matrix(c(2,1,1,1,1,1,1,1,3),3,3))   
    if(method=="RegScore"){
      plot(size,Reg.proj,xlab=xlab, ylab="Shape (Regression Score)",pch=21,bg="black",cex=1.25, asp=asp)
      if(!is.null(f2)){points(size,Reg.proj,pch=21,bg=pt.cols,cex=1.25)}
      if (length(label!=0)) {
        if(isTRUE(label)){text(size,Reg.proj,seq(1, n),adj=c(-0.7,-0.7)) }
        else{text(size,Reg.proj,label,adj=c(-0.1,-0.1))}
      }
      if(is.null(f2)){
        if(warpgrids==T && dim(A)[2]==2){
          arrows(min(size), (0.7*max(Reg.proj)), min(size), 0, length = 0.1,lwd = 2)
          arrows(max(size), (0.7 * min(Reg.proj)), max(size), 0, length = 0.1,lwd = 2)
        }
      }
    } 
    if(method=="PredLine"){
      plot(size,pred.val,xlab=xlab, ylab="Shape (Predicted)",pch=21,bg="black",cex=1.25, asp=asp)
      if(!is.null(f2)){points(size,pred.val,pch=21,bg=pt.cols,cex=1.25)}
      if (length(label!=0)) {
        if(isTRUE(label)){text(size,pred.val,seq(1, n),adj=c(-0.7,-0.7)) }
        else{text(size,pred.val,label,adj=c(-0.1,-0.1))}
      }
      if(is.null(f2)){
        if(warpgrids==T && dim(A)[2]==2){
          arrows(min(size), (0.7*max(pred.val)), min(size), 0, length = 0.1,lwd = 2)
          arrows(max(size), (0.7 * min(pred.val)), max(size), 0, length = 0.1,lwd = 2)
        }
      }
    }
    if(is.null(f2)){
      if(warpgrids==T && dim(A)[2]==2){
        tps(ref,Ahat[,,which.min(size)],20)
        tps(ref,Ahat[,,which.max(size)],20)
      }
    }
    layout(1)    
  }
  if(method=="CAC"){
    layout(matrix(c(3,1,1,1,1,1,1,1,4,2,2,2,2,2,2,2,2,2),3,6))   
    plot(size,CAC,xlab=xlab, ylab="CAC",pch=21,bg="black",cex=1.25, asp=asp)
    if(is.null(f2)){
      if(warpgrids==T && dim(A)[2]==2){
        arrows(min(size), (0.7*max(CAC)), min(size), 0, length = 0.1,lwd = 2)
        arrows(max(size), (0.7 * min(CAC)), max(size), 0, length = 0.1,lwd = 2)
      }
    }
    if(!is.null(f2)){points(size,CAC,pch=21,bg=pt.cols,cex=1.25)}
    if (length(label!=0)) {
      if(isTRUE(label)){text(size,CAC,seq(1, n),adj=c(-0.7,-0.7)) }
      else{text(size,CAC,label,adj=c(-0.1,-0.1))}
    }
    plot(CAC,RSC[,1], xlab="CAC",ylab="RSC 1", pch=21,bg="black",cex=1.25, asp=asp)
    if(!is.null(f2)){points(CAC,RSC[,1],pch=21,bg=pt.cols,cex=1.25)}
    if (length(label!=0)) {
      if(isTRUE(label)){text(CAC,RSC,seq(1, n),adj=c(-0.7,-0.7)) }
      else{text(CAC,RSC,label,adj=c(-0.1,-0.1))}
    }
    if(is.null(f2)){
      if(warpgrids==T && dim(A)[2]==2){
        tps(ref,Ahat[,,which.min(size)],20)
        tps(ref,Ahat[,,which.max(size)],20)
      }
    }
    layout(1)
  }
  if(warpgrids==T && dim(A)[2]==3){
    if (is.null(mesh)==TRUE){
      open3d()
      plot3d(Ahat[,,which.min(size)],type="s",col="gray",main="Shape at minimum size",size=1.25,aspect=FALSE)
      open3d()
      plot3d(Ahat[,,which.max(size)],type="s",col="gray",main="Shape at maximum size",size=1.25,aspect=FALSE)
    }
    if(is.null(mesh)==FALSE){
      plotRefToTarget(ref, Ahat[,,which.min(size)], mesh, method = "surface")
      title3d(main="Shape at minimum size")
      plotRefToTarget(ref, Ahat[,,which.max(size)], mesh, method = "surface")
      title3d(main="Shape at maximum size")
    }
  }
  
  if(RRPP == TRUE) {anova.title = "\nRandomized Residual Permutation Procedure used\n"
  } else {anova.title = "\nRandomization of Raw Values used\n"}
  attr(anova.tab, "heading") <- paste("\nType I (Sequential) Sums of Squares and Cross-products\n",anova.title)
  class(anova.tab) <- c("anova", class(anova.tab))
  if(!is.null(f2)){
    attr(anova.tab2, "heading") <- paste("\nType I (Sequential) Sums of Squares and Cross-products\n",anova.title)
    class(anova.tab2) <- c("anova", class(anova.tab))
  }
  
  if(verbose==TRUE){ 
    if(method=="CAC") return(list(allom.score=CAC,resid.shape=RSC,logSize=log(Size),ProcDist.lm=anova.tab,pred.shape=Ahat))
    if(method=="RegScore") return(list(allom.score=Reg.proj,logSize=log(Size),ProcDist.lm=anova.tab,pred.shape=Ahat))
    if(method=="PredLine") return(list(allom.score=pred.val,logSize=log(Size),ProcDist.lm=anova.tab2,pred.shape=Ahat))
  }
  if(verbose==FALSE) return(ProcDist.lm=anova.tab)
}
