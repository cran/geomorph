#' Plot shape differences between a reference and target specimen
#'
#' Function plots shape differences between a reference and target specimen
#'
#' The function generates a plot of the shape differences of a target specimen relative to a reference 
#'  specimen. The option {mag} allows the user to indicates the degree of magnification to be used when 
#'  displaying the shape difference. The function will plot either two- or three-dimensional data. 
#'  Four methods for plots are available:
#'  \enumerate{
#'  \item {TPS} a thin-plate spline deformation grid is generated. For 3D data, 
#'  this method will generate thin-plate spline deformations in the x-y and x-z planes.
#'  \item {vector}: a plot showing the vector displacements between corresponding landmarks in the reference 
#'  and target specimen is shown. 
#'  \item {points} a plot is displayed with the landmarks in the target (black) 
#'  overlaying those of the reference (gray). Additionally, if a matrix of links is provided, the 
#'  landmarks of the mean shape will be connected by lines.  The link matrix is an M x 2 matrix, where 
#'  M is the desired number of links. Each row of the link matrix designates the two landmarks to be 
#'  connected by that link. 
#'  \item {surface} a mesh3d surface is warped using thin-plate spline (for 3D data only). 
#'  Requires mesh3d object in option {mesh}, made using \code{\link{warpRefMesh}}. 
#'  }
#'  This function combines numerous plotting functions found in Claude (2008).
#'
#' @param M1 Matrix of landmark coordinates for the first (reference) specimen
#' @param M2 Matrix of landmark coordinates for the second (target) specimen
#' @param mesh A mesh3d object for use with {method="surface"}
#' @param method Method used to visualize shape difference; see below for details
#' @param mag The desired magnification to be used when visualizing the shape difference (e.g., mag=2)
#' @param links An optional matrix defining for links between landmarks
#' @param ... Additional parameters to be passed to \code{\link{plot}}, \code{\link{plot3d}} or \code{\link{shade3d}}.
#' @return If using {method="surface"}, function will return the warped mesh3d object.
#' @keywords visualization
#' @export
#' @author Dean Adams & Emma Sherratt
#' @references Claude, J. 2008. Morphometrics with R. Springer, New York. 
#' @examples
#' data(plethodon) 
#' Y.gpa<-gpagen(plethodon$land)    #GPA-alignment
#' ref<-mshape(Y.gpa$coords)
#' 
#' # Different plotting options
#' plotRefToTarget(ref,Y.gpa$coords[,,39])
#'
#' plotRefToTarget(ref,Y.gpa$coords[,,39],mag=3)   #magnify difference by 3X
#'
#' plotRefToTarget(ref,Y.gpa$coords[,,39],method="vector")
#'
#' plotRefToTarget(ref,Y.gpa$coords[,,39],method="points")
#'
#' # Three dimensional data
#' 
#' data(scallops)
#' Y.gpa<-gpagen(A=scallops$coorddata, curves=scallops$curvslide, surfaces=scallops$surfslide)
#' ref<-mshape(Y.gpa$coords)
#' plotRefToTarget(ref,Y.gpa$coords[,,1],method="points")
plotRefToTarget<-function(M1,M2,mesh= NULL,method=c("TPS","vector","points","surface"),
                          mag=1.0,links=NULL,...){
  method <- match.arg(method)
  if(any(is.na(M1))==T){
    stop("Data contains missing values. Estimate these first (see 'estimate.missing').")  }
  if(any(is.na(M2))==T){
    stop("Data contains missing values. Estimate these first (see 'estimate.missing').")  }
  k<-dim(M1)[2]
  mag<-(mag-1)
  M2<-M2+(M2-M1)*mag
  limits = function(x,s){ 
    r = range(x)
    rc=scale(r,scale=F)
    l=mean(r)+s*rc
  }
  if(k==2){
    if(method=="TPS"){
      tps(M1,M2,20)
      if(is.null(links)==FALSE){
        for (i in 1:nrow(links)){
          segments(M2[links[i,1],1],M2[links[i,1],2],M2[links[i,2],1],M2[links[i,2],2],lwd=2)
        }
      }
    }
    if(method=="vector"){
      plot(M1,asp=1,pch=21,bg="gray",xlim=limits(M1[,1],1.25),ylim=limits(M1[,2],1.25),cex=1,xlab="x",ylab="y",...)
      arrows(M1[,1],M1[,2],M2[,1],M2[,2],length=0.075,lwd=2)
      if(is.null(links)==FALSE){
        for (i in 1:nrow(links)){
          segments(M2[links[i,1],1],M2[links[i,1],2],M2[links[i,2],1],M2[links[i,2],2],lwd=1,lty=2)
        }
      }
    }
    if(method=="points"){
      plot(M1,asp=1,pch=21,bg="gray",xlim=limits(M1[,1],1.25),ylim=limits(M1[,2],1.25),cex=1,xlab="x",ylab="y",...)
      points(M2,asp=1,pch=21,bg="black",cex=1)
      if(is.null(links)==FALSE){
        for (i in 1:nrow(links)){
          segments(M2[links[i,1],1],M2[links[i,1],2],M2[links[i,2],1],M2[links[i,2],2],lwd=1,lty=2)
        }
      }
    }
    if(method=="surface"){
      stop("Surface plotting for 3D landmarks only.")
    }      
  }
  if(k==3){
    if(method=="TPS"){
      layout(matrix(c(1,2),1,2))
      par(mar=c(1,1,1,1))
      tps(M1[,1:2],M2[,1:2],20)
      if(is.null(links)==FALSE){
        for (i in 1:nrow(links)){
          segments(M2[links[i,1],1],M2[links[i,1],2],M2[links[i,2],1],M2[links[i,2],2],lwd=1)
        }
      }
      title("X,Y tps grid")
      b<-c(1,3)
      tps(M1[,b],M2[,b],20)
      if(is.null(links)==FALSE){
        for (i in 1:nrow(links)){
          segments(M2[links[i,1],1],M2[links[i,1],3],M2[links[i,2],1],M2[links[i,2],3],lwd=1)
        }
      }
      title("Y,Z tps grid")
      layout(1)
    }
    if(method=="vector"){
      plot3d(M1,type="s",col="gray",size=1.25,aspect=FALSE,...)
      for (i in 1:nrow(M1)){
        segments3d(rbind(M1[i,],M2[i,]),lwd=2)
      }
      if(is.null(links)==FALSE){
        for (i in 1:nrow(links)){
          segments3d(rbind(M2[links[i,1],],M2[links[i,2],]),lwd=.75)
        }
      }
    }
    if(method=="points"){
      plot3d(M1,type="s",col="gray",size=1.25,aspect=FALSE,...)
      points3d(M2,color="black",size=5)
      if(is.null(links)==FALSE){
        for (i in 1:nrow(links)){
          segments3d(rbind(M2[links[i,1],],M2[links[i,2],]),lwd=.75)
        }
      }
    }
    if (method == "surface") {
      if(is.null(mesh)==T){
        stop("Surface plotting requires a template mesh3d object (see 'warpRefMesh').")
      }
      warp.PLY <- mesh
      vb <- as.matrix(t(mesh$vb)[,-4])
      warp <- tps2d3d(vb, M1, M2)
      warp.PLY$vb <- rbind(t(warp), 1)
      open3d(); shade3d(warp.PLY, ...)
      return(warp.PLY)
    }
  }
}