#' Edit 3D template
#'
#' An interactive function to remove landmarks from a 3D template file. 
#' 
#' Function edits a 'template.txt' file made by \code{\link{buildtemplate}}, which must be in current working directory.
#'  Function overwrites 'template.txt' in working directory with edited version. Use read.table("template.txt", header = T).
#' \subsection{Selection}{ 
#' Choosing which landmarks will be deleted involves landmark selection using a mouse in the rgl plot window. 
#' With a standard 3-button (PC) buildtemplate uses:
#' \enumerate{
#'  \item the RIGHT mouse button (primary) to choose points to be deleted (click-drag a box around landmark),
#'  \item the LEFT mouse button (secondary) is used to rotate mesh, 
#'  \item the mouse SCROLLER (third/middle) is used to zoom in and out.
#' }
#' NOTE: Digitizing functions on MACINTOSH computers using a standard 3-button mice works as specified. Macs using platform 
#' specific single button mice, XQuartz must be configured: go to Preferences > Input > tick "Emulate three button mouse":
#' \enumerate{
#'  \item press button to rotate 3D mesh,
#'  \item press button while pressing COMMAND key to select points to be digitized (click-drag a box around a vertex to select as landmark),
#'  \item press button while pressing OPTION key to adjust mesh perspective.
#'  \item the mouse SCROLLER or trackpad two finger scroll is used to zoom in an out.
#'  }
#' }
#' @param template Matrix of template 3D coordinates.
#' @param fixed Number of "fixed" landmark points (non surface sliding points)
#' @param n Number of points to be removed 
#' @return Function returns a matrix containing the x,y,z coordinates of the new template landmarks.
#' Function also writes to working directory 'template.txt' containing the x,y,z coordinates of the template 
#' @export
#' @keywords digitizing
#' @author Erik Otarola-Castillo & Emma Sherratt
editTemplate<-function(template, fixed, n){
  if (is.null(dim(template))) stop ("File is not a matrix of 3D coordinates.")
  if (dim(template)[2]!=3) stop ("File is not a matrix of 3D coordinates.")
  spec.name<-deparse(substitute(template))
  clear3d();plot3d(template[-(1:fixed[1]),],size=7,col="blue",aspect=FALSE)
  points3d(template[(1:fixed),],size=10,color="red",add=TRUE)  
  cat("Remove Template Points","\n")
  selected2<-NULL
  for (i in 1:n)        {
    selected2.temp<-NULL
    f<-NULL
    keep<-NULL
    f<-select3d(button="right")
    cat( i, "of", n, "points have been removed" , "\n")
    keep<-f(template[,1],template[,2],template[,3])
    selected2.temp<-which(keep==TRUE)
    selected2<-c(selected2,selected2.temp)
    points3d(template[selected2.temp,][1],template[selected2.temp,][2],template[selected2.temp,][3],size=12,color="dark green",add=TRUE)
    points3d(template[-(1:fixed[1]),],size=7,col="blue",add=TRUE)
    points3d(template[(1:fixed),],size=10,color="red",add=TRUE)
  }  
  template<-cbind(template[-(selected2),][,1],template[-(selected2),][,2],template[-(selected2),][,3])
  write.table(template,file="template.txt",col.names=TRUE)
  clear3d();plot3d(template[-(1:fixed[1]),],size=7,col="blue",aspect=FALSE)
  points3d(template[(1:fixed),],size=10,color="red",add=TRUE)
  return(template)
}