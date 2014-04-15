#' Digitize 2D landmarks.
#'
#' An interactive function to digitize two-dimensional(2D) landmarks from .jpg files.
#'
#' Function for digitizing 2D landmarks on specimen images (.jpg). "nlandmarks" 
#' is the number of landmark points to be digitized by the user. Landmarks should include
#' "true" landmarks and semi-landmarks to be "sliders". For best results, digitizing sequence should proceed 
#' by selecting all true landmark points first, followed by selection of sliding semi-landmarks. 
#' Use function \code{\link{define.sliders.2d}} to select sliding semi-landmarks.
#' 
#' \subsection{Digitizing}{ 
#' Digitizing landmarks involves landmark selection using a mouse in the plot window, 
#' using the LEFT mouse button (or regular button for Mac users):
#' \enumerate{
#'  \item Digitize the scale bar by selecting the two end points (single click for start and end),
#'  \item Digitize each landmark with single click and the landmark is shown in red,
#'  When selection of n landmarks is completed, an ".nts" file is created in working directory using the specimen name, adding ".nts" as a suffix.
#' }
#' If verbose = TRUE, digitizing is interactive between landmark selection using a mouse and the R console. 
#' Once a landmark is selected, the user is asked if the system should keep or discard the 
#' selection (y/n). If "y", the user is asked to continue to select the next landmark. If "n", the user is asked to select it again. 
#' This can be repeated until the user is comfortable with the landmark chosen. verbose = FALSE is silent, and digitizing of each 
#' landmark is continuous and uninterupted until all landmarks are
#' chosen. Landmark coordinates are returned scaled.}
#' 
#' @param file Name of jpeg file to be digitized. File names can be 
#' written in manually, including paths, or obtained using directory/file manipulation functions 
#' e.g., \code{\link{list.files}}
#' @seealso \code{\link{list.files}}
#' @param nlandmarks Number of landmarks to be digitized.
#' @param scale Length of scale placed in image.
#' @param verbose logical. User decides whether to digitize in verbose or silent format (see details), default is verbose
#' @return Function returns an n-x-2 .nts file containing the 2d coordinates of digitized landmarks
#' @keywords digitizing
#' @export
#' @author Erik Otarola-Castillo and Emma Sherratt
digitize2d<-function(file, nlandmarks,scale,verbose=TRUE){
  spec.name<-unlist(strsplit(basename(file), "\\."))[1]
  specimen<-readJPEG(file, native = T)  
  plot(seq(0,dim(specimen)[2],length.out=10),seq(0,dim(specimen)[1],length.out=10), type='n',xlab="x",ylab="y",asp=1,tck=0,xaxt="n",yaxt="n")
  rasterImage(specimen, 1, 1,dim(specimen)[2],  dim(specimen)[1])
  cat("Set scale =",scale,"\n")
  scalebar<-picscale(scale)
  selected <- matrix(NA, nrow = nlandmarks, ncol = 2)
  fix <- NULL
  if(verbose==TRUE){
    for (i in 1:nlandmarks) {
      cat("Select landmark ", i , "\n")
      keep <- ans <- NULL
      fix<- locator(n=1, type ="p",col="black",cex=4,pch=21,bg="red")
      cat(paste("Keep Landmark ", i, "(y/n)?"), "\n")
      ans <- readLines(n = 1)
      if (ans == "y") {
        selected[i,1] <- fix$x
        selected[i,2] <- fix$y
      }
      if (ans == "n") {
        cat(paste("Select Landmark ", i, " Again"), "\n")
      }
      while (ans == "n") {
        fix<- locator(n=1, type ="p",col="black",cex=4,pch=21,bg="red")
        cat(paste("Keep Landmark ", i, "(y/n)?"), "\n")
        ans <- readLines(n = 1)
        if (ans == "y") { 
          selected[i,1] <- fix$x
          selected[i,2] <- fix$y
        }
        if (ans == "n") {
          cat(paste("Select Landmark ", i, " Again"), "\n")
        }
      }
    }
  }
  if(verbose==FALSE){
    cat("Select landmarks 1:",nlandmarks,"\n",sep="")
    for (i in 1:nlandmarks) {
      fix<- locator(n=1, type ="p",col="black",cex=4,pch=21,bg="red")
      selected[i,1] <- fix$x
      selected[i,2] <- fix$y  
    }
  }
  output<-selected/scalebar 
  writeland.nts(output, spec.name, comment=NULL)  
  output<-matrix(output,dimnames=list(paste("LM",seq(1,nlandmarks)),c("x","y")),ncol=2)
  return(LANDMARKS=output)
}