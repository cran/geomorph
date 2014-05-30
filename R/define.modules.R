#' Define modules (landmark partitions)
#' 
#' An interactive function to define which landmarks should be assigned to each module (landmark partition).
#' 
#' Function takes a matrix of two-dimensional digitized landmark coordinates and allows the user to assign 
#' landmarks to each module. The output is a list of which 
#' landmarks belong in which partition, to be used by \code{\link{compare.modular.partitions}}. 
#' The number of modules is chosen by the user (up to five). 
#'  
#'  \subsection{Selection}{ 
#' Choosing which landmarks will be included in each module involves landmark selection using a mouse in 
#' the plot window. The user is prompted to select each landmarks for module 1: using the LEFT mouse button 
#' (or regular button for Mac users), click on the hollow circle to choose the landmark. Selected landmarks 
#' will be filled in. When all landmarks for module 1 are chosen, press 'esc', and then start selecting
#' landmarks for module 2. Repeat until all modules are defined.
#' }
#' Note: Function currently only implemented for 2D landmark data.
#' 
#' @param spec Name of specimen, as an object matrix containing 2D landmark coordinates
#' @param nmodules Number of modules to be defined
#' @return Function returns a vector of which landmarks belong in which module (e.g. A,A,A,B,B,B,C,C,C) to be used
#' with \code{\link{compare.modular.partitions}} option 'landgroups'.
#' @export
#' @keywords utilities
#' @seealso  \code{\link{compare.modular.partitions}}
#' @author Emma Sherratt
#' 
define.modules <- function(spec, nmodules){
  spec.name <- deparse(substitute(spec))
  checkmat <- is.matrix(spec)
  if (checkmat == FALSE) {
    stop("Input must be a p-x-k matrix of landmark coordinates") }
  checkdim <- dim(spec)[2]
  if (checkdim == 3) {
    stop("Input must be a p-x-k matrix of two-dimensional landmark coordinates")
  }
  plot(spec[, 1], spec[, 2], cex = 1, pch = 21, bg = "white", 
       xlim = range(spec[, 1]), ylim = range(spec[, 2]), asp = 1)
  text(spec[, 1], spec[, 2], label = paste(1:dim(spec)[1]), 
       adj = 0.5, pos = 4)
  selected <- matrix(NA, nrow=nrow(spec), ncol=1)
  col <- c("red", "blue", "green", "purple", "brown")
  epit <- c("A","B","C","D", "E")
  for (i in 1:nmodules){
    cat("Select landmarks in module ",i,"\n",sep="")
    cat("Press esc when finished ","\n",sep="")
    select <- identifyPch(spec, col=col[i])
    selected[select] <- epit[i]
  }
  return(as.vector(selected))
}