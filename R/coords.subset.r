#' Subset landmark coordinates via a factor
#'
#' Subset (split) landmark coordinates via a grouping factor
#'
#' This function splits a set of landmark coordinates into subsets, as described by a factor.  The 
#' results is a list of separate sets of landmarks.
#'
#' @param A A 3D array (p x k x n) containing landmark coordinates for a set of specimens
#' @param group A grouping factor of length n, for splitting the array into sub-arrays
#' @keywords utilities
#' @export
#' @author Michael Collyer
#' @examples
#' data(pupfish) 
#' group <- factor(paste(pupfish$Pop, pupfish$Sex))
#' levels(group)
#' new.coords <- coords.subset(A = pupfish$coords, group = group)
#' names(new.coords) # see the list levels
#' # access any element by:
#  # new.coords$`Marsh F` # Can be used in any analysis, just like pupfish$coords
#  # NOTE: `` surrounds level name because it has a space in it
#'
coords.subset <- function(A, group){
  dims <- dim(A)
  if(length(dims) != 3) stop("coordinates must be in the form of a 3D array - consider using arrayspecs")
  p <- dims[1]; k <- dims[2]; n <- dims[3]
  group <- as.factor(group)
  if(length(group) != n) stop("number of specimens do not match between coords and grouping factor")
  Y <- as.data.frame(two.d.array(A))
  X <- split(Y, group)
  redo <- function(x) arrayspecs(x, p, k)
  out <- lapply(X, redo)
  out
}	