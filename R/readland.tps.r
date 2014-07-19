#' Read landmark data from tps file
#'
#' Read *.tps file to obtain landmark coordinates
#'
#' This function reads a *.tps file containing two- or three-dimensional landmark coordinates. 
#'   Tps files are text files in one of the standard formats for geometric morphometrics (see Rohlf 2010). 
#'   Two-dimensional landmarks coordinates are designated by the identifier "LM=", while three-dimensional 
#'   data are designated by "LM3=". Landmark coordinates are multiplied by their scale factor if this is 
#'   provided for all specimens. If one or more specimens are missing the scale factor, landmarks are treated 
#'   in their original units.  
#'   
#'   Missing data may be present in the file. In this case, they must be designated by 'NA'. The 
#'   positions of missing landmarks may then be estimated using estimate.missing.
#' 
#' The user may specify whether specimen names are to be extracted from the 'ID=' field or 'IMAGE=' field 
#' and included in the resulting 3D array. 
#' e.g., for 'ID=' use (file, specID = "ID") and for 'IMAGE=' use (file, specID = "imageID").
#' The default is specID="None".
#' 
#' NOTE: At present, all other information that can be contained in tps files (curves, comments, variables, radii, etc.)
#'   is ignored. 
#'
#' @param file A *.tps file containing two- or three-dimensional landmark data
#' @param specID a character specifying whether to extract the specimen ID names from the ID or IMAGE lines (default is "None").
#' @param warnmsg A logical value stating whether warnings should be printed
#' @export
#' @keywords IO
#' @author Dean Adams & Emma Sherratt
#' @return Function returns a (p x k x n) array, where p is the number of landmark points, k is the number 
#'   of landmark dimensions (2 or 3), and n is the number of specimens. The third dimension of this array 
#'   contains names for each specimen, which are obtained from the image names in the *.tps file. 
#' @references  Rohlf, F. J. 2010. tpsRelw: Relative warps analysis. Version 1.49. Department of Ecology 
#'   and Evolution, State University of New York at Stony Brook, Stony Brook, NY.

readland.tps <- function (file, specID = c("None", "ID", "imageID"), warnmsg = T) 
{
  specID <- match.arg(specID)
  tpsfile <- scan(file = file, what = "char", sep = "\n", quiet = TRUE)
  lmdata <- grep("LM=", tpsfile)
  if (length(lmdata) == 0) {
    lmdata <- grep("LM3=", tpsfile)
    nland <- as.numeric(sub("LM3=", "", tpsfile[lmdata]))
    k <- 3
  }
  else {
    nland <- as.numeric(sub("LM=", "", tpsfile[lmdata]))
    k <- 2
  }
  n <- nspecs <- length(lmdata)
  if (max(nland) - min(nland) != 0) {
    stop("Number of landmarks not the same for all specimens.")
  }
  p <- nland[1]
  imscale <- as.numeric(sub("SCALE=", "", tpsfile[grep("SCALE", 
                                                       tpsfile)]))
  if (is.null(imscale)) {
    imscale = array(1, nspecs)
  }
  if (warnmsg == T) {
    if (length(imscale) != nspecs) {
      print("Not all specimens have scale. Assuming landmarks have been previously scaled.")
    }
  }
  if (length(imscale) != nspecs) {
    imscale = array(1, nspecs)
  }
  tmp <- tpsfile[-(grep("=", tpsfile))]
  options(warn = -1)
  tmp <- matrix(as.numeric(unlist(strsplit(tmp, split = " +")), 
                           ncol = k, byrow = T))
  if (warnmsg == T) {
    if (sum(which(is.na(tmp) == TRUE)) > 0) {
      print("NOTE.  Missing data identified.")
    }
  }
  coords <- aperm(array(t(tmp), c(k, p, n)), c(2, 1, 3))
  imscale <- aperm(array(rep(imscale, p * k), c(n, k, p)), 
                   c(3, 2, 1))
  coords <- coords * imscale
  
  if (specID == "None") {
      if (warnmsg == T) {print("No Specimen names extracted")
    }
  }
  if (specID == "imageID") {
    imageID <- (sub("IMAGE=", "", tpsfile[grep("IMAGE", tpsfile)]))
    if (length(imageID) != 0) {
      imageID <- sub(".jpg", "", imageID)
      imageID <- sub(".tif", "", imageID)
      imageID <- sub(".bmp", "", imageID)
      imageID <- sub(".tiff", "", imageID)
      imageID <- sub(".jpeg", "", imageID)
      imageID <- sub(".jpe", "", imageID)
      dimnames(coords)[[3]] <- as.list(imageID)
      if (warnmsg == T) {
        print("Specimen names extracted from line IMAGE=")
      }
    }
    if (length(imageID) == 0) {
      if (warnmsg == T) {
        print("No name given under IMAGE=. Specimen names not extracted")
      }
    } 
  }
  
  if (specID == "ID") {
    ID <- sub("ID=", "", tpsfile[grep("ID", tpsfile)])
    if (length(ID) == 0) {
      if(warnmsg ==T){
        print("No name given under ID=. Specimen names not extracted")
        }
      }
    if (length(ID) != 0) {
      dimnames(coords)[[3]] <- as.list(ID)
      if (warnmsg == T) {
        print("Specimen names extracted from line ID=")
      }
    }
  }
  return(coords = coords)
}