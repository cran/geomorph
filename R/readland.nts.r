#' Read landmark data matrix from nts file
#'
#' Read single *.nts file containing landmark coordinates for one or more specimens
#'
#' Function reads a single *.nts file containing two- or three-dimensional landmark coordinates. 
#' 
#' NTS files are text files in one of the standard formats for geometric morphometrics (see Rohlf 2012).
#' Multiple specimen format: 
#'   The parameter line contains 5 or 6 elements, and must begin with a "1" to designate a rectangular 
#'   matrix. The second and third values designate how many specimens (n) and how many total variables 
#'   (p x k) are in the data matrix. The fourth value is a "0" if the data matrix is complete and a "1" 
#'   if there are missing values. If missing values are present, the '1' is followed by the arbitrary 
#'   numeric code used to represent missing values (e.g., -999). These values will be replaced with "NA" 
#'   in the output array. Subsequent analyses requires a full complement of data, see \code{\link{estimate.missing}}. 
#'   The final value of the parameter line denotes the dimensionality of the landmarks
#'   (2,3) and begins with "DIM=". If specimen and variable labels are included, these are designated placing 
#'   an "L" immediately following the specimen or variable values in the parameter file. The labels then 
#'   precede the data matrix.
#'   
#'   Missing data may also be represented by designating them using 'NA'. In
#'   this case, the standard NTSYS header is used with no numeric designation for missing data (i.e. the fourth value is '0').
#'   The positions of missing landmarks may then be estimated using estimate.missing.

#'
#' Special NTS files: *.dta files in the written by IDAV Landmark Editor, 
#' and *.nts files written by Stratovan Checkpoint have incorrect 
#' header notation; every header is 1 n p-x-k 1 9999 Dim=3, rather than 1 n p-x-k 0 Dim=3, which denotes
#' that missing data is in the file even when it is not. Users must change manually the header (in a text editor) before using this function
#'
#' @param file the name of a *.nts file containing two- or three-dimensional landmark data to be read in
#' @keywords IO
#' @export
#' @author Dean Adams & Emma Sherratt
#' @seealso \code{\link{readmulti.nts}}
#' @return Function returns a 3D array (p x k x n), where p is the number of landmark points, k is the number 
#'   of landmark dimensions (2 or 3), and n is the number of specimens. The third dimension of this array 
#'   contains names for each specimen, which are obtained from the names in the *.nts file (if included). 
#' @references  Rohlf, F. J. 2012 NTSYSpc: Numerical taxonomy and multivariate analysis system. Version 
#'   2.2. Exeter Software, New York.


readland.nts <- function(file){    	
  ntsfile <- scan(file=file, what="char", quote="", sep="\n", strip.white=TRUE, comment.char="\"", quiet=TRUE)
  comment <- grep("\'", ntsfile)
  if (length(comment) != 0){
    ntsfile <- scan(file=file, what="char", quote="", sep="\n", strip.white=TRUE, comment.char="\'", quiet=TRUE)
  }
  
  header <- unlist(strsplit(ntsfile[1], "\\s+"))
  if(header[1]!=1){
    stop("NTS file not a rectangular matrix. First value in parameter line must be '1'.") }
  header <- casefold(header, upper=TRUE)
  
  dimval <- unlist(grep("DIM=", header))
  if(length(dimval)==0){
    stop("Header does not contain 'DIM=' designator.") }  
  
  labval <- unlist(grep("L", header))
  r.lab <- ifelse(is.element("2", labval)==TRUE, T, F)
  c.lab <- ifelse(is.element("3", labval)==TRUE, T, F)
  
  header <- sub("L", "", header)
  header <- as.numeric(sub("DIM=","", header))

  missdata <- ifelse(header[4]!=0, T, F)
  if(missdata==TRUE) {
    missval <- ifelse(dimval==6, header[5], header[6]) 
  }
  
  if(header[3] == header[dimval]){
    n <- 1; k <- header[dimval]; p <- header[2]
    } else {
    n <- header[2]; k <- header[dimval]; p <- header[3]/k
    }
  
  tmp <- unlist(strsplit(ntsfile[-1],"\\s+"))
  if(r.lab) {
    speclab <- tmp[1:n]
    tmp <- tmp[-(1:n)]
  } else speclab <- NULL
  
  if(c.lab) tmp <- tmp[-(1:(p*k))]
  
  if(missdata==TRUE) {tmp[grep(missval, as.integer(tmp))] <- NA}
  options(warn=-1)
  landdata <- matrix(as.numeric(tmp), ncol=k, byrow=TRUE)
  if(sum(which(is.na(landdata)==TRUE))>0){cat("NOTE.  Missing data identified.")}
  
  if(nrow(landdata) != p) {
    coords <- aperm(array(t(landdata), c(k,p,n)), c(2,1,3))
  } else {
    coords <- array(landdata, c(p, k, n))
  }
  
  if(length(speclab)==1) {
    dimnames(coords)[[3]] <- list(speclab) 
  } else { dimnames(coords)[[3]] <- speclab }
  return(coords)
}
