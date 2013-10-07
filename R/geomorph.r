#' @name geomorph-package
#' @docType package
#' @aliases geomorph
#' @title Geometric morphometric analyses for 2D/3D data
#' @author Dean C. Adams, Erik Otarola-Castillo, & Emma Sherratt
#'
#' Functions in this package allow one to read, manipulate, and digitize landmark data; generate shape
#'  variables via Procrustes analysis for points, curves and surface data, perform statistical analyses
#'  of shape variation and covariation, and provide graphical depictions of shapes and patterns of
#'  shape variation.
#' 
#' @import rgl 
#' @importFrom geiger sim.char
#' @importFrom calibrate textxy
#' @importFrom jpeg readJPEG
#' @importFrom vegan adonis
#' @importFrom ape is.binary.tree ace compute.brlen vcv.phylo summary.phylo Ntip

NULL

#' Landmark data from Plethodon salamander heads
#'
#' @name plethodon
#' @docType data
#' @author Dean Adams
#' @references Adams, D. C. 2004. Character displacement via aggressive interference in Appalachian salamanders. 
#' Ecology. 85:2664-2670.
#' @references Adams, D.C. 2010. Parallel evolution of character displacement driven by competitive selection 
#' in terrestrial salamanders. BMC Evolutionary Biology. 10(72)1-10.
#' @keywords data
NULL

#' Landmark data from dataset rat
#'
#' @name ratland
#' @docType data
#' @author Dean Adams
#' @references Bookstein, F. L. 1991. Morphometric tools for landmark data: Geometry and Biology. 
#'  Cambridge Univ. Press, New York.
#' @keywords data
NULL

#' Landmark data from hummingbird bills (includes sliding semilandmarks on curves)
#'
#' @name hummingbirds
#' @docType data
#' @author Chelsea Berns and Dean Adams
#' @references Berns, C.M., and Adams, D.C. 2010. Bill shape and sexual shape dimorphism between two species 
#' of temperate hummingbirds: Archilochus alexandri (black-chinned hummingbirds) and Archilochus colubris 
#' (ruby-throated hummingbirds). The Auk. 127:626-635.
#' @keywords data
NULL

#' Average head shape and phylogenetic relationships for several Plethodon salamander species
#'
#' @name plethspecies
#' @docType data
#' @author Dean Adams
#' @references Phylogeny pruned from: Wiens et al. (2006). Evol.
#' @references Data from: Adams and Rohlf (2000); Adams et al. (2007); Arif et al. (2007) Myers and Adams (2008)
#' @keywords data
NULL

#' Landmark data from scallop shells
#'
#' @name scallops
#' @docType data
#' @author Dean Adams and Erik Otarola-Castillo
#' @references Serb et al. (2011). "Morphological convergence of shell shape in distantly related
#' scallop species (Mollusca: Pectinidae)." Zoological Journal of the Linnean Society 163: 571-584.
#' @keywords data
NULL

#' 3D scan of a scallop shell from a .ply file in mesh3d format
#'
#' @name scallopPLY
#' @docType data
#' @author Emma Sherratt
#' @references Serb et al. (2011). "Morphological convergence of shell shape in distantly related
#' scallop species (Mollusca: Pectinidae)." Zoological Journal of the Linnean Society 163: 571-584.
#' @keywords data
NULL

#' Simulated motion paths
#'
#' @name motionpaths
#' @docType data
#' @author Dean Adams
#' @references Adams, D. C., and M. L. Collyer. 2009. A general framework for the analysis of phenotypic 
#'   trajectories in evolutionary studies. Evolution 63:1143-1154.
#' @keywords data
NULL

#' landmarks on mosquito wings
#'
#' @name mosquito
#' @docType data
#' @author Dean Adams
#' @keywords data
NULL

#' Convert landmark data matrix into array (p x k x n)
#'
#' Convert a matrix of landmark coordinates into a 3-dimensional array 
#'
#' This function converts a matrix of landmark coordinates into a (p x k x n) 
#'  array, which is the required input format for many functions in geomorph. 
#'   Use \code{byLand}=TRUE if the input matrix is arranged such that the coordinates
#'   of each landmark are found on a separate row. Use \code{byLand}=FALSE if the 
#'   input matrix is arranged such that each row contains all landmark 
#'   coordinates for a single specimen.
#'
#' @param A A matrix containing landmark coordinates for a set of specimens
#' @param p Number of landmarks
#' @param k Number of dimensions (2 or 3)
#' @param byLand A logical value indicating whether rows of the input matrix contain
#' coordinates for each landmark separately, or all coordinates for each specimen
#' @export
#' @keywords arrayspecs
#' @author Dean Adams
#' @return Function returns a (p x k x n) array, where p is the number of landmark points, k is the number 
#'   of landmark dimensions (2 or 3), and n is the number of specimens. The third dimension of this array 
#'   contains names for each specimen if specified in the original input matrix.
#'  @examples 
#' x<-matrix(rnorm(18),nrow=3)  # Random triangles (all coordinates on same row for each triangle)
#' arrayspecs(x,3,2,byLand=FALSE) 
#'  
#' x2<-matrix(rnorm(18),ncol=2) # Random triangles (each landmark on its own row)
#' arrayspecs(x2,3,2,byLand=TRUE)
arrayspecs<-function(A,p,k,byLand=TRUE){	
  names<-rownames(A)
  n<-length(unlist(A))/(p*k)
  if(byLand==TRUE){
    n.m<-NULL 
    for(i in 1:n){
      temp<-as.matrix(A[((1+(i-1)*p):(i*p)),1:k])
      n.m<-cbind(n.m,temp)}
    specimens<-array(n.m,dim=c(p,k,n))
  }
  if(byLand==FALSE){
    n.m<-NULL 
    for(i in 1:n){
      temp<- matrix(A[i,],ncol=k,byrow=TRUE)
      n.m<-cbind(n.m,temp)}
    specimens<-array(n.m,dim=c(p,k,n))
  }
  if(length(names)==n){
    dimnames(specimens)[[3]]<-names}
  return(specimens)
}


#' Convert (p x k x n) data array into 2D data matrix
#'
#' Convert a three-dimensional array of landmark coordinates into a 2-dimensional matrix 
#'
#' This function converts a (p x k x n) array of landmark coordinates into a 2-dimensional 
#'  matrix. The latter format of the shape data is useful for performing subsequent statistical 
#'  analyses in R (e.g., PCA, MANOVA, PLS, etc.). Row labels are preserved if included in 
#'  the original array. 
#'
#' @param A An array (p x k x n) containing landmark coordinates for a set of specimens
#' @keywords two.d.array
#' @export
#' @author Dean Adams
#' @return Function returns a two-dimensional matrix of dimension (n x [p x k]), where rows 
#'   represent specimens and columns represent variables. 
#' @examples
#' data(plethodon) 
#' plethodon$land    #original data in the form of 3D array
#' 
#' two.d.array(plethodon$land)   # Convert to a 2D data matrix
two.d.array<-function(A){	
  newdata<-NULL
  for(i in 1:(dim(A)[3])){
    temp<-array(t(A[,,i]))
    newdata<-rbind(newdata,temp) }
  rownames(newdata)<-dimnames(A)[[3]]  
  return(newdata)
}

#' Read mesh data (vertices and faces) from ply files
#'
#' A function to read ply files to obtain landmark coordinates
#'
#' Function reads three-dimensional surface data in the form of a single ply file
#'  (Polygon File Format; ASCII format only, from 3D scanners such as NextEngine and David scanners). 
#'  Vertices of the surface may then be used to digitize three-dimensional points, 
#'  and semilandmarks on curves and surfaces. The function opens the ply file and plots the mesh,
#'  with faces rendered if file contains faces, and colored if the file contains vertex color.
#'
#' @param file An ASCII ply file
#' @param ShowSpecimen logical A logical value indicating whether or not the ply file should be displayed
#' @export
#' @keywords read.ply
#' @author Dean Adams & Emma Sherratt
#' @return Function returns the following components:
#'   \item{mesh3d}{list of class mesh3d- see rgl for details}
#' @examples 
#' # If the file has no mesh color, or color is undesirable, user can assign this as follows:
#' # Using the example scallop PLY
#' data(scallopPLY) 
#' myply <- scallopPLY$ply
#' myply$material <- "gray" # using color word
#' myply$material <- "#FCE6C9" # using RGB code
read.ply <- function (file, ShowSpecimen = TRUE) 
{
  plyfile <- scan(file = file, what = "char", sep = "\n", strip.white = TRUE, 
                  quiet = TRUE)
  is.ply <- grep("ply", plyfile)
  if ((length(is.ply) ==0) ) stop ("File is not a PLY file")
  format <- unlist(strsplit(grep(c("format "), plyfile, value = TRUE), " "))
  if (format[2] != "ascii")
    stop("PLY file is not ASCII format: ","format = ", format[2:length(format)])
  poly <- NULL
  material <- NULL
  xline <- unlist(strsplit(grep(c("vertex "), plyfile, value = TRUE), " "))
  npoints <- as.numeric(xline[grep(c("vertex"), xline) + 1])
  yline <- unlist(strsplit(grep(c("element face"), plyfile, value = TRUE), " "))
  npoly <- as.numeric(yline[grep(c("face"), yline) + 1])
  headerend <- grep(c("end_header"), plyfile)
  ncolpts <- (length(grep(c("property"), plyfile)) - 1)
  cols <- grep(c("property"), plyfile, value = TRUE)
  x <- grep(c(" x"), cols)
  y <- grep(c(" y"), cols)
  z <- grep(c(" z"), cols)
  points <- as.matrix(as.numeric(unlist(strsplit(plyfile[(headerend + 
                                                            1):(headerend + npoints)], " "))))
  dim(points) <- c(ncolpts, npoints)
  xpts <- points[x ,]
  ypts <- points[y ,]
  zpts <- points[z ,]
  vertices <- rbind(xpts, ypts, zpts, 1)
  if (yline[3] == 0)
    print("Object has zero faces")
  if (yline[3] != 0) {
    poly <- as.matrix(as.numeric(unlist(strsplit(plyfile[(headerend + 
                                                            npoints + 1):(headerend + npoints + npoly)], " "))))
    dim(poly) <- c((poly[1] + 1), npoly)
    poly <- poly[-1, ]
    poly = poly +1 }
  colinfo <- grep("property uchar red", plyfile)
  if (length(colinfo) !=0) {
    color <- rgb(points[4,], points[5,], points[6,], maxColorValue = 255)
    material$color <- matrix(color[poly], dim(poly))}   
  mesh <- list(vb = vertices, it = poly, primitivetype = "triangle", 
               material = material)
  class(mesh) <- c("mesh3d", "shape3d")
  if(ShowSpecimen==TRUE){ 
    clear3d()
    if (length(poly) == 0) { dot3d(mesh) }
    if (length(material) != 0){ shade3d(mesh) }
    shade3d(mesh, color="gray") }
  return(mesh)
}

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
#' The user may specify whether specimen names are to be included in the resulting 3D array. If so, names may
#' be read from either the 'ID' field or the 'IMAGE' field of the TPS file.
#'
#' NOTE: At present, all other information that can be contained in tps files (curves, comments, variables, radii, etc.)
#'   is ignored. 
#'
#' @param file A *.tps file containing two- or three-dimensional landmark data
#' @param specID a character specifying whether to extract the specimen ID names from the ID or IMAGE lines (default is "None").
#' @export
#' @keywords readland.tps
#' @author Dean Adams & Emma Sherratt
#' @return Function returns a (p x k x n) array, where p is the number of landmark points, k is the number 
#'   of landmark dimensions (2 or 3), and n is the number of specimens. The third dimension of this array 
#'   contains names for each specimen, which are obtained from the image names in the *.tps file. 
#' @references  Rohlf, F. J. 2010. tpsRelw: Relative warps analysis. Version 1.49. Department of Ecology 
#'   and Evolution, State University of New York at Stony Brook, Stony Brook, NY.
readland.tps<-function (file, specID=c("None","ID","imageID")){
  specID <- match.arg(specID)
  tpsfile <- scan(file = file, what = "char", sep = "\n", quiet = TRUE)
  lmdata <- grep("LM=", tpsfile)
    if (length(lmdata)==0){ 
      lmdata <- grep("LM3=", tpsfile) 
      nland <- as.numeric(sub("LM3=", "", tpsfile[lmdata]))
      k <- 3
      } else {nland <- as.numeric(sub("LM=", "", tpsfile[lmdata])) 
      k <- 2}
  n <- nspecs <- length(lmdata)
  if (max(nland) - min(nland) != 0) {stop("Number of landmarks not the same for all specimens.") }
  p <- nland[1]
  imscale <- as.numeric(sub("SCALE=", "", tpsfile[grep("SCALE", tpsfile)]))
    if (is.null(imscale)) {imscale = array(1, nspecs)}
    if (length(imscale) != nspecs) {print("Not all specimens have scale. Using scale = 1.0")}
    if (length(imscale) != nspecs) {imscale = array(1, nspecs)}
  landdata <- matrix(NA, nrow = (p * n), ncol = k)
  for (i in 1:n) {
    for (j in 1:p) {
      tmp <- gsub("\\t", " ", tpsfile[lmdata[i] + j])
      tmp <- as.numeric(unlist(strsplit(tmp, split = " +"))) *imscale[i]
      landdata[((i - 1) * p) + j, ] <- tmp; rownames(landdata) <- NULL
    }
  }
  coords <- arrayspecs(landdata, p, k)
  if(specID=="None"){   print("No Specimen names extracted") }
  if(specID=="imageID"){
    imageID <- (sub("IMAGE=", "", tpsfile[grep("IMAGE", tpsfile)]))
    if (length(imageID) ==0) {print("No name given under IMAGE=. Specimen names not extracted") }
    else if (length(imageID) !=0){
      imageID <- sub(".jpg", "", imageID); imageID <- sub(".tif", "", imageID)
      imageID <- sub(".bmp", "", imageID); imageID <- sub(".tiff", "", imageID)
      imageID <- sub(".jpeg", "", imageID); imageID <- sub(".jpe", "", imageID)
      dimnames(coords)[[3]] <- as.list(imageID) 
      print("Specimen names extracted from line IMAGE=") }
  }
  if(specID=="ID"){
    ID <-sub("ID=", "", tpsfile[grep("ID", tpsfile)])  
    if (length(ID) ==0) {print("No name given under ID=. Specimen names not extracted") }
    else if (length(ID) !=0){
      dimnames(coords)[[3]] <- as.list(ID)  
      print("Specimen names extracted from line ID=") }
  } 
  return(coords = coords)
}

#' Write landmark data to tps file
#'
#' Write *.tps file from obtain landmark coordinates in a 3-dimensional array 
#'
#' This function writes a *.tps file from a 3-dimensional array (p x k x n) 
#'  of landmark coordinates. 
#'
#' @param A An array (p x k x n) containing landmark coordinates for a set of specimens
#' @param file Name of the *.tps file to be created
#' @export
#' @keywords writeland.tps
#' @author Dean Adams
writeland.tps<-function(A, file){
  n<-dim(A)[3]
  k<-dim(A)[2]
  p<-dim(A)[1]
  lmline<-ifelse(k==2,paste("LM=",p,sep=""), paste("LM3=",p,sep=""))  
  file.create(file, overwrite = TRUE)
  for(i in 1:n){
    write(lmline,file,append = TRUE)
    write.table(A[,,i],file,col.names = FALSE, row.names = FALSE,append=TRUE)
    if(is.null(dimnames(A)[[3]])==FALSE){
      idline<-paste("ID=",dimnames(A)[[3]][i],sep="")
      write(idline,file,append = TRUE)  
    }
    write("",file,append = TRUE)
  }
}

#' Read landmark data matrix from nts file
#'
#' Read *.nts file to obtain landmark coordinates for a set of specimens
#'
#' Function reads a *.nts file containing a matrix of two- or three-dimensional landmark coordinates. 
#'   NTS files are text files in one of the standard formats for geometric morphometrics (see Rohlf 2012). 
#'   The parameter line contains 5 or 6 elements, and must begin with a "1" to designate a rectangular 
#'   matrix. The second and third values designate how many specimens (n) and how many total variables 
#'   (p x k) are in the data matrix. The fourth value is a "0" if the data matrix is complete and a "1" 
#'   if there are missing values. If missing values are present, the '1' is followed by the arbitrary 
#'   numeric code used to represent missing values (e.g., -999). These values will be replaced with "-999" 
#'   in the output array. The final value of the parameter line denotes the dimensionality of the landmarks
#'   (2,3) and begins with "DIM=". If specimen and variable labels are included, these are designated placing 
#'   an "L" immediately following the specimen or variable values in the parameter file. The labels then 
#'   precede the data matrix.
#'
#' Function is for *.nts file containing landmark coordinates for multiple specimens. Function can
#' read *.dta files in the nts format written by Landmark Editor \url{http://graphics.idav.ucdavis.edu/research/projects/EvoMorph}.
#'
#' @param file A *.nts file containing two- or three-dimensional landmark data
#' @keywords readland.nts
#' @export
#' @author Dean Adams & Emma Sherratt
#' @return Function returns a (p x k x n) array, where p is the number of landmark points, k is the number 
#'   of landmark dimensions (2 or 3), and n is the number of specimens. The third dimension of this array 
#'   contains names for each specimen, which are obtained from the image names in the *.nts file. 
#' @references  Rohlf, F. J. 2012 NTSYSpc: Numerical taxonomy and multivariate analysis system. Version 
#'   2.2. Exeter Software, New York.
readland.nts<-function(file){  		
  ntsfile<-scan(file=file,what="char",quote="",sep="\n",strip.white=TRUE,comment.char="\"",quiet=TRUE)
  comment <- grep("\'", ntsfile)
  if (length(comment) != 0){
    ntsfile<-scan(file=file,what="char",quote="",sep="\n",strip.white=TRUE,comment.char="\'",quiet=TRUE)
  }
  header<-unlist(strsplit(ntsfile[1]," "))
  if(header[1]!=1){
    stop("NTS file not a rectangular matrix. First value in parameter line must be '1'.") }
  header<-casefold(header,upper=TRUE)
  dimval<-unlist(grep("DIM=",header))
  if(length(dimval)==0){
    stop("Header does not contain 'DIM=' designator.") }  
  labval<-unlist(grep("L",header))
  r.lab<-ifelse(is.element("2",labval)==TRUE,T,F)
  c.lab<-ifelse(is.element("3",labval)==TRUE,T,F)
  header<-sub("L","",header)
  header<-as.numeric(sub("DIM=","",header))
  missdata<-ifelse(header[4]!=0,T,F)
  if(missdata==TRUE){missval<-ifelse(dimval==6,header[5],header[6]) } 
  n<-header[2];k<-header[dimval];p<-header[3]/k;   
  tmp<-gsub("\\t"," ",ntsfile[-1])  #replace tabs with spaces
  tmp<-unlist(strsplit(tmp,split=" +"))
  speclab<-NULL; 
  if(r.lab==TRUE){
    speclab<-tmp[1:n]
    tmp<-tmp[-(1:length(speclab))]   }
  if(c.lab==TRUE){ tmp<-tmp[-(1:(p*k))] }
  if(missdata==TRUE){tmp<-sub(missval,"-999",tmp)}
  landdata<-matrix(as.numeric(tmp),ncol=k,byrow=TRUE)
  coords<-arrayspecs(landdata,p,k,byLand=TRUE)
  dimnames(coords)[[3]]<-speclab
  return(coords=coords)
}

#' Read landmark data from multiple nts files
#'
#' Read a list of names for several *.nts files to obtain landmark coordinates for a set of specimens
#'
#' This function reads a list containing the names of multiple *.nts files, where each contains the 
#'   landmark coordinates for a single specimen. For these files, the number of variables (columns) of 
#'   the data matrix will equal the number of dimensions of the landmark data (k=2 or 3). When the function 
#'   is called a dialog box is opened, from which the user may select multiple *.nts files. These are then read 
#'   and concatenated into a single matrix for all specimens. 
#'
#' @param filelist A list of names for the *.nts files to be read by the function. The names in the list
#'   require quotes (") and .nts/.NTS suffix.
#' @keywords readmulti.nts
#' @export
#' @author Dean Adams & Emma Sherratt
#' @return Function returns a (p x k x n) array, where p is the number of landmark points, k is 
#'   the number of landmark dimensions (2 or 3), and n is the number of specimens. The third dimension 
#'   of this array contains names for each specimen, which are obtained from the original file names. 
readmulti.nts<-function(filelist){   
  n<-length(filelist)
  names<-gsub(".nts","",filelist)
  names<-gsub(".NTS","",names)
  landdata<-nind<-NULL
  for (i in 1:n){
    ntsfile<-scan(file=filelist[i],what="char",quote="",sep="\n",strip.white=TRUE,comment.char="\"",quiet=TRUE)
    comment <- grep("\'", ntsfile)
    if (length(comment) != 0){
      ntsfile<-scan(file=file,what="char",quote="",sep="\n",strip.white=TRUE,comment.char="\'",quiet=TRUE)
    }
    header<-unlist(strsplit(ntsfile[1]," "))
    if(header[1]!=1){
      stop("NTS file not a rectangular matrix. First value in parameter line must be '1'.") }
    header<-casefold(header,upper=TRUE)
    dimval<-unlist(grep("DIM=",header))
    if(length(dimval)==0){
      stop("Header does not contain 'DIM=' designator.") }  
    labval<-unlist(grep("L",header))
    r.lab<-ifelse(is.element("2",labval)==TRUE,T,F)
    c.lab<-ifelse(is.element("3",labval)==TRUE,T,F)
    header<-sub("L","",header)
    header<-as.numeric(sub("DIM=","",header))
    missdata<-ifelse(header[4]!=0,T,F)
    if(missdata==TRUE){missval<-ifelse(dimval==6,header[5],header[6]) } 
    p<-header[2];k<-header[3]
    nind<-rbind(nind,p)
    if (min(nind)!=max(nind)) {
      stop("Number of landmarks not the same in all files.") } 
    tmp<-gsub("\\t"," ",ntsfile[-1])  #replace tabs with spaces
    tmp<-unlist(strsplit(tmp,split=" +"))
    rowlab<-NULL; 
    if(r.lab==TRUE){
      rowlab<-tmp[1:p]
      tmp<-tmp[-(1:length(rowlab))]   }
    if(c.lab==TRUE){ tmp<-tmp[-(1:k)] }
    if(missdata==TRUE){tmp<-sub(missval,"-999",tmp)}
    data<-matrix(as.numeric(tmp),ncol=k,byrow=TRUE)
    landdata<-rbind(landdata,data)
  }
  coords<-arrayspecs(landdata,p,k,byLand=TRUE)
  dimnames(coords)[[3]]<-names
  return(coords=coords)
}

#' Read 3D landmark data from Morphologika files
#'
#' Read Morphologika files to obtain 3D landmark coordinates and specimen information
#'
#' This function reads the commonly used Morphologika file format for 2D and 3D landmarks.
#'  If semilandmarks are present in the file, after import these should be defined 
#'  using \code{\link{define.sliders.2d}} or \code{\link{define.sliders.3d}}.
#'  Landmark coordinates and specimen information may then be used to conduct GPA using \code{\link{gpagen}}, 
#' 
#' If the Morphologika headers "[labels]" and "[labelvalues]" are present, then a data matrix containing all 
#'  individual specimen information is returned. If "[labels]" and "[labelvalues]" are not present in the file, 
#'  then only landmark coordinates are returned.
#' 
#' @param file A Morphologika text file. File name can be written in manually, including path, or obtained using 
#' directory/file manipulation functions e.g., \code{\link{list.files}}
#' @seealso \code{\link{list.files}} 
#' @param plot Logical should (unaligned) specimens be plotted? Defaults to FALSE.
#' @export
#' @keywords read.morphologika
#' @author Erik Otarola-Castillo
#' @return Function returns a list with the following components:
#'   \item{coords}{If Morphologika headers "[labels]" and "[labelvalues]" are not present in the file, function 
#'   only returns a (p x k x n) array, where p is the number of landmark points, k is the number 
#'   of landmark dimensions, and n is the number of specimens}
#'   \item{dataframe}{If Morphologika headers "[labels]" and "[labelvalues]" are present read.morphologika returns 
#'   the above p x k x n array and a dataframe containing specimen specific information stored in file}
read.morphologika<-function(file,plot=FALSE){
  mfile<-scan(file,what="character",sep = "\n",strip.white = TRUE,quiet=TRUE)
  tab<-length(grep("\t",mfile))
  com<-length(grep(",",mfile))
  sem<-length(grep(";",mfile))
  spac1 <- length(grep("  ", mfile))
  spac2 <- length(grep("   ", mfile))
  spac3 <- length(grep("    ", mfile))
  spac4 <- length(grep("     ", mfile))
  if(tab>0){mfile<-gsub("\t"," ",mfile)}
  if(com>0){mfile<-gsub(","," ",mfile)}
  if(sem){mfile<-gsub(";"," ",mfile)} 
  if (spac1 > 0) {mfile <- gsub("  ", " ", mfile)}
  if (spac2 > 0) {mfile <- gsub("   ", " ", mfile)}
  if (spac3 > 0) {mfile <- gsub("    ", " ", mfile)}
  if (spac4 > 0) {mfile <- gsub("     ", " ", mfile)}
  spac5 <- length(grep("  ", mfile))
  if (spac5 > 0) {mfile <- gsub("  ", " ", mfile)}
  tags.obs<-mfile[grep(']', mfile)]
  if(length(which(tags.obs=="[wireframe]"))>0){tags.obs<-tags.obs[-which(tags.obs=="[wireframe]")]}
  if(length(which(tags.obs=="[polygon]"))>0){tags.obs<-tags.obs[-which(tags.obs=="[polygon]")]}
  if(sum(c("[labels]" ,"[labelvalues]") %in%tags.obs)!=2){attributes<-FALSE} else {attributes<-TRUE}
  inds<-as.numeric(mfile[which(mfile=="[individuals]" | mfile=="[Individuals]") + 1])
  lms<-as.numeric(mfile[which(mfile=="[landmarks]" | mfile=="[Landmarks]") +1])
  dims<-as.numeric(mfile[which(mfile=="[dimensions]" | mfile=="[Dimensions]")+1])
  coords<-as.numeric(which(mfile=="[rawpoints]" | mfile=="[Rawpoints]") +1)
  names<-as.numeric(which(mfile=="[names]" | mfile=="[Names]") +1)
  strfun<-function(dat,a){unlist(strsplit(dat[a]," "))}
  if (attributes==TRUE){
    label<-as.numeric(which(mfile=="[labels]" | mfile=="[Labels]") + 1)
    labvals<-as.numeric(which(mfile=="[labelvalues]" | mfile=="[Labelvalues]") +1)
    labs<-unlist(strsplit(mfile[label]," "))
    dat<-mfile[labvals:(labvals + (inds-1))]  
    info<-as.data.frame(t(sapply(1:inds,strfun,dat=dat)))
    colnames(info)<-labs  
    incoords<-(coords + 1)
  }
  if(length(which(tags.obs=="[names]"))>0){
    coords2<-array(NA,c(lms,dims,inds),dimnames=list(NULL,NULL,mfile[names:length(mfile)][1:(grep(']',mfile[names:length(mfile)])[1]-1)]))
  } else {
    coords2<-array(NA,c(lms,dims,inds),dimnames=list(NULL,NULL,rep(paste("Specimen",1:inds))))
  }  
  if(suppressWarnings(sum(is.na(as.numeric(as.matrix(unlist(strsplit(mfile[coords]," ")))))))>0){
    incoords<-(coords + 1)
    coords3<-matrix(NA,nrow=inds,ncol=dims*lms)
    for(i in 1:inds){
      incoords<-c(incoords,incoords[i] + lms + 1)
      coords2[,,i]<-matrix(as.numeric(sapply(1:lms,strfun,dat=as.matrix(mfile[incoords[i]:((incoords[i]+lms-1))]))),ncol=dims,byrow=TRUE)    
      if(attributes==TRUE){coords3[i,]<-as.numeric(sapply(1:lms,strfun,dat=as.matrix(mfile[incoords[i]:((incoords[i]+lms-1))])))}       
    } 
  }  
  if(suppressWarnings(sum(is.na(as.numeric(as.matrix(unlist(strsplit(mfile[coords]," ")))))))==0){
    incoords<-(coords)
    coords3<-matrix(NA,nrow=inds,ncol=dims*lms)
    for(i in 1:inds){
      incoords<-c(incoords,incoords[i] + lms)
      coords2[,,i]<-matrix(as.numeric(sapply(1:lms,strfun,dat=as.matrix(mfile[incoords[i]:((incoords[i]+lms-1))]))),ncol=dims,byrow=TRUE)    
      if(attributes==TRUE){coords3[i,]<-as.numeric(sapply(1:lms,strfun,dat=as.matrix(mfile[incoords[i]:((incoords[i]+lms-1))])))}        
    } 
  }  
  if(plot==TRUE){
    coords23<-coords2
    coords23[which(coords23==9999.99)]<-NA
    plot3d(coords23[,1,1],coords23[,2,1],coords23[,3,1],type="n",aspect=FALSE)
    for(i in 1:dim(coords2)[3]){
      points3d(coords23[,1,i],coords23[,2,i],coords23[,3,i],size=5)
    }
  }
  if (attributes==TRUE){
    info.frame<-cbind(info,coords3)
    return(list(coords=coords2,dataframe=info.frame))
  } else {return(coords=coords2)}
}

#' Estimate locations of missing landmarks using the thin-plate spline
#'
#' A function for estimating the locations of missing landmarks 
#' 
#' The function estimates the locations of missing landmarks on a target specimen based on the locations of 
#'  the corresponding landmarks on a complete, reference specimen. Missing landmarks in a target specimen are 
#'  designated by '-999' in place of the x,y,z coordinates.  First, a complete reference specimen is 
#'  aligned to a target specimen (missing one or more landmarks), using the set of landmarks common to both specimens. 
#'  Next, the thin-plate spline interpolating function is used to estimate the locations of the missing landmarks in the target specimen. 
#' 
#' @param ref Landmark coordinates of a reference specimen
#' @param target Landmark coordinates of a target specimen (contains -999 for missing landmarks)
#' @author Dean Adams
#' @return Function returns a n x p matrix of coordinates for the target specimen that includes the original landmarks
#' plus the estimated coordinates for the missing landmarks. 
#' @export
#' @references  Bookstein, F. L., K. Schafer, H. Prossinger, H. Seidler, M. Fieder, G. Stringer, G. W. Weber, 
#' J.-L. Arsuaga, D. E. Slice, F. J. Rohlf, W. Recheis, A. J. Mariam, and L. F. Marcus. 1999. Comparing 
#' frontal cranial profiles in archaic and modern Homo by morphometric analysis. Anat. Rec. (New Anat.) 257:217-224.
#' @references Gunz, P., P. Mitteroecker, S. Neubauer, G. W. Weber, and F. L. Bookstein. 2009. Principles for 
#' the virtual reconstruction of hominin crania. J. Hum. Evol. 57:48-62.
#' @examples
#' data(plethodon) 
#' ref<-plethodon$land[,,1]
#' target<-plethodon$land[,,11]
#' target[2,]<-target[6,]<- -999    #create some missing landmarks
#' estimate.missing(ref,target)
estimate.missing<-function(ref,target){ 
  M1<-ref;M2<-target
  missing<-which(M2[,1]==-999)
  M2<-tps2d3d(M1[,],M1[-missing,],M2[-missing,])
  return(M2)
}

#' Rotate a subset of 2D landmarks to common articulation angle
#'
#' A function for rotating a subset of landmarks so that the articulation angle between subsets is constant
#' 
#' This function standardizes the angle between two subsets of landmarks for a set of specimens. The approach assumes a simple
#'  hinge-point articulation between the two subsets, and rotates all specimens such that the angle between landmark subsets 
#'  is equal across specimens (see Adams 1999).  As a default, the mean angle is used, though the user may specify an additional amount by which 
#'  this may be augmented. 
#' 
#' Presently, the function is only implemented for two-dimensional landmark data. 
#' @param A An array (p x k x n) containing landmark coordinates for a set of specimens
#' @param art.pt A number specifying which landmark is the articulation point between the two landmark subsets
#' @param angle.pts A vector containing numbers specifying which two points used to define the angle (one per subset)
#' @param rot.pts A vector containing numbers specifying which landmarks are in the subset to be rotated
#' @param angle An optional value specifying the additional amount by which the rotation should be augmented (in radians)
#' @param degrees A logical value specifying whether the additional rotation angle is expressed in degrees or radians (radians is default)
#' @author Dean Adams
#' @return Function returns a (p x k x n) array of landmark coordinates. 
#' @export
#' @references Adams, D. C. 1999. Methods for shape analysis of landmark data from articulated structures. Evolutionary Ecology Research. 1:959-970.
#' @examples
#' #Example using Plethodon
#' #Articulation point is landmark 1, rotate mandibular landmarks (2-5) relative to cranium
#'
#' data(plethspecies) 
#' fixed.angle(plethspecies$land,art.pt=1,angle.pts=c(5,6),rot.pts=c(2,3,4,5))
fixed.angle<-function(A,art.pt=NULL,angle.pts=NULL,rot.pts=NULL,angle=0,degrees = FALSE){
  if (length(dim(A))!=3){
    stop("Data matrix 1 not a 3D array (see 'arrayspecs').")  }
  if(length(grep("-999",A))!=0){
    stop("Data matrix 1 contains missing values. Estimate these first(see 'estimate.missing').")  }
  n<-dim(A)[3];   k<-dim(A)[2];  p<-dim(A)[1]
  if (k!=2){
    stop("Method presently implemented for two-dimensional data only.")}
  if (degrees == FALSE & angle>pi*2){
    stop("Angle was specified in degrees (use degrees = T).")}
  if (degrees == FALSE & angle< -pi*2){
    stop("Angle was specified in degrees (use degrees = T).")}
  angl<-array(NA,dim=n)
  for (i in 1:n){   
    A[,,i]<-t(t(A[,,i])-A[art.pt,,i])  
    angl[i]<-  acos((A[angle.pts[1],,i]/sqrt(sum(A[angle.pts[1],,i]^2)))%*%
                      (A[angle.pts[2],,i]/sqrt(sum(A[angle.pts[2],,i]^2))))
  }  
  if(degrees ==TRUE)  angle = angle*pi/180
  dev.angle<- (angl-mean(angl))+angle
  if(A[angle.pts[1],1,1]<0){dev.angle<- -1*dev.angle}
  for (i in 1:n){   
    r = matrix(c(cos(dev.angle[i]),-sin(dev.angle[i]),
                 sin(dev.angle[i]),cos(dev.angle[i])),2) 
    A[rot.pts,,i] = A[rot.pts,,i]%*%r 
  }  
  return(A)
}


#' Generalized Procrustes analysis of points, curves, and surfaces
#'
#' A general function to perform Procrustes analysis of two- or three-dimensional landmark data that 
#'  can include both fixed landmarks and sliding semilandmarks
#'
#' The function performs a Generalized Procrustes Analysis (GPA) on two-dimensional or three-dimensional
#'  landmark coordinates. The analysis can be performed on fixed landmark points, semilandmarks on 
#'  curves, semilandmarks on surfaces, or any combination. To include semilandmarks on curves, one 
#'  must specify a matrix defining which landmarks are to be treated as semilandmarks using the "curves=" 
#'  option. Likewise, to include semilandmarks 
#'  on surfaces, one must specify a vector listing which landmarks are to be treated as surface semilandmarks 
#'  using the "surfaces=" option. The "ProcD=TRUE" option will slide the semilandmarks along their tangent 
#'  directions using the Procrustes distance criterion, while "ProcD=FALSE" will slide the semilandmarks 
#'  based on minimizing bending energy. The aligned Procrustes residuals can be projected into tangent 
#'  space using the "Proj=TRUE" option. NOTE: Large datasets may exceed the memory limitations of R. 
#'
#'  Generalized Procrustes Analysis (GPA: Gower 1975, Rohlf and Slice 1990) is the primary means by which 
#'   shape variables are obtained from landmark data (for a general overview of geometric morphometrics see 
#'   Bookstein 1991, Rohlf and Marcus 1993, Adams et al. 2004, Zelditch et al. 2004, Mitteroecker and 
#'   Gunz 2009, Adams et al. 2012). GPA translates all specimens to the origin, scales them to unit-centroid 
#'   size, and optimally rotates them (using a least-squares criterion) until the coordinates of corresponding
#'   points align as closely as possible. The resulting aligned Procrustes coordinates represent the shape 
#'   of each specimen, and are found in a curved space related to Kendall's shape space (Kendall 1984). 
#'   Typically, these are projected into a linear tangent space yielding Kendall's tangent space coordinates 
#'   (Dryden and Mardia 1993, Rohlf 1999), which are used for subsequent multivariate analyses. Additionally, 
#'   any semilandmarks on curves and are slid along their tangent directions or tangent planes during the 
#'   superimposition (see Bookstein 1997; Gunz et al. 2005). Presently, two implementations are possible: 
#'   1) the locations of semilandmarks can be optimized by minimizing the bending energy between the 
#'   reference and target specimen (Bookstein 1997), or by minimizing the Procrustes distance between the two 
#'   (Rohlf 2010).
#'
#' @param A An array (p x k x n) containing landmark coordinates for a set of specimens
#' @param Proj A logical value indicating whether or not the aligned Procrustes residuals should be projected 
#'   into tangent space 
#' @param ProcD A logical value indicating whether or not Procrustes distance should be used as the criterion
#'   for optimizing the positions of semilandmarks
#' @param curves An optional matrix defining which landmarks should be treated as semilandmarks on boundary 
#'   curves, and which landmarks specify the tangent directions for their sliding
#' @param pointscale An optional value defining the size of the points for all specimens
#' @param surfaces An optional vector defining which landmarks should be treated as semilandmarks on surfaces
#' @param ShowPlot A logical value indicating whether or not a plot of Procrustes residuals should be displayed
#' @keywords gpagen
#' @export
#' @author Dean Adams
#' @return Function returns a list with the following components: 
#'   \item{coords}{A (p x k x n) array of aligned Procrustes coordinates, where p is the number of landmark 
#'     points, k is the number of landmark dimensions (2 or 3), and n is the number of specimens. The third 
#'     dimension of this array contains names for each specimen if specified in the original input array}
#'   \item{Csize}{A vector of centroid sizes for each specimen, containing the names for each specimen if 
#'     specified in the original input array}
#' @references  Adams, D. C., F. J. Rohlf, and D. E. Slice. 2004. Geometric morphometrics: ten years of 
#'    progress following the 'revolution'. It. J. Zool. 71:5-16.
#' @references Adams, D. C., F. J. Rohlf, and D. E. Slice. 2012. In Press. A field comes of age: Geometric 
#'   morphometrics in the 21st century. Hystrix.
#' @references Bookstein, F. L. 1991. Morphometric tools for landmark data: Geometry and Biology. 
#'  Cambridge Univ. Press, New York.
#' @references Bookstein, F. L. 1997. Landmark methods for forms without landmarks: morphometrics of 
#'   group differences in outline shape.  1:225-243.
#' @references Dryden, I. L., and K. V. Mardia. 1993. Multivariate shape analysis. Sankhya 55:460-480.
#' @references Gower, J. C. 1975. Generalized Procrustes analysis. Psychometrika 40:33-51.
#' @references Gunz, P., P. Mitteroecker, and F. L. Bookstein. 2005. semilandmarks in three dimensions. 
#'   Pp. 73-98 in D. E. Slice, ed. Modern morphometrics in physical anthropology. Klewer Academic/Plenum, New York.
#' @references Kendall, D. G. 1984. Shape-manifolds, Procrustean metrics and complex projective spaces. 
#'   Bulletin of the London Mathematical Society 16:81-121.
#' @references Mitteroecker, P., and P. Gunz. 2009. Advances in geometric morphometrics. Evol. Biol. 36:235-247.
#' @references Rohlf, F. J., and D. E. Slice. 1990. Extensions of the Procrustes method for the optimal 
#'   superimposition of landmarks. Syst. Zool. 39:40-59.
#' @references Rohlf, F. J., and L. F. Marcus. 1993. A revolution in morphometrics. Trends Ecol. Evol. 8:129-132.
#' @references Rohlf, F. J. 1999. Shape statistics: Procrustes superimpositions and tangent spaces. 
#'   Journal of Classification 16:197-223.
#' @references Rohlf, F. J. 2010. tpsRelw: Relative warps analysis. Version 1.49. Department of Ecology and 
#'   Evolution, State University of New York at Stony Brook, Stony Brook, NY.
#' @references Zelditch, M. L., D. L. Swiderski, H. D. Sheets, and W. L. Fink. 2004. Geometric morphometrics 
#'   for biologists: a primer. Elsevier/Academic Press, Amsterdam.
#' @examples
#' #Example 1: fixed points only
#' data(plethodon) 
#' gpagen(plethodon$land)
#' points(mshape(gpagen(plethodon$land)$coords),pch=22,col="red",bg="red",cex=1.2)    
#' 
#' #Example 2: points and semilandmarks on curves
#' data(hummingbirds)
#'
#' #Matrix defining which points are semilandmarks (middle column) 
#' #and in which directions they slide (columns 1 vs. 3)
#'
#' hummingbirds$curvepts    
#'
#' #Using Procrustes Distance for sliding
#' gpagen(hummingbirds$land,curves=hummingbirds$curvepts,)   
#' 
#' #Using bending energy for sliding
#' gpagen(hummingbirds$land,curves=hummingbirds$curvepts,ProcD=FALSE,)   
#'
#' #Example 3: points, curves and surfaces
#' data(scallops)
#' #Using Procrustes Distance for sliding
#' gpagen(A=scallops$coorddata, curves=scallops$curvslide, surfaces=scallops$surfslide) 
#' @useDynLib geomorph
gpagen<-function(A, Proj=TRUE,ProcD=TRUE,ShowPlot=TRUE,curves = NULL, surfaces = NULL,pointscale=1){
  if (length(dim(A))!=3){
    stop("Data matrix not a 3D array (see 'arrayspecs').")  }
  if(length(grep("-999",A))!=0){
    stop("Data matrix contains missing values. Estimate these first(see 'estimate.missing').")  }
  n<-dim(A)[3];   k<-dim(A)[2];  p<-dim(A)[1]
  nsurf<-0; ncurve<-0
  slided<-ifelse(ProcD==T,1,0)
  if(!is.null(curves)){ncurve<-nrow(curves) }
  if(!is.null(surfaces)){nsurf<-nrow(surfaces) }
  specs.size<-NULL
    for (i in 1:n)
      {specs.size[i]<-csize(A[,,i])[[1]]}
  if (is.null(curves) && is.null(surfaces)){	
    temp<-.C("DoGPA", as.integer(p),as.integer(k),as.integer(n),as.double(A),double(p*k*n),PACKAGE = "geomorph")[[5]]
    temp<-arrayspecs(matrix(temp,ncol=k,byrow=T),p,k,byLand=T)
  }
 else{ 
   temp<-.C("DoGPA1", as.integer(p),as.integer(k),as.integer(n),as.double(A),double(p*k*n),PACKAGE = "geomorph")[[5]]
   temp<-arrayspecs(matrix(temp,ncol=k,byrow=T),p,k,byLand=T)
   ref.gpa<-mshape(temp)	
   new.gpa<-.C("DoSlide", as.integer(slided), as.integer(p),as.integer(k),as.integer(n),
    as.double(temp),as.double(ref.gpa), double(p*k*n),as.integer(ncurve),as.integer(curves),as.integer(nsurf),
    as.integer(surfaces),PACKAGE = "geomorph" )[[7]]
    new.gpa<-arrayspecs(matrix(new.gpa,ncol=k,byrow=T),p,k,byLand=T)
    temp<-.C("DoGPA1", as.integer(p),as.integer(k),as.integer(n),as.double(new.gpa),double(p*k*n),PACKAGE = "geomorph")[[5]]
    temp<-arrayspecs(matrix(temp,ncol=k,byrow=T),p,k,byLand=T)
  }
  if(Proj==TRUE){temp<-orp(temp)  }
  dimnames(temp)[[3]]<-dimnames(A)[[3]]
  names(specs.size)<-dimnames(A)[[3]]
  ptsz<-pointscale
  if(ShowPlot==TRUE){ plotAllSpecimens(temp,pointscale=ptsz)}
  return(list(coords=temp,Csize=specs.size))
}

#' Quantify morphological integration between two modules
#'
#' Function quantifies the degree of morphological integration between two modules of Procrustes-aligned 
#'   coordinates
#'
#' The function quantifies the degree of morphological integration between two modules of shape data as 
#'   defined by landmark coordinates. It is assumed that the landmarks have previously been aligned using 
#'   Generalized Procrustes Analysis (GPA) [e.g., with \code{\link{gpagen}}]. The function may be used to assess
#'   the degree of morphological integration between two separate structures or between two modules defined within 
#'   the same landmark configuration. 
#'   
#'   Two analytical approaches are currently implemented to assess the degree of morphological 
#'   integration. If "method=PLS" (the default) the function estimates the degree of morphological 
#'   integration using two-block partial least squares, or PLS. When used with landmark data, this analysis 
#'   is referred to as singular warps analysis (Bookstein et al. 2003). When "method=PLS", the scores along 
#'   the X & Y PLS axes are also returned.  Alternatively, if "method=RV" the 
#'   function estimates the degree of morphological integration using the RV coefficient (Klingenberg 2009). 
#'   Significance testing for both approaches is found by permuting the objects in one data matrix relative 
#'   to those in the other. A histogram of coefficients obtained via resampling is presented, with the 
#'   observed value designated by an arrow in the plot.
#'   
#'    If evaluating an a priori hypothesis of modularity within a structure is of interest, one may use the 
#'    average RV coefficient as implemented in the function \code{\link{compare.modular.partitions}}. 
#'
#' @param A An array (p1 x k x n) containing landmark coordinates for the first module
#' @param A2 An array (p2 x k x n) containing landmark coordinates for the second module 
#' @param method Method to estimate morphological integration; see below for details
#' @param iter Number of iterations for significance testing
#' @export
#' @keywords morphol.integr
#' @author Dean Adams
#' @return Function returns a list with the following components: 
#'   \item{value}{The estimate of morphological integration: PLS.corr or RV}
#'   \item{pvalue}{The significance level of the observed signal}
#'   \item{Xscores}{PLS scores for the first block of landmarks (PLS method only)}
#'   \item{Yscores}{PLS scores for the second block of landmarks (PLS method only)}
#' @references  Bookstein, F. L., P. Gunz, P. Mitteroecker, H. Prossinger, K. Schaefer, and H. Seidler. 
#'   2003. Cranial integration in Homo: singular warps analysis of the midsagittal plane in ontogeny and 
#'   evolution. J. Hum. Evol. 44:167-187.
#' @references Klingenberg, C. P. 2009. Morphometric integration and modularity in configurations of 
#'   landmarks: tools for evaluating a priori hypotheses. Evol. Develop. 11:405-421.
#' @examples
#' data(plethodon) 
#' Y.gpa<-gpagen(plethodon$land)    #GPA-alignment    
#'
#' #Morphological integration using PLS between two modules within a structure
#' morphol.integr(Y.gpa$coords[1:5,,],Y.gpa$coords[6:12,,],method="PLS",iter=99)
#'
#' #Morphological integration using RV between two modules within a structure
#' morphol.integr(Y.gpa$coords[1:5,,],Y.gpa$coords[6:12,,],method="RV",iter=99)
morphol.integr<-function(A,A2,method=c("PLS","RV"),iter=999){
  method <- match.arg(method)
  if (length(dim(A))!=3){
    stop("Data matrix 1 not a 3D array (see 'arrayspecs').")  }
  if(length(grep("-999",A))!=0){
    stop("Data matrix 1 contains missing values. Estimate these first(see 'estimate.missing').")  }
  if (length(dim(A2))!=3){
    stop("Data matrix 2 not a 3D array (see 'arrayspecs').")  }
  if(length(grep("-999",A2))!=0){
    stop("Data matrix 2 contains missing values. Estimate these first(see 'estimate.missing').")  }
  if(is.null(dimnames(A)[[3]])){
    print("No specimen names in data matrix 1. Assuming specimens in same order.")  }
  if(is.null(dimnames(A2)[[3]])){
    print("No specimen names in data matrix 2. Assuming specimens in same order.")  }
  x<-two.d.array(A)
  y<-two.d.array(A2)
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
  XScores<-x%*%U[,1]; YScores<-y%*%V[,1]
  PLS.obs<-cor(XScores,YScores)
  RV.obs<- sum(diag(S12%*%S21))/sqrt(sum(diag(S11%*%S11))*sum(diag(S22%*%S22))) 
  integ.obs<-ifelse(method=="PLS",PLS.obs,RV.obs)
  P.val<-1
  integ.val<-rep(0,iter)
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
    RV.r<- sum(diag(S12.r%*%S21.r))/sqrt(sum(diag(S11.r%*%S11.r))*sum(diag(S22.r%*%S22.r))) 
    integ.r<-ifelse(method=="PLS",PLS.r,RV.r)
    integ.val[i]<-integ.r
    P.val<-ifelse(integ.r>=integ.obs, P.val+1,P.val) 
  }  
  integ.val[iter+1]=integ.obs
  P.val<-P.val/(iter+1)
  if(method=="PLS"){
    par(mfrow=c(1,2))
    hist(integ.val,30,freq=TRUE,col="gray",xlab="PLS Correlation")
      arrows(integ.obs,50,integ.obs,5,length=0.1,lwd=2)
    plot(XScores,YScores,pch=21,bg="black",asp=1,main="PLS Plot")
    layout(1)
    return(list(x.scores=XScores,y.scores=YScores,PLS.corr=integ.obs,pvalue=P.val))
  }
  if(method=="RV"){
    hist(integ.val,30,freq=TRUE,col="gray",xlab="RV Coefficient")
    arrows(integ.obs,50,integ.obs,5,length=0.1,lwd=2)
    return(list(RV=integ.obs,pvalue=P.val))
  }
}

#' Compare modular signal to alternative landmark subsets
#'
#' Function quantifies the degree of morphological integration between two or more modules of Procrustes-aligned 
#'   landmark coordinates and compares this to patterns found by randomly assigning landmarks into subsets
#'
#' The function quantifies the degree of morphological integration between two or more modules of shape data as 
#'   defined by landmark coordinates, and compares this to modular signals found by randomly assigning landmarks 
#'   to modules. It is assumed that the landmarks have previously been aligned using Generalized 
#'   Procrustes Analysis (GPA) [e.g., with \code{\link{gpagen}}]. The degree of morphological integration 
#'   is quantified using the RV coefficient (Klingenberg 2009). If more than two modules are defined, the average
#'   RV coefficient is utilized (see Klingenberg 2009). The RV coefficient for the observed modular 
#'   hypothesis is then compared to a distribution of values obtained by randomly assigning landmarks into 
#'   subsets, with the restriction that the number of landmarks in each subset is identical to that observed 
#'   in each of the original partitions. A significant modular signal is found when the observed RV coefficient 
#'   is small relative to this distribution (see Klingenberg 2009). A histogram of coefficients obtained via 
#'   resampling is presented, with the observed value designated by an arrow in the plot. 
#'
#' @param A An array (p x k x n) containing GPA-aligned coordinates for all specimens
#' @param landgroups A list of which landmarks belong in which partition (e.g. A,A,A,B,B,B,C,C,C)
#' @param iter Number of iterations for significance testing
#' @export
#' @keywords compare.modular.partitions
#' @author Dean Adams
#' @return Function returns a list with the following components: 
#'   \item{RV}{The estimate of morphological integration}
#'   \item{pvalue}{The significance level of the observed signal}
#'   \item{RV.min}{The minimal RV coefficient found via landmark permutation}
#'   \item{RV.min.partitions}{A list of landmarks assigned to partitions that yields the minimal RV coefficient}
#' @references Klingenberg, C. P. 2009. Morphometric integration and modularity in configurations of 
#'   landmarks: tools for evaluating a priori hypotheses. Evol. Develop. 11:405-421.
#' @examples
#' data(plethodon) 
#' Y.gpa<-gpagen(plethodon$land)    #GPA-alignment    
#'  #landmarks on the skull and mandible assigned to partitions
#' land.gps<-c("A","A","A","A","A","B","B","B","B","B","B","B") 
#'
#' compare.modular.partitions(Y.gpa$coords,land.gps,iter=99)
#' #Result implies that the skull and mandible are not independent modules
compare.modular.partitions<-function(A,landgroups,iter=999){
  if (length(dim(A))!=3){
    stop("Data matrix 1 not a 3D array (see 'arrayspecs').")  }
  if(length(grep("-999",A))!=0){
    stop("Data matrix 1 contains missing values. Estimate these first(see 'estimate.missing').")}
  p<-dim(A)[1]; k<-dim(A)[2]
  landgroups<-as.factor(landgroups)
  if(length(landgroups)!=p){stop("Not all landmarks are assigned to a partition.")}
  ngps<-nlevels(landgroups)
  gps<-as.factor(rep(landgroups,k,each = k, length=p*k))
  x<-two.d.array(A)
  S<-cov(x)
  RV.gp<-array(0,dim=c(ngps,ngps))
  for (i in 1:(ngps-1)){
    for (j in 2:ngps){
      S11<-S[which(gps==levels(gps)[i]),which(gps==levels(gps)[i])]
      S22<-S[which(gps==levels(gps)[j]),which(gps==levels(gps)[j])]
      S12<-S[which(gps==levels(gps)[i]),which(gps==levels(gps)[j])]
      S21<-t(S12)
      RV.gp[i,j]<- sum(diag(S12%*%S21))/sqrt(sum(diag(S11%*%S11))*sum(diag(S22%*%S22)))
      diag(RV.gp)<-0
    }
  }
  RV.obs<-sum(RV.gp)/(ngps/2*(ngps-1))
  RV.min<-RV.obs; partition.min<-landgroups
  P.val<-1
  RV.val<-rep(0,iter)
  for(ii in 1:iter){
    landgroups.r<-sample(landgroups)
    gps.r<-as.factor(rep(landgroups.r,k,each = k, length=p*k))    
    RV.gp.r<-array(0,dim=c(ngps,ngps))
    for (i in 1:(ngps-1)){
      for (j in 2:ngps){
        S11.r<-S[which(gps.r==levels(gps.r)[i]),which(gps.r==levels(gps.r)[i])]
        S22.r<-S[which(gps.r==levels(gps.r)[j]),which(gps.r==levels(gps.r)[j])]
        S12.r<-S[which(gps.r==levels(gps.r)[i]),which(gps.r==levels(gps.r)[j])]
        S21.r<-t(S12.r)
        RV.gp.r[i,j]<- sum(diag(S12.r%*%S21.r))/sqrt(sum(diag(S11.r%*%S11.r))*sum(diag(S22.r%*%S22.r)))
        diag(RV.gp.r)<-0
      }
    }
    RV.r<-sum(RV.gp.r)/(ngps/2*(ngps-1))
    RV.val[ii]<-RV.r
    if (RV.r< RV.min) {RV.min<-RV.r; partition.min<-landgroups.r}
    P.val<-ifelse(RV.r<=RV.obs, P.val+1,P.val) 
  }
  RV.val[iter+1]=RV.obs
  P.val<-P.val/(iter+1)
  hist(RV.val,30,freq=TRUE,col="gray",xlab="RV Coefficient")
  arrows(RV.obs,50,RV.obs,5,length=0.1,lwd=2)
  return(list(RV=RV.obs,pvalue=P.val,RV.min=RV.min,RV.min.partitions=partition.min))
}

#' Assessing phylogenetic signal in morphometric data
#'
#' Function calculates the degree of phylogenetic signal from a set of Procrustes-aligned specimens
#'
#' The function estimates the degree of phylogenetic signal present in shape data for a given phylogeny based 
#'   on a Brownian motion model of evolution. It is assumed that the landmarks have previously been aligned 
#'   using Generalized Procrustes Analysis (GPA) [e.g., with \code{\link{gpagen}}]. To assess phylogenetic 
#'   signal under alternative models, first consider branch-length transformations of the phylogeny.  
#'   Phylogenetic signal is estimated  as the sum of squared changes (SSC) in 
#'   shape along all branches of the phylogeny (Klingenberg and Gidaszewski 2010). Significance testing 
#'   is found by permuting the shape data among the tips of the phylogeny. A plot of the specimens in tangent 
#'   space with the phylogeny superimposed is included, and ancestral states are returned. Note that the 
#'   method can be slow as ancestral states must be estimated for every iteration.
#'
#' @param phy A phylogenetic tree of type 'phylo'
#' @param A An array (p x k x n) containing GPA-aligned coordinates for a set of specimens
#' @param iter Number of iterations for significance testing
#' @keywords physignal
#' @author Dean Adams
#' @export
#' @return Function returns a list with the following components: 
#'   \item{phy.signal}{The estimate of phylogenetic signal}
#'   \item{pvalue}{The significance level of the observed signal}
#' @references Klingenberg, C. P., and N. A. Gidaszewski. 2010. Testing and quantifying phylogenetic signals 
#'   and homoplasy in morphometric data. Syst. Biol. 59:245-261.
#' @examples
#' data(plethspecies) 
#' Y.gpa<-gpagen(plethspecies$land)    #GPA-alignment    
#'
#' physignal(plethspecies$phy,Y.gpa$coords,iter=5)
physignal<-function(phy,A,iter=249){
  if (length(dim(A))!=3){
    stop("Data matrix not a 3D array (see 'arrayspecs').")  }
  if(length(grep("-999",A))!=0){
    stop("Data matrix contains missing values. Estimate these first(see 'estimate.missing').")  }
  if (class(phy) != "phylo") 
    stop("tree must be of class 'phylo.'")
  if (!is.binary.tree(phy)) 
    stop("tree is not fully bifurcating.")
  if(is.null(dimnames(A)[[3]])){
    stop("Data matrix does not include taxa names as dimnames for 3rd dimension.")  }
  N<-length(phy$tip.label)
  x<-two.d.array(A)
  if(N!=dim(x)[1]){
    stop("Number of taxa in data matrix and tree are not not equal.")  }
  if(length(match(rownames(x), phy$tip.label))!=N) 
    stop("Data matrix missing some taxa present on the tree.")
  if(length(match(phy$tip.label,rownames(x)))!=N) 
    stop("Tree missing some taxa in the data matrix.")
  x<-x[phy$tip.label, ]  
  SSC.o<-NULL
  anc.states<-matrix(NA, nrow=(nrow(x)-1), ncol=ncol(x))
  for (i in 1:ncol(x)){
    options(warn=-1)  
    anc.states[,i]<-as.vector(ace(x[,i], compute.brlen(phy,1), type="continuous", method="ML")$ace)}
  dist.mat<-as.matrix(dist(rbind(as.matrix(x),as.matrix(anc.states)))^2)   
  SSC.o<-0
  for (i in 1:nrow(phy$edge)){
    SSC.o<-SSC.o+dist.mat[phy$edge[i,1],phy$edge[i,2]]    }
  P.val<-1
  SSC.val<-rep(0,iter)
  for(ii in 1:iter){
    x.r<-x[sample(nrow(x)),]	
    row.names(x.r)<-row.names(x)
    SSC.r<-NULL
    anc.states.r<-matrix(NA, nrow=(nrow(x)-1), ncol=ncol(x))
    for (i in 1:ncol(x.r)){
      options(warn=-1)  
      anc.states.r[,i]<-as.vector(ace(x.r[,i], compute.brlen(phy,1), type="continuous", method="ML")$ace)}
    dist.mat.r<-as.matrix(dist(rbind(as.matrix(x.r),as.matrix(anc.states.r)))^2)   
    SSC.r<-0
    for (i in 1:nrow(phy$edge)){
      SSC.r<-SSC.r+dist.mat.r[phy$edge[i,1],phy$edge[i,2]]    }
    P.val<-ifelse(SSC.r<=SSC.o, P.val+1,P.val) 
    SSC.val[ii]<-SSC.r
  }  
  P.val<-P.val/(iter+1)
  SSC.val[iter+1]=SSC.o
  hist(SSC.val,30,freq=TRUE,col="gray",xlab="Phylogenetic Signal")
  arrows(SSC.o,50,SSC.o,5,length=0.1,lwd=2)
  plotGMPhyloMorphoSpace(phy,A,ancStates=FALSE)
  return(list(phy.signal=SSC.o,pvalue=P.val,anc.states=anc.states))
}

#' Comparing rates of shape evolution on phylogenies
#'
#' Function calculates rates of shape evolution for two or more groups of species on a phylogeny from a set of Procrustes-aligned specimens
#'
#' The function compares rates of morphological evolution for two or more groups of species on a phylogeny, under a 
#'  Brownian motion model of evolution. It is assumed that the landmarks have previously been aligned 
#'  using Generalized Procrustes Analysis (GPA) [e.g., with \code{\link{gpagen}}]. The approach is based on the distances
#'  between species in morphospace after phylogenetic transformation (Adams 2014). From the data the rate of shape evolution
#'  for each group is calculated, and a ratio of rates is obtained. If three or more groups of species are used, the ratio of 
#'  the maximum to minimum rate is used as a test statistic (see Adams 2014). Significance testing 
#'  is accomplished by phylogenetic simulation in which tips data are obtained under Brownian motion using a single 
#'  evolutionary rate for all species on the phylogeny. A histogram of evolutionary rate ratios obtained via phylogenetic simulation is presented, 
#'  with the observed value designated by an arrow in the plot. 
#'
#' @param phy A phylogenetic tree of type 'phylo'
#' @param A An array (p x k x n) containing GPA-aligned coordinates for a set of specimens
#' @param gp A factor array designating group membership
#' @param iter Number of iterations for significance testing
#' @keywords compare.phylo.shape.rate
#' @author Dean Adams
#' @export
#' @return Function returns a list with the following components: 
#'   \item{sigmad.all}{The phylogenetic evolutionary rate for all species on the phylogeny}
#'   \item{sigmad.ratio}{The ratio of maximum to minimum evolutionary rates}
#'   \item{pvalue}{The significance level of the observed ratio}
#' @references Adams, D.C. 2014. Quantifying and comparing phylogenetic evolutionary rates for 
#'  shape and other high-dimensional phenotypic data. Syst. Biol. (Accepted).
#' @examples
#' data(plethspecies) 
#' Y.gpa<-gpagen(plethspecies$land)    #GPA-alignment    
#'
#'  gp.end<-factor(c(0,0,1,0,0,1,1,0,0))  #endangered species vs. rest
#'  names(gp.end)<-plethspecies$phy$tip
#' 
#' compare.evol.rates(plethspecies$phy,Y.gpa$coords,gp=gp.end,iter=49)
compare.evol.rates<-function(phy,A,gp,iter=999 ){
  if (length(dim(A))!=3){
    stop("Data matrix not a 3D array (see 'arrayspecs').")  }
  if(length(grep("-999",A))!=0){
    stop("Data matrix contains missing values. Estimate these first(see 'estimate.missing').")  }
  if(!is.factor(gp)){
    stop("gp is not a factor.")}
  if (is.null(names(gp))){
    stop("Factor contains no names. Use names() to assign specimen names to group factor.")}
  if (class(phy) != "phylo") 
    stop("tree must be of class 'phylo.'")
  if(is.null(dimnames(A)[[3]])){
    stop("Data matrix does not include taxa names as dimnames for 3rd dimension.")  }
  ntaxa<-length(phy$tip.label)
  x<-two.d.array(A)
  N<-nrow(x)  
  if(N!=dim(x)[1]){
    stop("Number of taxa in data matrix and tree are not not equal.")  }
  if(length(match(rownames(x), phy$tip.label))!=N) 
    stop("Data matrix missing some taxa present on the tree.")
  if(length(match(phy$tip.label,rownames(x)))!=N) 
    stop("Tree missing some taxa in the data matrix.")
  p<-ncol(x)           
  sigma.d<-function(phy,x,N,gp){
    x<-prcomp(x)$x
    gp<-gp[rownames(x)]
    ngps<-nlevels(gp)
    gpsz<-table(gp)
    ones<-array(1,N)
    C<-vcv.phylo(phy)
    C<-C[rownames(x),rownames(x)]
    a.obs<-colSums(solve(C))%*%x/sum(solve(C))  
    eigC <- eigen(C)
    D.mat<-solve(eigC$vectors %*% diag(sqrt(eigC$values)) %*% t(eigC$vectors)) 
    dist.adj<-as.matrix(dist(rbind((D.mat%*%(x-(ones%*%a.obs))),0))) 
    vec.d2<-dist.adj[N+1,1:N]^2
    sigma.d<-tapply(vec.d2,gp,sum)/gpsz/p   
    sigma.d.all<-sum(vec.d2)/N/p            
    if (ngps==2){sigma.d.rat<-max(sigma.d)/min(sigma.d)}
    if(ngps>2){
      sigma.d.rat.gp<- array(0, dim = c(ngps, ngps))
      for (i in 1:(ngps - 1)) {
        for (j in 2:ngps) {
          tmp<-c(sigma.d[i],sigma.d[j])
          sigma.d.rat.gp[i, j] <-max(tmp)/min(tmp)
            diag(sigma.d.rat.gp) <- 0
        }
      }
      sigma.d.rat <- sum(sigma.d.rat.gp)/(ngps/2 * (ngps - 1))
    }
    return(list(sigma.all = sigma.d.all, ratio = sigma.d.rat,sigma.d.all=sigma.d))
  }
  sigmad.obs<-sigma.d(phy,x,ntaxa,gp)  
  rate.mat<-diag(sigmad.obs$sigma.all,p)
  x.sim<-sim.char(phy,rate.mat,nsim=iter) 
  sig.sim<-1
  rate.val<-rep(0,iter)
  for(ii in 1:iter){
      sigmad.sim<-sigma.d(phy,x.sim[,,ii],ntaxa,gp)
      sig.sim<-ifelse(sigmad.sim$ratio>=sigmad.obs$ratio, sig.sim+1,sig.sim)
      rate.val[ii]<-sigmad.sim$ratio
  }
  sig.sim<-sig.sim/(iter+1)
  rate.val[iter+1]=sigmad.obs$ratio
  hist(rate.val,30,freq=TRUE,col="gray",xlab="SigmaD ratio")
  arrows(sigmad.obs$ratio,50,sigmad.obs$ratio,5,length=0.1,lwd=2)
  return(list(sigmad.all = sigmad.obs$sigma.d.all, sigmad.ratio = sigmad.obs$ratio, pvalue=sig.sim))
}

#' Procrustes ANOVA/regression for shape data
#'
#' Function performs Procrustes ANOVA with permutation procedures to assess statistical hypotheses describing 
#'   patterns of shape variation and covariation for a set of Procrustes-aligned coordinates
#'
#' The function quantifies the relative amount of shape variation attributable to one or more factors in a 
#'   linear model and assesses this variation via permutation. Data input is specified by a formula (e.g., 
#'   y~X), where 'y' specifies the response variables (shape data), and 'X' contains one or more independent 
#'   variables (discrete or continuous). The response matrix 'y' must be in the form of a two-dimensional data 
#'   matrix of dimension (n x [p x k]), rather than a 3D array.  It is assumed that the landmarks have previously 
#'   been aligned using Generalized Procrustes Analysis (GPA) [e.g., with \code{\link{gpagen}}]. The function
#'   \code{\link{two.d.array}} can be used to obtain a two-dimensional data matrix from a 3D array of landmark
#'   coordinates. The names specified for the independent (x) variables in the formula represent one or more 
#'   vectors containing continuous data or factors. It is assumed that the order of the specimens in the 
#'   shape matrix matches the order of values in the independent variables.
#'
#'   The function performs statistical assessment of the terms in the model using Procrustes distances among 
#'   specimens, rather than explained covariance matrices among variables. With this approach, the sum-of-squared 
#'   Procrustes distances are used as a measure of SS (see Goodall 1991). The observed SS are evaluated through 
#'   permutation, where the rows of the shape matrix are randomized relative to the design matrix. Procedurally, 
#'   Procrustes ANOVA is identical to permutational-MANOVA as used in other fields (Anderson 2001). For several 
#'   reasons, Procrustes ANOVA is particularly useful for shape data. First, covariance matrices from GPA-aligned 
#'   Procrustes coordinates are singular, and thus standard approaches such as MANOVA cannot be accomplished 
#'   unless generalized inverses are utilized. This problem is accentuated when using sliding semilandmarks. 
#'   Additionally, GM datasets often have more variables than specimens (the 'small N large P' problem). In 
#'   these cases, distance-based procedures can still be utilized to assess statistical hypotheses, whereas 
#'   standard linear models cannot. 
#'
#' @param f1 A formula for the linear model (e.g., y~x1+x2)
#' @param data An optional value specifying a data frame containing all data (not required)
#' @param iter Number of iterations for significance testing
#' @keywords procD.lm
#' @export
#' @author Dean Adams
#' @return Function returns an ANOVA table of statistical results for all factors: df (for each factor), SS, MS, Prand.
#' @references Anderson MJ. 2001. A new method for non-parametric multivariate analysis of variance. 
#'    Austral Ecology 26: 32'46.
#' @references Goodall, C. R. 1991. Procrustes methods in the statistical analysis of shape. Journal of the 
#'    Royal Statistical Society B 53:285-339.
#' @examples
#' ### MANOVA example for Goodall's F test (multivariate shape vs. factors)
#' data(plethodon) 
#' Y.gpa<-gpagen(plethodon$land)    #GPA-alignment    
#' y<-two.d.array(Y.gpa$coords)
#'
#' procD.lm(y~plethodon$species*plethodon$site,iter=99)
#'
#' ### Regression example
#' data(ratland)
#' rat.gpa<-gpagen(ratland)         #GPA-alignment
#'
#' procD.lm(two.d.array(rat.gpa$coords)~rat.gpa$Csize,iter=99)
procD.lm<-function(f1,data=NULL,iter=999){
  form.in<-formula(f1)
  Terms<-terms(form.in,keep.order=TRUE)
  Y<-eval(form.in[[2]],parent.frame())
  if (length(dim(Y))!=2){
    stop("Response matrix (shape) not a 2D array. Use 'two.d.array' first.")  }
  if(length(grep("-999",Y))!=0){
    stop("Response data matrix (shape) contains missing values. Estimate these first(see 'estimate.missing').")  }
  if(is.null(dimnames(Y)[[1]])){
    print("No specimen names in response matrix. Assuming specimens in same order.")  }
  df<-df.tmp<-SS.tmp<-SS.obs<-array()
  dat<-model.frame(form.in,data)
  for (i in 1:ncol(attr(Terms, "factors"))){
    mod.mat<-model.matrix(Terms[1:i],data=dat)
    SS.tmp[i]<-sum(dist(predict(lm(Y~mod.mat)))^2)/(nrow(Y))
    ifelse(i==1, SS.obs[i]<-SS.tmp[i], SS.obs[i]<-(SS.tmp[i]-SS.tmp[i-1]))
    df.tmp[i]<-ifelse(ncol(mod.mat)==1,1,(ncol(mod.mat)-1))
    ifelse(i==1, df[i]<-df.tmp[i], df[i]<-(df.tmp[i]-df.tmp[i-1]))
  }
  MS<-SS.obs/df
  mod.mat<-model.matrix(Terms[1])
  SS.tot<-sum(dist(predict(lm(Y~mod.mat)))^2)/(nrow(Y))+
    sum(dist(resid(lm(Y~mod.mat)))^2)/(nrow(Y))
  df.tot<-nrow(Y)-1
  MS.tot<-SS.tot/df.tot
  P.val<-array(1,dim=length(SS.obs))
  for(i in 1:iter){
    SS.tmp<-SS.r<-array()
    Y.r<-Y[sample(nrow(Y)),]	
    for (ii in 1:ncol(attr(Terms, "factors"))){
      mod.mat<-model.matrix(Terms[1:ii])
      SS.tmp[ii]<-sum(dist(predict(lm(Y.r~mod.mat)))^2)/(nrow(Y))
      ifelse(ii==1, SS.r[ii]<-SS.tmp[ii], SS.r[ii]<-(SS.tmp[ii]-SS.tmp[ii-1]))
    }
    P.val<-ifelse(SS.r>=SS.obs, P.val+1,P.val) 
  }
  P.val<-P.val/(iter+1)
  anova.tab<-cbind(df,SS.obs,MS,P.val)
  anova.tab<-rbind(anova.tab,c(df.tot,SS.tot,MS.tot,NA))
  rownames(anova.tab)<-c(colnames(attr(Terms, "factors")), "Total")
  return(anova.tab)
}

#' Pairwise Group Comparisons
#'
#' Function performs pairwise comparisons among groups using the Euclidean distances among group means.
#'
#' The function performs pairwise comparisons to identify shape among groups. The function is designed as a post-hoc
#'  test to Procrustes ANOVA, where the latter has identified significant shape variation explained by a grouping factor. 
#'  The function takes as input the shape data (y), and a grouping factor (x). It then estimates the Euclidean distances 
#'  among group means, which are used as test values. These are then statistically evaluated through permutation, 
#'  where the rows of the shape matrix are randomized relative to the grouping variable. 
#'
#' The input for the shape data (y) must be in the form of a two-dimensional data matrix of dimension (n x [p x k]), 
#'   rather than a 3D array.  It is assumed that the landmarks have previously 
#'   been aligned using Generalized Procrustes Analysis (GPA) [e.g., with \code{\link{gpagen}}]. The function
#'   \code{\link{two.d.array}} can be used to obtain a two-dimensional data matrix from a 3D array of landmark
#'   coordinates. 
#'
#' @param y A two-dimensional array of shape data
#' @param x A factor defining groups
#' @param iter Number of iterations for permutation test
#' @keywords PairwiseD.test
#' @export
#' @author Dean Adams
#' @return Function returns a list with the following components: 
#'   \item{Dist.obs}{A matrix of Euclidean distances among group means}
#'   \item{Prob.Dist}{A matrix of pairwise significance levels based on permutation}
#' @examples
#' data(plethodon) 
#' Y.gpa<-gpagen(plethodon$land)    #GPA-alignment    
#' y<-two.d.array(Y.gpa$coords)
#' ### Procrustes ANOVA
#' procD.lm(y~plethodon$species,iter=99)
#'
#' ### Pairwise comparisons
#' pairwiseD.test(y,plethodon$species,iter=99)
pairwiseD.test<-function(y,x,iter=999){
  if (length(dim(y))!=2){
    stop("Response matrix (shape) not a 2D array. Use 'two.d.array' first.")  } 
  if(length(grep("-999",y))!=0){
    stop("Response data matrix (shape) contains missing values. Estimate these first(see 'estimate.missing').")  }
  if(!is.factor(x)){
    stop("X-variable is not a factor.")}
  yhat.obs<-predict(lm(y~x))
  lsmeans.obs <- rowsum(yhat.obs, x)/as.vector(table(x))     
  D.obs<-as.matrix(dist(lsmeans.obs))   
  PDist<-array(1, dim=c(dim(lsmeans.obs)[1],dim(lsmeans.obs)[1]))
  for(i in 1:iter){
    y.r<-y[sample(nrow(y)),] 
    yhat.rand<-predict(lm(y.r~x))
    lsmeans.rand <- rowsum(yhat.rand, x)/as.vector(table(x))     
    D.rand<-as.matrix(dist(lsmeans.rand))   
    PDist<-ifelse(D.rand>=D.obs, PDist+1, PDist)
  }
  PDist<-PDist/(iter+1)
  return(list(Dist.obs=D.obs,Prob.Dist=PDist))
}

#' Analysis of bilateral symmetry
#'
#' Function performs an analysis of directional and fluctuating asymmetry for bilaterally symmetric objects 
#'
#' The function quantifies components of shape variation for a set of specimens as described by their patterns of symmetry
#'  and asymmetry. Here, shape variation is decomposed into variation among individuals, variation among sides (directional 
#'  asymmetry), and variation due to an individual x side interaction (fluctuating symmetry). These components are then 
#'  statistically evaluated using Procrustes ANOVA and Goodall's F tests (i.e. an isotropic model of shape variation). Methods for both 
#'  matching symmetry and object symmetry can be implemented. Matching symmetry is when each object contains mirrored 
#'  pairs of structures (e.g., right and left hands) while object symmetry is when a single object is symmetric 
#'  about a midline (e.g., right and left sides of human faces). Analytical and computational details concerning the 
#'  analysis of symmetry in geometric morphometrics can be found in Mardia et al. 2000; Klingenberg et al. 2002.
#'
#' Analyses of symmetry for matched pairs of objects is implemented when "object.sym=FALSE". Here, a 3D array [p x k x 2n] 
#'  contains the landmark coordinates for all pairs of structures (2 structures for each of n specimens). Because the two sets of 
#'  structures are on opposite sides, they represent mirror images, and one set must be reflected prior to the analysis to 
#'  allow landmark correspondence. IT IS ASSUMED THAT THE USER HAS DONE THIS PRIOR TO PERFORMING THE SYMMETRY ANALYSIS. 
#'  Reflecting a set of specimens may be accomplished by multiplying one coordinate dimension 
#'  by '-1' for these structures (either the x-, the y-, or the z-dimension). A vector containing information on individuals 
#'  and sides must also be supplied. Replicates of each specimen may also be included in the dataset, and when specified will be 
#'  used as measurement error (see Klingenberg and McIntyre 1998). 
#' 
#' Analyses of object symmetry is implemented when "object.sym=TRUE". Here, a 3D array [p x k x n] contains the landmark 
#'  coordinates for all n specimens. To obtain information about asymmetry, the function generates a second set of objects 
#'  by reflecting them about one of their coordinate axes. The landmarks across the line of symmetry are then relabeled to obtain
#'  landmark correspondence. The user must supply a list of landmark pairs. A vector containing information on individuals 
#'  must also be supplied. Replicates of each specimen may also be included in the dataset, and when specified will be 
#'  used as measurement error. 
#'
#' @param A An array (p x k x n) containing GPA-aligned coordinates for a set of specimens [for "object.sym=FALSE, A is of dimension (n x k x 2n)]
#' @param ind A vector containing labels for each individual. For matching symmetry, the matched pairs receive the same 
#' label (replicates also receive the same label).
#' @param side An optional vector (for matching symmetry) designating which object belongs to which 'side-group'
#' @param replicate An optional vector designating which objects belong to which group of replicates
#' @param object.sym A logical value specifying whether the analysis should proceed based on object symmetry (TRUE) or matching symmetry (FALSE)
#' @param land.pairs An optional matrix (for object symmetry) containing numbers for matched pairs of landmarks across the line of symmetry 
#' @param GPAoutput A logical value indicating whether results from GPA should be returned
#' @param plot.method A logical value indicating which type of plot should be used for the directional and fluctuating components
#' of asymmetry (see \code{\link{plotRefToTarget}})
#' @param mag The desired magnification to be used when visualizing the shape differences for the directional and fluctuating
#' components of shape variation (see \code{\link{plotRefToTarget}})
#' @param ... Additional parameters to be passed to 'gpagen'
#' @keywords bilat.symmetry
#' @export
#' @author Dean Adams
#' @return Function returns a list with the following components: 
#'   \item{aligned}{Procrustes aligned coordinates (when GPAoutput=TRUE)}
#'   \item{csize}{Centroid size of aligned specimens (when GPAoutput=TRUE)}
#'   \item{symm.shape}{The symmetric component of shape variation of the aligned specimens}
#'   \item{asymm.shape}{The asymmetric component of shape variation of the aligned specimens}
#'   \item{ANOVA.shape}{Procrustes ANOVA table assessing patterns of shape asymmetry}
#'   \item{ANOVA.size}{Procrustes ANOVA table assessing patterns of shape asymmetry (when object.sym=FALSE)}
#' @references Klingenberg, C.P. and G.S. McIntyre. 1998. Quantitative genetics of geometric shape in the mouse mandible. Evolution. 55:2342-2352.
#' @references Mardia, K.V., F.L. Bookstein, and I.J. Moreton. 2000. Statistical assessment of bilateral symmetry of shapes. Biometrika. 87:285-300.
#' @references Klingenberg, C.P., M. Barluenga, and A. Meyer. 2002. Shape analysis of symmetric structures: quantifying variation among
#' individuals and asymmetry. Evolution. 56:1909-1920.
#' @examples
#' #Example of matching symmetry
#'
#' data(mosquito)
#' bilat.symmetry(mosquito$wingshape,ind=mosquito$ind,side=mosquito$side,
#' replicate=mosquito$replicate,object.sym=FALSE)
#'
#' #Example of object symmetry
#'
#' data(scallops)
#' bilat.symmetry(scallops$coorddata,ind=scallops$ind,object.sym=TRUE,land.pairs=scallops$land.pairs)
bilat.symmetry<-function(A,ind=NULL,side=NULL,replicate=NULL,object.sym=FALSE,land.pairs=NULL,GPAoutput=FALSE,plot.method="TPS",mag=1.0,...){
  if (length(dim(A))!=3){
    stop("Data matrix 1 not a 3D array (see 'arrayspecs').")  }
  if(length(grep("-999",A))!=0){
    stop("Data matrix 1 contains missing values. Estimate these first(see 'estimate.missing').")  }
  if(is.null(ind)){stop("Individuals not specified.")}
  ind<-as.factor(ind)
  n<-dim(A)[3];   k<-dim(A)[2];  p<-dim(A)[1]; shpsp<-k*p-k-k*(k-1)/2-1; nind<-nlevels(ind) 
  if(!is.null(replicate)){replicate<-as.factor(replicate); nrep<-nlevels(replicate) }
  if(object.sym==FALSE){
    if(is.null(side)){stop("Sides not specified.")}  
    side<-as.factor(side)
    gpa.res<-gpagen(A,ShowPlot=FALSE)
    shape<-two.d.array(gpa.res$coords)
    f1<-"shape~ind*side"; if(!is.null(replicate)){f1<-paste(f1,"ind:side:replicate",sep="+")}
    f1<-as.formula(f1)
    f2<-"gpa.res$Csize~ind*side"; if(!is.null(replicate)){f2<-paste(f2,"ind:side:replicate",sep="+")}
    f2<-as.formula(f2)
    res.shape<-procD.lm(f1,iter=1)
    res.shape[1,1]<-(nind-1)*shpsp;res.shape[2,1]<-shpsp; res.shape[3,1]<-(nind-1)*shpsp
    if(!is.null(replicate)){res.shape[4,1]<-(nrep-1)*nind*2*shpsp}
    res.shape[,3]<-res.shape[,2]/res.shape[,1] 
    F<-array(NA,nrow(res.shape)); F[1]<-res.shape[1,3]/res.shape[3,3]
    F[2]<-res.shape[2,3]/res.shape[3,3]; if(!is.null(replicate)){F[3]<-res.shape[3,3]/res.shape[4,3]}
    P<-array(NA,nrow(res.shape)); P[1]<-1-pf(F[1],res.shape[1,1],res.shape[3,1])
    P[2]<-1-pf(F[2],res.shape[2,1],res.shape[3,1]); if(!is.null(replicate)){
      P[3]<-1-pf(F[3],res.shape[3,1],res.shape[4,1])}
    res.shape<-cbind(res.shape[,-4],F,P)
    colnames(res.shape)[4]<-"F.Goodall";colnames(res.shape)[5]<-"P.param"
    res.shape<-res.shape[-nrow(res.shape),]
    res.size<-summary(aov(f2))[[1]]; res.size<-res.size[,(1:3)]   
    F<-array(NA,nrow(res.size)); F[1]<-res.size[1,3]/res.size[3,3]
    F[2]<-res.size[2,3]/res.size[3,3]; if(!is.null(replicate)){F[3]<-res.size[3,3]/res.size[4,3]}
    P<-array(NA,nrow(res.size)); P[1]<-1-pf(F[1],res.size[1,1],res.size[3,1])
    P[2]<-1-pf(F[2],res.size[2,1],res.size[3,1]); if(!is.null(replicate)){
      P[3]<-1-pf(F[3],res.size[3,1],res.size[4,1])}
    res.size<-cbind(res.size,F,P);  if(!is.null(replicate)){res.size<-res.size[(1:4),];
                                                            rownames(res.size)[4]<-"replicate"}
    DA.mns <- arrayspecs((rowsum(predict(lm(shape~side)), 
                                 side)/as.vector(table(side))),p,k,byLand=FALSE)
    MSCP.FA<-summary(manova(lm(f1)))$SS[[3]]/res.shape[3,1]
    eig.FA<-eigen(MSCP.FA); PC1.eigval<-eig.FA$values[1]/sum(eig.FA$values)
    PC1<-shape%*%eig.FA$vec[,1]
    FA.mns<-arrayspecs((rbind(shape[which.min(PC1),],shape[which.max(PC1),])),p,k,byLand=FALSE) 
    mn.shape<-mshape(gpa.res$coords)
    symm.component<-arrayspecs(predict(lm(shape~ind)),p,k,byLand=FALSE)
    asymm.component<-array(data=NA,dim=c(p,k,n))     
    for (i in 1:n){ asymm.component[,,i]<-(gpa.res$coords[,,i]-symm.component[,,i]) + mn.shape}
    if(k==2){  
      #2d plots
      par(mfrow=c(2,2),oma=c(1.5,0,1.5,0))
      plotAllSpecimens(symm.component)
      plotAllSpecimens(asymm.component)
      plotRefToTarget(DA.mns[,,1],DA.mns[,,2],method=plot.method,mag=mag,main="Directional Asymmetry")
      plotRefToTarget(FA.mns[,,1],FA.mns[,,2],method=plot.method,mag=mag,main="Fluctuating Asymmetry")
      mtext("Symmetric Shape Component (left) and Asymmetric Shape Component (right)",outer = TRUE,side=3)
      mtext("Mean directional (left) and fluctuating (right) asymmetry",side = 1, outer = TRUE)
    }
    if (k==3){
     #3d plots
     par(mfrow=c(2,2),oma=c(1.5,0,1.5,0)) 
     plotRefToTarget(DA.mns[,,1],DA.mns[,,2],method=plot.method,mag=mag,main="Directional Asymmetry")
     plotRefToTarget(FA.mns[,,1],FA.mns[,,2],method=plot.method,mag=mag,main="Fluctuating Asymmetry")
     mtext("Mean directional (first plot) and fluctuating (second plot) asymmetry",side = 1, outer = TRUE)
    }    
    if(GPAoutput==FALSE){return(list(symm.shape=symm.component,asymm.shape=asymm.component, 
		ANOVA.size=res.size,ANOVA.Shape=res.shape))}
    if(GPAoutput==TRUE){return(list(aligned=gpa.res$coords,csize=gpa.res$Csize,symm.shape=symm.component,asymm.shape=asymm.component,
		ANOVA.size=res.size,ANOVA.Shape=res.shape))}    
  }
  if(object.sym==TRUE){
    npairs<-nrow(land.pairs); nl<-p-2*npairs
    A2<-A; 
    if(is.null(land.pairs)){stop("Landmark pairs not specified.")}  
    for (i in 1:n){
      for (j in 1:nrow(land.pairs)){
        A2[land.pairs[j,1],,i]<-A[land.pairs[j,2],,i]
        A2[land.pairs[j,2],,i]<-A[land.pairs[j,1],,i]            
      }
    }
    A2[,1,]<-A2[,1,]*-1
    A<-array(c(A,A2), c(p,k, 2*n))
    ind<-rep(ind,2);side<-gl(2,n); if(!is.null(replicate)){replicate<-rep(replicate,2)}
    gpa.res<-gpagen(A,ShowPlot = FALSE)
    shape<-two.d.array(gpa.res$coords)    
    f1<-"shape~ind*side"; if(!is.null(replicate)){f1<-paste(f1,"ind:side:replicate",sep="+")}
    f1<-as.formula(f1)
    res.shape<-procD.lm(f1,iter=1)
    res.shape[,2]<-res.shape[,2]/2 
    res.shape[2,1]<-ifelse(k==2,((2*npairs+nl-2)),((3*npairs+nl-3)))
    res.shape[1,1]<-res.shape[3,1]<-(nind-1)*res.shape[2,1]
    if(k==3){res.shape[1,1]=res.shape[1,1]+((nind-1)*(nl-1))}
    if(!is.null(replicate)){res.shape[4,1]<-(nrep-1)*nind*shpsp}
    res.shape[,3]<-res.shape[,2]/res.shape[,1]  
    F<-array(NA,nrow(res.shape)); F[1]<-res.shape[1,3]/res.shape[3,3]
    F[2]<-res.shape[2,3]/res.shape[3,3]; if(!is.null(replicate)){F[3]<-res.shape[3,3]/res.shape[4,3]}
    P<-array(NA,nrow(res.shape)); P[1]<-1-pf(F[1],res.shape[1,1],res.shape[3,1])
    P[2]<-1-pf(F[2],res.shape[2,1],res.shape[3,1]); if(!is.null(replicate)){
      P[3]<-1-pf(F[3],res.shape[3,1],res.shape[4,1])}
    res.shape<-cbind(res.shape[,-4],F,P)
    res.shape<-res.shape[-nrow(res.shape),]
    colnames(res.shape)[4]<-"F.Goodall";colnames(res.shape)[5]<-"P.param"
    if(!is.null(replicate)){rownames(res.shape)[4]<-"replicate=error"}  
    DA.mns <- arrayspecs((rowsum(predict(lm(shape~side)), 
                                 side)/as.vector(table(side))),p,k,byLand=FALSE)
    MSCP.FA<-summary(manova(lm(f1)))$SS[[3]]/res.shape[3,1]
    eig.FA<-eigen(MSCP.FA); PC1.eigval<-eig.FA$values[1]/sum(eig.FA$values)
    PC1<-shape%*%eig.FA$vec[,1]
    FA.mns<-arrayspecs((rbind(shape[which.min(PC1),],shape[which.max(PC1),])),p,k,byLand=FALSE) 
    symm.component<-arrayspecs(predict(lm(shape~ind)),p,k,byLand=FALSE)
    asymm.component<-array(data=NA,dim=c(p,k,n))  
    mn.shape<-mshape(gpa.res$coords)
    for (i in 1:n){ asymm.component[,,i]<-(gpa.res$coords[,,i]-symm.component[,,i]) + mn.shape}
    if(k==2){  
      #2d plots
      par(mfrow=c(2,2),oma=c(1.5,0,1.5,0))
      plotAllSpecimens(symm.component)
      plotAllSpecimens(asymm.component)
      plotRefToTarget(DA.mns[,,1],DA.mns[,,2],method=plot.method,mag=mag,main="Directional Asymmetry")
      plotRefToTarget(FA.mns[,,1],FA.mns[,,2],method=plot.method,mag=mag,main="Fluctuating Asymmetry")
      mtext("Symmetric Shape Component (left) and Asymmetric Shape Component (right)",outer = TRUE,side=3)
      mtext("Mean directional (left) and fluctuating (right) asymmetry",side = 1, outer = TRUE)
    }
    if (k==3){
      #3d plots
      par(mfrow=c(1,2),oma=c(1.5,0,1.5,0)) 
      plotRefToTarget(DA.mns[,,1],DA.mns[,,2],method=plot.method,mag=mag,main="Directional Asymmetry")
      plotRefToTarget(FA.mns[,,1],FA.mns[,,2],method=plot.method,mag=mag,main="Fluctuating Asymmetry")
      mtext("Mean directional (first plot) and fluctuating (second plot) asymmetry",side = 1, outer = TRUE)
    }
    if(GPAoutput==FALSE){return(list(symm.shape=symm.component,asymm.shape=asymm.component,ANOVA.Shape=res.shape))}
    if(GPAoutput==TRUE){return(list(aligned=gpa.res$coords,csize=gpa.res$Csize,symm.shape=symm.component,asymm.shape=asymm.component,
		ANOVA.Shape=res.shape))}
  }
}

#'  Quantify and compare shape change trajectories
#'
#'  Function estimates attributes of shape change trajectories or motion trajectories for a set of 
#'  Procrustes-aligned specimens and compares them statistically
#'
#'  The function quantifies phenotypic shape change trajectories from a set of specimens, and assesses variation 
#'  in these parameters via permutation. A shape change trajectory is defined by a sequence 
#'  of shapes in tangent space. These trajectories can be quantified various attributes (their size, orientation, 
#'  and shape), and comparisons of these attribute enables the statistical comparison of shape change 
#'  trajectories (see Collyer and Adams 2007; Adams and Collyer 2007; Adams and Collyer 2009). 
#'
#'  Data input is specified by a formula (e.g., Y~X), where 'Y' specifies the response variables (trajectory data), 
#'  and 'X' contains one or more independent variables (discrete or continuous). The response matrix 'Y' must be 
#'  in the form of a two-dimensional data matrix of dimension (n x [p x k]), rather than a 3D array. The function
#'  \code{\link{two.d.array}} can be used to obtain a two-dimensional data matrix from a 3D array of landmark
#'  coordinates. It is assumed that the order of the specimens 'Y' matches the order of specimens in 'X'. 
#' 
#'  There are two primary modes of analysis through this function. If "estimate.traj=TRUE" the function 
#'  estimates shape trajectories using the least-squares means for groups, based on a two-factor model
#'  (e.g., Y~A+B+A:B). Under this implementation, the last factor in 'X' must be the interaction term, and
#'  the preceding two factors must be the effects of interest. Covariates may be included in 'X', and must
#'  precede the factors of interest (e.g., Y~cov+A*B). In this implementation, 'Y' contains a matrix of landmark
#'  coordinates. It is assumed that the landmarks have previously been aligned using Generalized Procrustes 
#'  Analysis (GPA) [e.g., with \code{\link{gpagen}}]. 
#'
#'  If "estimate.traj=FALSE" the trajectories are assembled directly from the set of shapes provided in 'Y'. 
#'  With this implementation, the user must specify the number of shapes that comprise each trajectory. This 
#'  approach is useful when the set of shapes forming each trajectory have been quantified directly 
#'  (e.g., when motion paths are compared: see Adams and Cerney 2007). With this implementation, variation in 
#'  trajectory size, shape, and orientation are evaluated for each term in 'X'.(see Adams and Cerney 2007). 
#'
#'  Once the function has performed the analysis, it generates a plot of the trajectories as visualized in the 
#'  space of principal components (PC1 vs. PC2). The first point in each trajectory is displayed as white, the 
#' last point is black, and any middle points on the trajectories are in gray.  The colors of trajectories follow
#'  the order in which they are found in the dataset, using R's standard color palette: black, red, green3,
#'  blue, cyan, magenta, yellow, and gray. 
#'
#' @param f1 A formula for the linear model (e.g., y~x1+x2)
#' @param data An optional value specifying a data frame containing all data (not required)
#' @param estimate.traj A logical value indicating whether trajectories are estimated from original data; 
#'   described below
#' @param iter Number of iterations for significance testing
#' @param traj.pts An optional value specifying the number of points in each trajectory (if estimate.traj=FALSE)
#' @export
#' @keywords trajectory.analysis
#' @author Dean Adams
#' @return If "estimate.traj=TRUE", the function returns a list with the following components: 
#'   \item{procDist.lm}{Procrustes ANOVA table}
#'   \item{traj.size}{A matrix of pairwise differences in trajectory size}
#'   \item{p.size}{A matrix of pairwise significance levels for trajectory size}
#'   \item{traj.orient}{A matrix of pairwise differences in trajectory orientation}
#'   \item{p.orient}{A matrix of pairwise significance levels for trajectory orientation}
#'   \item{traj.shape}{A matrix of pairwise differences in trajectory shape (if applicable)}
#'   \item{p.shape}{A matrix of pairwise significance levels for trajectory shape}
#' @return If "estimate.traj=FALSE", the function returns a list with the following components: 
#'   \item{MANOVA.location.covariation}{Procrustes ANOVA table}
#'   \item{ANOVA.Size}{Results of permutational-ANOVA assessing variation in trajectory size}
#'   \item{ANOVA.Dir}{Results of permutational-ANOVA assessing variation in trajectory orientation}
#'   \item{ANOVA.Shape}{Results of permutational-ANOVA assessing variation in trajectory shape (if applicable)}
#' @references Collyer, M.L., and D.C. Adams. 2013. Phenotypic trajectory analysis: Comparison of 
#'  shape change patterns in evolution and ecology. Hystrix. 24:75-83.
#' @references Adams, D. C. 2010. Parallel evolution of character displacement driven by competitive 
#'   selection in terrestrial salamanders. BMC Evol. Biol. 10:1-10.
#' @references Adams, D. C., and M. M. Cerney. 2007. Quantifying biomechanical motion using Procrustes 
#'   motion analysis. J. Biomech. 40:437-444.
#' @references Adams, D. C., and M. L. Collyer. 2007. The analysis of character divergence along environmental 
#'   gradients and other covariates. Evolution 61:510-515.
#' @references Adams, D. C., and M. L. Collyer. 2009. A general framework for the analysis of phenotypic 
#'   trajectories in evolutionary studies. Evolution 63:1143-1154.
#' @references Collyer, M. L., and D. C. Adams. 2007. Analysis of two-state multivariate phenotypic change 
#'   in ecological studies. Ecology 88:683-692.
#' @examples
#' #1: Estimate trajectories from LS means in 2-factor model
#' data(plethodon) 
#' Y.gpa<-two.d.array(gpagen(plethodon$land)$coords)    
#'
#' trajectory.analysis(Y.gpa~plethodon$species*plethodon$site,iter=15)
#'
#' #2: Compare motion trajectories
#' data(motionpaths) 
#'
#' #Motion paths represented by 5 time points per motion 
#'
#' trajectory.analysis(motionpaths$trajectories~motionpaths$groups,
#' estimate.traj=FALSE, traj.pts=5,iter=15)
trajectory.analysis<-function(f1,data=NULL,estimate.traj=TRUE,traj.pts=NULL,iter=99){
  form.in<-formula(f1)
  Terms<-terms(form.in)
  y<-eval(form.in[[2]],parent.frame())
  dat<-model.frame(form.in,data)
  ncol.x<-length(attr(Terms,"term.labels"))  
  all.terms<-attr(Terms,"term.labels")
  if (length(dim(y))!=2){
    stop("Response matrix (shape) not a 2D array. Use 'two.d.array' first.")  }
  if(length(grep("-999",y))!=0){
    stop("Response data matrix (shape) contains missing values. Estimate these first(see 'estimate.missing').")  }
  lm.res<-procD.lm(form.in,iter=iter)
  
  if(estimate.traj==TRUE){
    y<-prcomp(y)$x
    if(ncol.x<3){
      stop("X-matrix does not specify enough model factors (see help file).") }
    int.term<-grep(":",attr(Terms,"term.labels")[ncol.x])
    if(int.term!=1){
      stop("Last col of X-matrix does not contain interaction between main effects (see help file).") }          
    nterms<-dim(dat)[2]
    k1<-ncol(y) 
    n1<-length(levels(dat[,nterms-1]))
    p1<-length(levels(dat[,nterms]))
    fac12<-as.factor(paste(dat[,nterms-1],dat[,nterms]))  
    form.full<-as.formula(paste("y ~", paste(all.terms,collapse="+")))
    yhat.full<-predict(lm(form.full))
    lsmeans.obs <- rowsum(yhat.full, fac12)/as.vector(table(fac12))
    form.red<-as.formula(paste("y ~", paste(all.terms[-ncol.x],collapse="+")))
    yhat.red<-predict(lm(form.red))
    res.red<-resid(lm(form.red))
    traj.specs.obs<-arrayspecs(lsmeans.obs,p1,k1) 
    trajsize.obs<-trajsize(traj.specs.obs,n1,p1) 
    trajdir.obs<-trajorient(traj.specs.obs,n1,k1); diag(trajdir.obs)<-0 
    trajshape.obs<-trajshape(traj.specs.obs) 
    PSize<-POrient<-PShape<-array(1,dim=c(n1,n1))
    for(i in 1:iter){
      res.rand<-res.red[sample(nrow(res.red)),]	  
      y.r<-yhat.red+res.rand	
      form.full.r<-as.formula(paste("y.r ~", paste(all.terms,collapse="+")))
      yhat.r<-predict(lm(form.full.r))
      lsmeans.r<-rowsum(yhat.r, fac12)/as.vector(table(fac12))
      traj.specs.r<-arrayspecs(lsmeans.r,p1,k1) 
      trajsize.r<-trajsize(traj.specs.r,n1,p1) 
      trajdir.r<-trajorient(traj.specs.r,n1,k1); diag(trajdir.r)<-0 
      trajshape.r<-trajshape(traj.specs.r) 
      PSize<-ifelse(trajsize.r>=trajsize.obs, PSize+1,PSize) 
      POrient<-ifelse(trajdir.r>=trajdir.obs,POrient+1,POrient) 
      PShape<-ifelse(trajshape.r>=trajshape.obs,PShape+1,PShape) 
    }  
    PSize<-PSize/(iter+1)
    POrient<-POrient/(iter+1)
    PShape<-PShape/(iter+1)
    trajplot(y,traj.specs.obs)
    if(p1>2){
      return(list(ProcDist.lm=lm.res,traj.size=trajsize.obs,p.size=PSize,traj.orient=trajdir.obs,
                  p.orient=POrient,traj.shape=trajshape.obs,p.shape=PShape))
    }
    if(p1<3){
      return(list(ProcDist.lm=lm.res,traj.size=trajsize.obs,p.size=PSize,traj.orient=trajdir.obs,
                  p.orient=POrient))
    }
  }
  
  if(estimate.traj==FALSE){
    if(is.null(traj.pts)==TRUE){
      stop("Number of points in the trajectory not specified.") }
    p1<-traj.pts
    n1<-nrow(y)
    k1<-ncol(y)/p1  
    if (k1>2){
      y.2d<-matrix(t(y),ncol=k1,byrow=TRUE)
      y.2d<-prcomp(y.2d)$x
      y<-two.d.array(arrayspecs(y.2d,p1,k1))
    }        
    traj.specs.obs<-arrayspecs(y,p1,k1,byLand=F) 
    size.obs<-as.dist(trajsize(traj.specs.obs,n1,p1)) 
    dir.obs<-trajorient(traj.specs.obs,n1,k1) 
      diag(dir.obs)<-0; dir.obs<-as.dist(dir.obs)
    shape.obs<-as.dist(trajshape(traj.specs.obs)) 
    size.form<-as.formula(paste("size.obs ~", paste(all.terms,collapse="+")))    
    shape.form<-as.formula(paste("shape.obs ~", paste(all.terms,collapse="+"))) 
    dir.form<-as.formula(paste("dir.obs ~", paste(all.terms,collapse="+")))    
    size.res<-adonis(size.form,permutations=iter)[[1]][1:6] 
    shape.res<-adonis(shape.form,permutations=iter)[[1]][1:6]
    dir.res<-adonis(dir.form,permutations=iter)[[1]][1:6]
    y.plot<-matrix(t(two.d.array(traj.specs.obs)),ncol=k1,byrow=TRUE)
    trajplot(y.plot,traj.specs.obs)
    if(p1>2){
      return(list(MANOVA.location.covariation=lm.res,ANOVA.Size=size.res,ANOVA.Dir=dir.res))
    }
    if(p1<3){
      return(list(MANOVA.location.covariation=lm.res,ANOVA.Size=size.res,ANOVA.Dir=dir.res,ANOVA.Shape=shape.res))
    }
  }
}

#' Plot landmark coordinates for all specimens
#'
#' Function plots landmark coordinates for a set of specimens
#'
#' The function creates a plot of the landmark coordinates for all specimens. This is useful for examining 
#'  patterns of shape variation after GPA. If "mean=TRUE", the mean shape will be calculated and added to the plot.
#'  Additionally, if a matrix of links is provided, the landmarks of the mean shape will be connected by lines.  
#'  The link matrix is an m x 2 matrix, where m is the desired number of links. Each row of the link matrix 
#'  designates the two landmarks to be connected by that link. The function will plot either two- or 
#'  three-dimensional data.
#'
#' @param A An array (p x k x n) containing GPA-aligned coordinates for a set of specimens
#' @param mean A logical value indicating whether the mean shape should be included in the plot
#' @param links An optional matrix defining for links between landmarks
#' @param pointscale An optional value defining the size of the points for all specimens
#' @param meansize An optional value defining the size of the points representing the average specimen
#' @export
#' @keywords plotAllSpecimens
#' @author Dean Adams
#' @examples
#' data(plethodon) 
#' Y.gpa<-gpagen(plethodon$land)    #GPA-alignment
#'
#' plotAllSpecimens(Y.gpa$coords,links=plethodon$links)
plotAllSpecimens<-function(A,mean=TRUE,links=NULL,pointscale=1,meansize=2){
  if (length(dim(A))!=3){
    stop("Data matrix not a 3D array (see 'arrayspecs').")  }
  if(length(grep("-999",A))!=0){
    stop("Data matrix contains missing values. Estimate these first(see 'estimate.missing').")  }
  k<-dim(A)[2]
  if(mean==TRUE){
    mn<-mshape(A)
  }
  if(k==2){
    plot(A[,1,],A[,2,],asp=1, pch=21,bg="gray",xlab="x",ylab="y",cex=pointscale*1) 
    if(mean==TRUE){ 
      points(mn,pch=21,bg="black",cex=meansize) 
      if(is.null(links)==FALSE){
        for (i in 1:nrow(links)){
          segments(mn[links[i,1],1],mn[links[i,1],2],mn[links[i,2],1],mn[links[i,2],2],lwd=2)
        }
      }
    }
  }
  if(k==3){
    A3d<-NULL
    for (i in 1:dim(A)[[3]]){
      A3d<-rbind(A3d,A[,,i])
    }
    plot3d(A3d,type="s",col="gray",xlab="x",ylab="y",zlab="z",size=pointscale*1.5,aspect=FALSE)
    if(mean==TRUE){ 
      points3d(mn,color="black",size=meansize*2)
      if(is.null(links)==FALSE){
        for (i in 1:nrow(links)){
          segments3d(rbind(mn[links[i,1],],mn[links[i,2],]),lwd=3)
        }
      }
    }
  }
}

#' Plot shape differences between a reference and target specimen
#'
#' Function plots shape differences between a reference and target specimen
#'
#' The function generates a plot of the shape differences of a target specimen relative to a reference 
#'  specimen. The option "mag" allows the user to indicates the degree of magnification to be used when 
#'  displaying the shape difference. If "method=TPS" a thin-plate spline deformation grid is generated. 
#'  If "method=vector" a plot showing the vector displacements between corresponding landmarks in the reference 
#'  and target specimen is shown. If "method=points" a plot is displayed with the landmarks in the target (black) 
#'  overlaying those of the reference (gray). Additionally, if a matrix of links is provided, the 
#'  landmarks of the mean shape will be connected by lines.  The link matrix is an M x 2 matrix, where 
#'  M is the desired number of links. Each row of the link matrix designates the two landmarks to be 
#'  connected by that link. The function will plot either two- or three-dimensional data (note: for 3D, 
#'  "method=TPS" will generate thin-plate spline deformations in the x-y and x-z planes). 
#'  This function combines numerous plotting functions found in Claude (2008).
#'
#' @param M1 Matrix of landmark coordinates for the first (reference) specimen
#' @param M2 Matrix of landmark coordinates for the second (target) specimen
#' @param method Method used to visualize shape difference; see below for details
#' @param mag The desired magnification to be used when visualizing the shape difference (e.g., mag=2)
#' @param links An optional matrix defining for links between landmarks
#' @param ... Additional parameters to be passed to 'plot' or 'plot3d'
#' @keywords plotRefToTarget
#' @export
#' @author Dean Adams
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
plotRefToTarget<-function(M1,M2,method=c("TPS","vector","points"),mag=1.0,links=NULL,...){
  method <- match.arg(method)
  if(length(grep("-999",M1))!=0){
    stop("Data contains missing values. Estimate these first(see 'estimate.missing').")  }
  if(length(grep("-999",M2))!=0){
    stop("Data contains missing values. Estimate these first(see 'estimate.missing').")  }
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
      plot(M1,asp=1,pch=21,bg="gray",xlim=limits(M1[,1],1.25),ylim=limits(M1[,2],1.25),cex=1,xlab="x",ylab="y",,...)
      points(M2,asp=1,pch=21,bg="black",cex=1)
      if(is.null(links)==FALSE){
        for (i in 1:nrow(links)){
          segments(M2[links[i,1],1],M2[links[i,1],2],M2[links[i,2],1],M2[links[i,2],2],lwd=1,lty=2)
        }
      }
    }
    if(method=="surface"){
      stop("Surface mesh plotting for 3D landmarks only.")
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
      plot3d(M1,type="s",col="gray",,size=1.25,aspect=FALSE,,...)
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
      plot3d(M1,type="s",col="gray",,size=1.25,aspect=FALSE,,...)
      points3d(M2,color="black",size=5)
      if(is.null(links)==FALSE){
        for (i in 1:nrow(links)){
          segments3d(rbind(M2[links[i,1],],M2[links[i,2],]),lwd=.75)
        }
      }
    }
  }
}

#' Plot specimens in tangent space
#'
#' Function plots a set of Procrustes-aligned specimens in tangent space along their principal axes
#'
#' The function performs a principal components analysis of shape variation and plots two 
#'  dimensions of tangent space for a set of Procrustes-aligned specimens (default is PC1 vs. PC2). The percent variation along each PC-axis 
#'  is returned. Additionally (and optionally), deformation grids can be requested, which display the shape of specimens at the ends 
#'  of the range of variability along PC1. The function returns a table summarizing the percent variation explained by each
#'  pc axis, and the set of principal component scores for all specimens. If groups are provided, specimens from 
#'  each group are plotted using distinct colors based on the order in which the groups are found in the dataset, 
#'  and using R's standard color palette: black, red, green, blue, cyan, magenta, yellow, and gray. 
#'
#' @param A An array (p x k x n) containing landmark coordinates for a set of specimens 
#' @param warpgrids A logical value indicating whether deformation grids for shapes along PC1 should be displayed
#' @param axis1 A value indicating which PC axis should be displayed as the X-axis (default = PC1)
#' @param axis2 A value indicating which PC axis should be displayed as the X-axis (default = PC2)
#' @param label A logical value indicating whether labels for each specimen should be displayed
#' @param groups An optional factor vector specifying group identity for each specimen
#' @export
#' @keywords plotTangentSpace
#' @author Dean Adams
#' @examples
#' data(plethodon) 
#' Y.gpa<-gpagen(plethodon$land)    #GPA-alignment
#' ref<-mshape(Y.gpa$coords)
#' 
#' plotTangentSpace(Y.gpa$coords, groups = paste(plethodon$species, plethodon$site)) 
plotTangentSpace<-function (A, axis1 = 1, axis2 = 2, warpgrids = TRUE, label = FALSE,groups=NULL){
  if (length(dim(A)) != 3) {
    stop("Data matrix not a 3D array (see 'arrayspecs').") }
  if (length(grep("-999", A)) != 0) {
    stop("Data matrix contains missing values. Estimate these first(see 'estimate.missing').") }
  if(!is.null(groups)){groups<-factor(groups)}
  k <- dim(A)[2]
  p <- dim(A)[1]
  n <- dim(A)[3]
  ref <- mshape(A)
  x <- two.d.array(A)
  pc.res <- prcomp(x)
  pcdata <- pc.res$x
  if (warpgrids == F) {
    plot(pcdata[, axis1], pcdata[, axis2], asp = 1, pch = 21,bg = "black", cex = 2, xlab = paste("PC ",axis1),
	 ylab = paste("PC ",axis2))
      if(!is.null(groups)){points(pcdata[, axis1], pcdata[, axis2],pch=21,bg=groups,cex=2)}
    segments(min(pcdata[, axis1]), 0, max(pcdata[, axis1]), 0, lty = 2, lwd = 1)
    segments(0, min(pcdata[, axis2]), 0, max(pcdata[, axis2]), lty = 2, lwd = 1)
    if (label == T) {text(pcdata[, axis1], pcdata[, axis2], seq(1, n), adj = c(-0.7, -0.7)) }
  }
  if (warpgrids == T) {
    if (k == 2) {
      layout(t(matrix(c(2, 1, 1, 1, 1, 1, 1, 1, 3), 3,3)))
    }
    plot(pcdata[, axis1], pcdata[, axis2], asp = 1, pch = 21,bg = "black", cex = 2, xlab = paste("PC ",axis1),
	 ylab = paste("PC ",axis2))
      if(!is.null(groups)){points(pcdata[, axis1], pcdata[, axis2],pch=21,bg=groups,cex=2)}
      segments(min(pcdata[, axis1]), 0, max(pcdata[, axis1]), 0, lty = 2, lwd = 1)
      segments(0, min(pcdata[, axis2]), 0, max(pcdata[, axis2]), lty = 2, lwd = 1)
      if (label == T) {text(pcdata[, axis1], pcdata[, axis2], seq(1, n), adj = c(-0.7, -0.7))}
      pcaxis.min <- min(pcdata[, axis1]); pcaxis.max <- max(pcdata[, axis1])
      pc.min<-pc.max<-rep(0, dim(pcdata)[2])
      pc.min[axis1]<-pcaxis.min; pc.max[axis1]<-pcaxis.max
      shape.min <- arrayspecs(as.matrix(pc.min %*% (t(pc.res$rotation))), 
                            p, k, byLand = FALSE)[, , 1] + ref
      shape.max <- arrayspecs(as.matrix(pc.max %*% (t(pc.res$rotation))), 
                            p, k, byLand = FALSE)[, , 1] + ref
    if (k == 2) {
      arrows(min(pcdata[, axis1]), (0.7 * max(pcdata[,axis2])), min(pcdata[, axis1]), 0, length = 0.1,lwd = 2)
      arrows(max(pcdata[, axis1]), (0.7 * min(pcdata[,axis2])), max(pcdata[, axis1]), 0, length = 0.1,lwd = 2)
      tps(ref, shape.min, 20)
      tps(ref, shape.max, 20)
    }
    if (k == 3) {
      open3d()
      plot3d(shape.min, type = "s", col = "gray", main = paste("PC ", axis1," negative"),size = 1.25, aspect = FALSE)
      open3d()
      plot3d(shape.max, type = "s", col = "gray", main = paste("PC ", axis1," positive"), size = 1.25, aspect = FALSE)
    }
    layout(1)
  }
  return(list(pc.summary = summary(pc.res), pc.scores = pcdata))
}


#' Plot phylogenetic tree and specimens in tangent space
#'
#' Function plots a phylogenetic tree and a set of Procrustes-aligned specimens in tangent space
#'
#' The function creates a plot of the first two dimensions of tangent space for a set of Procrustes-aligned 
#'   specimens. The phylogenetic tree for these specimens is superimposed in this plot revealing how shape 
#'   evolves (e.g., Rohlf 2002; Klingenberg and Gidaszewski 2010). The plot also displays the ancestral 
#'   states for each node of the phylogenetic tree, whose values can optionally be returned. 
#'
#' @param phy A phylogenetic tree of type 'phylo'
#' @param A An array (p x k x n) containing landmark coordinates for a set of specimens 
#' @param labels A logical value indicating whether taxa labels should be included
#' @param ancStates A logical value indicating whether ancestral state values should be returned
#' @export
#' @keywords plotGMPhyloMorphoSpace
#' @author Dean Adams
#' @references Klingenberg, C. P., and N. A. Gidaszewski. 2010. Testing and quantifying phylogenetic 
#'   signals and homoplasy in morphometric data. Syst. Biol. 59:245-261.
#' @references Rohlf, F. J. 2002. Geometric morphometrics and phylogeny. Pp. 175'193 in N. Macleod, and 
#'   P. Forey, eds. Morphology, shape, and phylogeny. Taylor & Francis, London.
#' @examples
#' data(plethspecies) 
#' Y.gpa<-gpagen(plethspecies$land)    #GPA-alignment    
#'
#' plotGMPhyloMorphoSpace(plethspecies$phy,Y.gpa$coords)
plotGMPhyloMorphoSpace<-function(phy,A,labels=TRUE,ancStates=T){
  if (length(dim(A))!=3){
    stop("Data matrix not a 3D array (see 'arrayspecs').")  }
  if(length(grep("-999",A))!=0){
    stop("Data matrix contains missing values. Estimate these first(see 'estimate.missing').")  }
  if (class(phy) != "phylo") 
    stop("tree must be of class 'phylo.'")
  if (!is.binary.tree(phy)) 
    stop("tree is not fully bifurcating.")
  if(is.null(dimnames(A)[[3]])){
    stop("Data matrix does not include taxa names as dimnames for 3rd dimension.")  }
  N<-length(phy$tip.label)
  x<-two.d.array(A)
  if(N!=dim(x)[1]){
    stop("Number of taxa in data matrix and tree are not not equal.")  }
  if(length(match(rownames(x), phy$tip.label))!=N) 
    stop("Data matrix missing some taxa present on the tree.")
  if(length(match(phy$tip.label,rownames(x)))!=N) 
    stop("Tree missing some taxa in the data matrix.")
  x<-x[phy$tip.label, ]	
  names<-row.names(x)
  anc.states<-NULL
  for (i in 1:ncol(x)){
    options(warn=-1)  
    tmp<-as.vector(ace(x[,i], compute.brlen(phy,1), type="continuous", method="ML")$ace)
    anc.states<-cbind(anc.states,tmp)   }
  colnames(anc.states)<-NULL
  all.data<-rbind(x,anc.states)  
  pcdata<-prcomp(all.data)$x  
  limits = function(x,s){ 
   r = range(x)
   rc=scale(r,scale=F)
   l=mean(r)+s*rc}
  if(labels==TRUE){
    plot(pcdata,type="n",xlim=limits(pcdata[,1],1.5),ylim=limits(pcdata[,2],1.5),asp=1) }
  if(labels==FALSE) {
    plot(pcdata,type="n",asp=1) }
  for (i in 1:nrow(phy$edge)){
    lines(pcdata[(phy$edge[i,]),1],pcdata[(phy$edge[i,]),2],type="l",pch=21,col="black",lwd=3)
  }
  points(pcdata[1:N,],pch=21,bg="black",cex=2)
  points(pcdata[(N+1):nrow(pcdata),],pch=21,bg="white",cex=1.25)
  if(labels==TRUE){
    textxy(pcdata[1:N,1],pcdata[1:N,2],rownames(pcdata),cx=.75)  }
  if(ancStates==TRUE){ return(anc.states)  }
}

#' Plot allometric patterns in landmark data
#'
#' Function plots allometry curves for a set of specimens
#'
#' The function generates a plot that describes the multivariate relationship between size and shape 
#'   derived from landmark data (ie. allometry). It is assumed that the landmarks have previously been 
#'   aligned using Generalized Procrustes Analysis (GPA) [e.g., with \code{\link{gpagen}}]. The abscissa 
#'   of the plot is log(centroid size) while the ordinate represents shape. Three complementary approaches 
#'   can be implemented to visualize allometry: 
#'  \enumerate{
#'   \item {If "method=CAC" (the default) the function calculates the 
#'   common allometric component of the shape data, which is an estimate of the average allometric trend 
#'   within groups (Mitteroecker et al. 2004). }
#'   \item {If "method=RegScore" the function calculates shape scores 
#'   from the regression of shape on size, and plots these versus size (Drake and Klingenberg 2008). 
#'   For a single group, these shape scores are mathematically identical to the CAC (Adams et al. 2012).}
#'   \item {If "method=PredLine" the function calculates predicted values from a regression of shape on size, and 
#'   plots the first principal component of the predicted values versus size as a stylized graphic of the 
#'   allometric trend (Adams and Nistri 2010). }
#'   }
#'   For all methods, both centroid size and allometry scores are returned. Optionally, deformation grids can be 
#'   requested, which display the shape of the smallest and largest specimens relative to the average specimen (using 
#'   'warpgrid=T' or 'warpgrid=F'). 
#'   Finally, if groups are provided, the above approaches are implemented while 
#'   accounting for within-group patterns of covariation (see references for explanation). In this case,
#'   specimens from each group are plotted using distinct colors based on the order in which the groups are
#'   found in the dataset, and using R's standard color palette: black, red, green, blue, cyan, magenta,
#'   yellow, and gray. 
#'
#' @param A An array (p x k x n) containing landmark coordinates for a set of specimens 
#' @param sz A vector of centroid size measures for all specimens 
#' @param groups An optional vector containing group labels for each specimen if available 
#' @param method Method for estimating allometric shape components; see below for details
#' @param warpgrids A logical value indicating whether deformation grids for small and large shapes 
#'  should be displayed
#' @param label A logical value indicating whether labels for each specimen should be displayed
#' @keywords plotAllometry
#' @export
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
#' plotAllometry(Y.gpa$coords,Y.gpa$Csize,method="CAC")
#'
#' #Using Regression Scores for plot
#' plotAllometry(Y.gpa$coords,Y.gpa$Csize,method="RegScore")
#'
#' #Using predicted allometry curve for plot
#' plotAllometry(Y.gpa$coords,Y.gpa$Csize,method="PredLine")
plotAllometry<-function(A,sz,groups=NULL,method=c("CAC","RegScore","PredLine"),warpgrids=TRUE,label=FALSE){
  method <- match.arg(method)
  if (length(dim(A))!=3){
    stop("Data matrix 1 not a 3D array (see 'arrayspecs').")  }
  if(length(grep("-999",A))!=0){
    stop("Data matrix contains missing values. Estimate these first(see 'estimate.missing').")  }
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
  Reg.proj<-y%*%B[2,]%*%sqrt(solve(t(B[2,])%*%B[2,])) 
  pred.val<-prcomp(yhat)$x[,1] 
  if(warpgrids==F){
    if(method=="CAC"){
      plot(csz,CAC,xlab="log(CSize)", ylab="CAC",pch=21,bg="black",cex=1.25)
      if(!is.null(groups)){points(csz,CAC,pch=21,bg=groups,cex=1.25)}
      if(label ==T){text(csz,CAC,seq(1,n),adj=c(-.7,-.7))}
    }
    if(method=="RegScore"){
      plot(csz,Reg.proj,xlab="log(CSize)", ylab="Shape (Regression Score)",pch=21,bg="black",cex=1.25)
      if(!is.null(groups)){points(csz,Reg.proj,pch=21,bg=groups,cex=1.25)}
      if(label ==T){text(csz,Reg.proj,seq(1,n),adj=c(-.7,-.7))}
    }
    if(method=="PredLine"){
      plot(csz,pred.val,xlab="log(CSize)", ylab="Shape (Predicted)",pch=21,bg="black",cex=1.25)
      if(!is.null(groups)){points(csz,pred.val,pch=21,bg=groups,cex=1.25)}
      if(label ==T){text(csz,pred.val,seq(1,n),adj=c(-.7,-.7))}
    }
  }
  if(warpgrids==T){
    k<-dim(A)[2]
    p<-dim(A)[1]
    ref<-mshape(A)
    if(k==2){  
     layout(t(matrix(c(2,1,1,1,1,1,1,1,3),3,3)))   
    }
    if(method=="CAC"){
      plot(csz,CAC,xlab="log(CSize)", ylab="CAC",pch=21,bg="black",cex=1.25)
      if(!is.null(groups)){points(csz,CAC,pch=21,bg=groups,cex=1.25)}
      if(label ==T){text(csz,CAC,seq(1,n),adj=c(-.7,-.7))}
      mypar<-par("usr")
      if(k==2){  
        text((mypar[1]+.2*(mypar[2]-mypar[1])),.5*mypar[4],"Shape at minimum size")
        text((mypar[2]-.2*(mypar[2]-mypar[1])),.5*mypar[3],"Shape at maximum size")
      }
    }
    if(method=="RegScore"){
      plot(csz,Reg.proj,xlab="log(CSize)", ylab="Shape (Regression Score)",pch=21,bg="black",cex=1.25)
      if(!is.null(groups)){points(csz,Reg.proj,pch=21,bg=groups,cex=1.25)}
      if(label ==T){text(csz,Reg.proj,seq(1,n),adj=c(-.7,-.7))}
      mypar<-par("usr")
      if(k==2){  
        text((mypar[1]+.2*(mypar[2]-mypar[1])),.5*mypar[4],"Shape at minimum size")
        text((mypar[2]-.2*(mypar[2]-mypar[1])),.5*mypar[3],"Shape at maximum size")      }
    }
    if(method=="PredLine"){
      plot(csz,pred.val,xlab="log(CSize)", ylab="Shape (Predicted)",pch=21,bg="black",cex=1.25)
      if(!is.null(groups)){points(csz,pred.val,pch=21,bg=groups,cex=1.25)}
      if(label ==T){text(csz,pred.val,seq(1,n),adj=c(-.7,-.7))}
      mypar<-par("usr")
      if(k==2){  
        text((mypar[1]+.2*(mypar[2]-mypar[1])),.5*mypar[4],"Shape at minimum size")
        text((mypar[2]-.2*(mypar[2]-mypar[1])),.5*mypar[3],"Shape at maximum size")
      }
    }
    if(k==2){
      tps(ref,A[,,which.min(csz)],20)
      tps(ref,A[,,which.max(csz)],20)
      text(1,9,"Shape of small specimen")
    }
    if(k==3){
     open3d()
      plot3d(A[,,which.min(csz)],type="s",col="gray",main="Smallest Specimen",size=1.25,aspect=FALSE)
      open3d()
      plot3d(A[,,which.max(csz)],type="s",col="gray",main="Largest Specimen",size=1.25,aspect=FALSE)
    }
  layout(1)  #reset plot layout
  }
  if(method=="CAC"){a.scr=CAC}
  if(method=="RegScore"){a.scr=Reg.proj}
  if(method=="PredLine"){a.scr=pred.val}
  return(list(allom.score=a.scr,logCsize=csz))
}

#### TPS and GPA routines (DCA and J Claude code) 

tps<-function(matr, matt, n){		#DCA: altered from J. Claude: 2D only	
  xm<-min(matt[,1])
  ym<-min(matt[,2])
  xM<-max(matt[,1])
  yM<-max(matt[,2])
  rX<-xM-xm; rY<-yM-ym
  a<-seq(xm-1/5*rX, xM+1/5*rX, length=n)
  b<-seq(ym-1/5*rX, yM+1/5*rX,by=(xM-xm)*7/(5*(n-1)))
  m<-round(0.5+(n-1)*(2/5*rX+ yM-ym)/(2/5*rX+ xM-xm))
  M<-as.matrix(expand.grid(a,b))
  ngrid<-tps2d(M,matr,matt)
  plot(ngrid, cex=0.2,asp=1,axes=FALSE,xlab="",ylab="")
  for (i in 1:m){lines(ngrid[(1:n)+(i-1)*n,])}
  for (i in 1:n){lines(ngrid[(1:m)*n-i+1,])}
  points(matt,pch=21,bg="black",cex=1.5)
}

tps2d<-function(M, matr, matt)
{p<-dim(matr)[1]; q<-dim(M)[1]; n1<-p+3
 P<-matrix(NA, p, p)
 for (i in 1:p)
 {for (j in 1:p){
   r2<-sum((matr[i,]-matr[j,])^2)
   P[i,j]<- r2*log(r2)}}
 P[which(is.na(P))]<-0
 Q<-cbind(1, matr)
 L<-rbind(cbind(P,Q), cbind(t(Q),matrix(0,3,3)))
 m2<-rbind(matt, matrix(0, 3, 2))
 coefx<-solve(L)%*%m2[,1]
 coefy<-solve(L)%*%m2[,2]
 fx<-function(matr, M, coef)
 {Xn<-numeric(q)
  for (i in 1:q)
  {Z<-apply((matr-matrix(M[i,],p,2,byrow=TRUE))^2,1,sum)
   Xn[i]<-coef[p+1]+coef[p+2]*M[i,1]+coef[p+3]*M[i,2]+sum(coef[1:p]*(Z*log(Z)))}
  Xn}
 matg<-matrix(NA, q, 2)
 matg[,1]<-fx(matr, M, coefx)
 matg[,2]<-fx(matr, M, coefy)
 matg}


tps2d3d<-function(M, matr, matt){		#DCA: altered from J. Claude 2008  
  p<-dim(matr)[1]; k<-dim(matr)[2];q<-dim(M)[1]
  Pdist<-as.matrix(dist(matr))
  ifelse(k==2,P<-Pdist^2*log(Pdist^2),P<-Pdist) 
  P[which(is.na(P))]<-0
  Q<-cbind(1, matr)
  L<-rbind(cbind(P,Q), cbind(t(Q),matrix(0,k+1,k+1)))
  m2<-rbind(matt, matrix(0, k+1, k))   
  coefx<-solve(L)%*%m2[,1]
  coefy<-solve(L)%*%m2[,2]
  if(k==3){coefz<-solve(L)%*%m2[,3]}
  fx<-function(matr, M, coef){
    Xn<-numeric(q)
    for (i in 1:q){
      Z<-apply((matr-matrix(M[i,],p,k,byrow=TRUE))^2,1,sum)  
      ifelse(k==2,Z1<-Z*log(Z),Z1<-sqrt(Z)); Z1[which(is.na(Z1))]<-0
      ifelse(k==2,Xn[i]<-coef[p+1]+coef[p+2]*M[i,1]+coef[p+3]*M[i,2]+sum(coef[1:p]*Z1),
             Xn[i]<-coef[p+1]+coef[p+2]*M[i,1]+coef[p+3]*M[i,2]+coef[p+4]*M[i,3]+sum(coef[1:p]*Z1))
    }    
    Xn}
  matg<-matrix(NA, q, k)
  matg[,1]<-fx(matr, M, coefx)
  matg[,2]<-fx(matr, M, coefy)
  if(k==3){matg[,3]<-fx(matr, M, coefz)}  
  matg
}

scan.to.ref<-function(scandata,specland,refland){  	#DCA
  ref.scan<-tps2d3d(scandata,specland,refland)
  ref.scan}

refscan.to.spec<-function(refscan,refland,specland){ 	#DCA
  unwarp.scan<-tps2d3d(refscan,refland,specland)
  unwarp.scan}

trans<-function(A){scale(A,scale=FALSE)} 		#J. Claude 2008

csize<-function(A)				#J. Claude 2008
{p<-dim(A)[1]
 size<-sqrt(sum(apply(A,2,var))*(p-1))
 list("centroid_size"=size,"scaled"=A/size)}

#' Estimate mean shape for a set of aligned specimens
#'
#' Estimate the mean shape for a set of aligned specimens
#'
#' The function estimates the average landmark coordinates for a set of aligned specimens. It is assumed 
#' that the landmarks have previously been aligned using Generalized Procrustes Analysis (GPA) 
#'  [e.g., with \code{\link{gpagen}}]. This function is described in Claude (2008).
#'
#' @param A An array (p x k x n) containing GPA-aligned coordinates for a set of specimens
#' @keywords mshape
#' @export
#' @author Julien Claude 
#' @references Claude, J. 2008. Morphometrics with R. Springer, New York.
#' @examples
#' data(plethodon) 
#' Y.gpa<-gpagen(plethodon$land)    #GPA-alignment   
#'
#' mshape(Y.gpa$coords)   #mean (consensus) configuration
mshape<-function(A){apply(A,c(1,2),mean)}	

pPsup<-function(M1,M2){				#J. Claude 2008# MODIFIED BY DCA: 2/2014
  k<-ncol(M1)
  Z1<-trans(csize(M1)[[2]])
  Z2<-trans(csize(M2)[[2]])
  sv<-svd(t(Z2)%*%Z1)
  U<-sv$v; V<-sv$u; Delt<-sv$d
  sig<-sign(det(t(Z2)%*%Z1))
  Delt[k]<-sig*abs(Delt[k]); V[,k]<-sig*V[,k]
  Gam<-U%*%t(V)
  # beta<-sum(Delt)   #commented out: retain size-scaling (DCA)
#  list(Mp1=beta*Z1%*%Gam,Mp2=Z2,rotation=Gam,scale=beta,
#       df=sqrt(1-beta^2))}
  list(Mp1=Z1%*%Gam,Mp2=Z2,rotation=Gam)}    #simplify output and remove size re-scaling (DCA) 

pgpa<-function(A)				#J. Claude 2008	
{p<-dim(A)[1]; k<-dim(A)[2]; n<-dim(A)[3]  
 temp2<-temp1<-array(NA,dim=c(p,k,n)); Siz<-numeric(n)#; Qm2<-numeric(n)
 for (i in 1:n)
 {Acs<-csize(A[,,i])
  Siz[i]<-Acs[[1]]
  temp1[,,i]<-trans(Acs[[2]])}
 Qm1<-dist(t(matrix(temp1,k*p,n)))
 Q<-sum(Qm1); iter<-0
 while (abs(Q)> 0.0001)
 {for(i in 1:n){
   M<-mshape(temp1[,,-i])
   temp2[,,i]<-pPsup(temp1[,,i],M)[[1]]}
  Qm2<-dist(t(matrix(temp2,k*p,n)))
  Q<-sum(Qm1)-sum(Qm2)
  Qm1<-Qm2
  iter=iter+1
  temp1<-temp2}
 list("rotated"=temp2,"it.number"=iter,"Q"=Q,"intereucl.dist"=Qm2,"mshape"=
   csize(mshape(temp2))[[2]],"cent.size"=Siz)
}

orp<-function(A){			#DCA: altered from J. Claude 2008		 
  n<-dim(A)[3]; k<-dim(A)[2]; p<-dim(A)[1]  
  Y1<-as.vector(csize(mshape(A))[[2]])
  oo<-as.matrix(rep(1,n))%*%Y1
  mat<-matrix(NA,n,k*p)
  for (i in 1:n){mat[i,]<-as.vector(A[,,i])}
  Xp<-(mat%*%(diag(1,p*k)- (Y1%*%t(Y1))))+oo
  array(t(Xp),dim=c(p,k,n))
}

# Trajectory Size: Pathlength Distance
pathdist<-function(M) {as.matrix(dist(M))} 
trajsize<-function(M,n,p){
  traj.pathdist<-array(0,dim=c(n,1))   		
  for (i in 1:n){
    temp<-pathdist(M[,,i])
    for (j in 1:(p-1)){
      traj.pathdist[i]<-traj.pathdist[i]+temp[j,j+1]
    }
  }
  traj.size.dist<-as.matrix(dist(traj.pathdist))		
}

# Trajectory Orientation
trajorient<-function(M,n,k){
  traj.orient<-array(NA,dim=c(n,k))   
  check.1<-array(NA,dim=c(n))
  for (i in 1:n){
    temp<-svd(var(M[,,i]))$v[1:k,1]
    traj.orient[i,]<-temp
    check.1[i]<-M[1,,i]%*%traj.orient[i,]  
    check.1[i]<-check.1[i]/abs(check.1[i])
    if(check.1[i]==-1) traj.orient[i,]<--1*traj.orient[i,]
  }
  options(warn=-1)				
  traj.ang.diff<-(180/pi)*acos(traj.orient%*%t(traj.orient))
}

# Trajectory Shape
trajshape<-function(M){
  x<-pgpa(M)
  traj.shape.dist<-as.matrix(x$intereucl.dist) 
}

# general plotting function for phenotypic trajectories
trajplot<-function(Data,M){
  n<-dim(M)[3]; p<-dim(M)[1]
  plot(Data[,1:2],type="n",xlab="Summary Axis I", ylab="Summary Axis II",main="Two Dimensional View  of Phenotypic Trajectories",asp=1)
  points(Data[,1:2],pch=21,bg="gray",cex=.75)
  for (i in 1:n){  	 	
    for (j in 1:(p-1)){		
      points(M[(j:(j+1)),1,i],M[(j:(j+1)),2,i],type="l",pch=21,col=i)  #was black    
    }
  }
  for (i in 1:n){		 	
    for (j in 2:(p-1)){		
      points(M[j,1,i],M[j,2,i],pch=21,bg="gray",col="black",cex=1.5)
    }
  }
  for (i in 1:n){
    points(M[1,1,i],M[1,2,i],pch=21,bg="white",col="black",cex=1.5)
  }
  for (i in 1:n){
    points(M[p,1,i],M[p,2,i],pch=21,bg="black",col="black",cex=1.5)
  }
}

# Write .nts file for output of digitize2d(), buildtemplate() digit.fixed() and digitsurface()
# A is an nx2 or nx3 matrix of the output coordinates. To be used internally only.

writeland.nts <- function(A, spec.name, comment=NULL){
  ntsfile=paste(spec.name,".nts",sep="")
  file.create(file=ntsfile, overwrite = TRUE)
  if(is.null(comment)){
    cat(paste('"',spec.name,sep=""),file= ntsfile,sep="\n",append=TRUE)
  }
  else if(!is.null(comment)){
    cat(paste('"',spec.name,sep=""),file= ntsfile,sep="\n")
    cat(paste('"',comment,sep=""),file= ntsfile,sep="\n",append=TRUE)
  }
  dims <- dim(A)
  if (dims[2] == 2){
    cat(paste(1,dims[1],2,0,"dim=2"),file= ntsfile,sep="\n",append=TRUE)
  }
  else if (dims[2] == 3){
    cat(paste(1,dims[1],3,0, "dim=3"),file= ntsfile,sep="\n",append=TRUE)
  }
  write.table(A ,file= ntsfile,col.names = FALSE, row.names = FALSE,sep="  ",append=TRUE)
}


#' Build 3D surface template 
#'
#' An interactive function to build template of three-dimensional (3D) surface sliding semilandmarks.
#' Input for the function is either a matrix of vertex coordinates defining a 3D surface object
#' or a mesh3d object as obtained from \code{\link{read.ply}}.
#'
#' Function constructs a template of fixed landmarks and n "surface sliders", semilandmarks that slide over a surface.
#'  The user digitizes the fixed points (see digitizing below), then the function finds n surface semilandmarks following 
#'  algorithm outlined in Gunz et al. (2005:90-92) and Mitteroecker and Gunz (2009:242). Surface semilandmarks are 
#'  roughly equidistant set of predetermined number of points, chosen over the mesh automatically using a nearest-neighbor approach.
#'  The set of fixed and surface slider landmarks are exported as a "template", which is used to extract a set of similarly numbered 
#'  landmarks on every specimen using function \code{\link{digitsurface}}. Some of the "fixed" landmarks can be later designated as
#'  "curve sliders" using function \code{\link{define.sliders.3d}} if required - see details in \code{\link{digit.fixed}}.
#'  
#' To ensure a strong match between the scan and the template, it is recommended that a reasonable number of fixed points be used.
#'  These fixed points can be designated as "curve sliders" later using function \code{\link{define.sliders.3d}}, 
#'  see the function \code{\link{digit.fixed}} for details.
#'  NOTE: Function centers the mesh before digitizing (required for further operations).
#' 
#' \subsection{Digitizing}{
#' Digitizing using buildtemplate is interactive between landmark selection using a mouse (see below for instructions), 
#'  and the R console. Once a point is selected, the user is asked if the system should keep or discard the 
#'  selection #'(y/n). If "y", the user is asked to continue to select the next landmark. If "n" the removes the last chosen
#'  landmark, and the user is asked to select it again. This can be repeated until the user is comfortable with the landmark
#'  chosen. 
#'
#' To digitize with a standard 3-button (PC):
#' \enumerate{
#'  \item the RIGHT mouse button (primary) to select points to be digitized (click-drag a box around a vertex to select as landmark)
#'  \item the LEFT mouse button (secondary) is used to rotate mesh, 
#'  \item the mouse SCROLLER (third/middle) is used to zoom in and out.
#' }
#' NOTE: Digitizing functions on MACINTOSH computers using a standard 3-button mice works as specified. 
#' Macs using platform specific single button mice: 
#' \enumerate{
#'  \item press button to rotate 3D mesh,
#'  \item press button while pressing CONTROL key to select points to be digitized (click-drag a box around a vertex to select as landmark), 
#'  \item press button while pressing OPTION key to adjust mesh perspective.
#' }
#' Mac mouse settings need adjustment to acquire zooming functions of the "middle/third" mouse button/scroller.
#' Some instructions can be found at \url{http://blog.boastr.net/}. In addition, freeware providing "Middle Click" 
#' functionality is available at \url{http://magicprefs.com/} for "magicmice" now standard 
#' on many Macintosh machines.
#' 
#' NOTE: there is no pan (translate) functionality in rgl library for all platforms at this time.
#' The template can be edited using function \code{\link{editTemplate}}. 
#' }
#'
#' @param spec Name of surface file, as either an object of class shape3d/mesh3d, or matrix of three-dimensional vertex coordinates.
#' @param fixed numeric: The number of fixed template landmarks
#' @param surface.sliders numeric: The number of template surface sliders desired 
#' @param ptsize numeric: Size to plot the mesh points (vertices), e.g. 0.1 for dense meshes, 3 for sparse meshes
#' @export
#' @keywords buildtemplate
#' @author Erik Otarola-Castillo & Emma Sherratt
#' @return Function writes to the working directory three files: an NTS file with the name of the specimen and .nts suffix containing 
#' the landmark coordinates, "template.txt" containing the same coordinates for use with the function \code{\link{digitsurface}}, 
#' and "surfslide.csv", a file containing the address of the landmarks defined as "surface sliders" for use with \code{\link{gpagen}}.
#' Function also returns to console an n x 3 matrix containing the x,y,z coordinates of the digitized landmarks. 
#' @references Gunz P, Mitteroecker P, & Bookstein FJ (2005) Semilandmarks in Three Dimensions. Modern Morphometrics in Physical 
#' Anthropology, ed Slice DE (Springer-Verlag, New York), pp 73-98.
#' @references Mitteroecker P & Gunz P (2009) Advances in Geometric Morphometrics. Evolutionary Biology 36(2):235-247.

buildtemplate<-function(spec, fixed, surface.sliders, ptsize = 1)    {
  spec.name<-deparse(substitute(spec))
  mesh <- NULL
  if (inherits(spec, "shape3d") == TRUE || inherits(spec, "mesh3d") == TRUE){
    specimen <- scale(as.matrix(t(spec$vb)[,-4]), scale = FALSE)
    spec$vb <- rbind(t(specimen), 1)
    mesh <- spec 
    if (is.null(mesh$material)) { mesh$material <- "gray" } 
  } else if (inherits(spec, "matrix") == FALSE) {
    stop ("File is not a shape3d/mesh3d object or xyz matrix")
  } else if (inherits(spec, "matrix") == TRUE && dim(spec)[2]==3) {
    specimen <- scale(spec, scale = FALSE)
  } else { stop ("File is not matrix in form: vertices by xyz")} 
  selected<-digit.fixed(spec,fixed,index=TRUE,ptsize)
  fix<-selected$fix
  selected<-selected$selected
  surfs<-specimen[-fix,]
  template<-rbind(selected,kmeans(x=surfs,centers=surface.sliders,iter.max=100)$centers)
  points3d(template[-(1:fixed),1],template[-(1:fixed),2],template[-(1:fixed),3],size=10,col="blue")
  writeland.nts(template, spec.name, comment="Landmark coordinates for digitized template")
  write.table(template,file="template.txt",row.names=F,col.names=TRUE)
  write.csv(seq(from=(fixed+1),to=(surface.sliders+fixed)),file="surfslide.csv",row.names=FALSE)  
  rownames(template)<-seq(from=1,to=nrow(template))
  return(template)
}

#' Choose points to "slide" along three-dimensional curves.
#'
#' An interactive function to define which digitized landmarks of an '*.nts' file
#' will "slide" along three-dimensional (3D) curves.
#' 
#' Function takes a matrix of digitized landmark coordinates, such as made by \code{\link{digit.fixed}}, or \code{\link{buildtemplate}}
#'  and helps user choose which landmarks will be treated as "curve sliders" in Generalized Procrustes analysis
#'  \code{\link{gpagen}}. This type of semilandmark "slides" along curves lacking known landmarks 
#'  (see Bookstein 1991:376-382, 1997 for algorithm details). 
#'  Each sliding semilandmark ("sliders") will slide between two designated points, along a line 
#'  tangent to the specified curvature.
#' 
#' \subsection{Selection}{ 
#' Choosing which landmarks will be sliders involves landmark selection using a mouse in the rgl plot window. 
#' With a standard 3-button (PC) buildtemplate uses:
#' \enumerate{
#'  \item the RIGHT mouse button (primary) to choose points to be defined as sliders (click-drag a box around a vertex to select),
#'  \item the LEFT mouse button (secondary) is used to rotate mesh, 
#'  \item the mouse SCROLLER (third/middle) is used to zoom in and out.
#' }
#' NOTE: Digitizing functions on MACINTOSH computers using a standard 3-button mice works as specified. 
#' Macs using platform specific single button mice: 
#' \enumerate{
#'  \item press button to rotate 3D mesh,
#'  \item press button while pressing CONTROL key to choose points to be defined as sliders,
#'  \item press button while pressing OPTION key to zoom in and out.}
#'
#' Mac mouse settings need adjustment to acquire zooming functions of the "middle/third" mouse button/scroller.
#' Some instructions can be found at \url{http://blog.boastr.net/}. In addition, freeware providing "Middle Click" 
#' functionality is available at \url{http://magicprefs.com/} for "magicmice" now standard 
#' on many Macintosh machines.
#' 
#' To define the sliders, for each sliding landmark along the curve in the format 'before-slider-after':
#' \enumerate{
#'  \item Click-drag to choose the first landmark between which semi-landmark will "slide",
#'  \item Click-drag to choose sliding landmark,
#'  \item Click-drag to choose the last landmark between which semi-landmark will "slide",
#' Screen will show lines connecting the three landmarks, and will highlight the sliding semilandmark in red. 
#' }
#' This procedure is overlapping, so for example a curve defined by a sequence of semilandmarks,
#' the user must select the 2nd point of the first three to be the 1st for the next
#' e.g. 1 2 3 then 2 3 4, etc.
#' }
#' @param spec Name of specimen, as an object matrix containing 3D landmark coordinates
#' @param nsliders Number of landmarks to be semilandmarks that slide along curves
#' @param surfsliders Logical If spec contains landmarks that are "surface sliders",
#' made by \code{\link{buildtemplate}}, "surfslide.csv" should be in working directory 
#' @export
#' @seealso  \code{\link{digit.fixed}}, \code{\link{digitsurface}}, \code{\link{gpagen}}
#' @keywords define.sliders.3d
#' @author Erik Otarola-Castillo & Emma Sherratt
#' @return Function returns a 'curves x 3' matrix containing the landmark address of the curve sliders, indicating the 
#' points between which the selected point will "slide". Written to the working directory is this matrix as "curveslide.csv".
#' Matrix (or "curveslide.csv") is designed for use by \code{\link{gpagen}} during GPA.
#' @references  Bookstein, F. J. 1991  Morphometric Tools for Landmark Data: Geometry and Biology. 
#' Cambridge University Press, New York.
#' @references Bookstein, F. J. 1997 Landmark Methods for Forms without Landmarks: Morphometrics of 
#' Group Differences in Outline Shape. Medical Image Analysis 1(3):225-243.
define.sliders.3d<-function(spec, nsliders,surfsliders=FALSE)    {
  spec.name <- deparse(substitute(spec))
  rownames(spec) <- c(1:nrow(spec)) 
  checkmat <- is.matrix(spec)
  if (checkmat==FALSE) { stop("Input must be a p-x-k matrix of landmark coordinates")}
  checkdim <- dim(spec)[2]
  if (checkdim==2) {stop("Input must be a p-x-k matrix of three-dimensional landmark coordinates") }
  if (surfsliders == TRUE){
    surf <- as.matrix(read.csv("surfslide.csv", header=T))
    spec <- spec[-surf,]
  }
  n <- dim(spec)[1]
  index <- as.numeric(rownames(spec))
  clear3d();plot3d(spec,size=5,col="black",xlab="x",ylab="y",zlab="z",aspect=FALSE)
  text3d(spec, texts=index,cex=1,adj=c(2,1))  
  curveslide<-NULL 
    for (i in 1:nsliders)      {
      selected<-NULL
      for (j in 1:3)    		{
        f<-keep<-NULL
        f<-select3d(button="right")
        keep<-f(spec[,1],spec[,2],spec[,3])
        selected.temp<-cbind(index[which(keep==TRUE)][1],spec[,1][which(keep==TRUE)][1],spec[,2][which(keep==TRUE)][1],spec[,3][which(keep==TRUE)][1])
        selected<-rbind(selected,selected.temp)
        if (j==2) {
          points3d(selected[2,2],selected[2,3],selected[2,4],size=10,color="red",add=TRUE) 
        } else {
          points3d(selected[,2],selected[,3],selected[,4],size=7,color="blue",add=TRUE)
        }  
        if (j>=2) {
          lines3d(selected[c(j-1,j),2],selected[c(j-1,j),3],selected[c(j-1,j),4],size=10,color="red",add=TRUE)
        } 
      }
      cat("semi-landmark",selected[2,1],"slides between landmarks",selected[1,1],"and",selected[3,1],"\n")
      curveslide<-rbind(curveslide,selected[,1])
    }
  write.table(curveslide,file="curveslide.csv",row.names=FALSE,col.names=c("before","slide","after"),sep=",")
  return(curveslide)
}



#' Defunct functions in geomorph
#'
#' The following functions are no longer supported in geomorph
#'
#' This function is now defunct. Below shows the original function and the replacement function.
#'
#' curves2d:  use define.sliders.2d
#'
#' digit.curves: use define.sliders.3d
#' @export
digit.curves<-function()    {
  .Defunct(define.sliders.3d)
}

#' Digitize 3D landmarks
#'
#' An interactive function to digitize three-dimensional (3D) landmarks.
#' Input for the function is either a matrix of vertex coordinates defining a 3D surface object
#' or a mesh3d object as obtained from \code{\link{read.ply}}. 
#'
#' Function for digitizing "n" three-dimensional landmarks. The landmarks are "fixed" (traditional landmarks). They can be later
#'  designated as "curve sliders" (semilandmarks, that will "slide" along curves lacking known landmarks (see Bookstein 1991:376-382, 1997 for algorithm details) 
#'  if required. A sliding semi-landmark ("sliders") will slide between two designated points, along a line tangent to the specified curvature,
#'  and must be defined as "sliders" using function \code{\link{define.sliders.3d}} or with similar format matrix made outside R. 
#' 
#' For 3D "surface sliders" (surface semilandmarks that slide over a surface) the function \code{\link{digitsurface}} should be used instead.
#'  NOTE: Function centers the mesh before digitizing (required for further operations).
#' 
#' \subsection{Digitizing}{ 
#' Digitizing is interactive between landmark selection using a mouse (see below for instructions), 
#' and the R console. Once a point is selected, the user is asked if the system should keep or discard the 
#' selection #'(y/n). If "y", the user is asked to continue to select the next landmark. If "n" the removes the last chosen
#' landmark, and the user is asked to select it again. This can be repeated until the user is comfortable with the landmark
#' chosen. 
#' 
#' To digitize with a standard 3-button (PC):
#' \enumerate{
#'  \item the RIGHT mouse button (primary) to select points to be digitized (click-drag a box around a vertex to select as landmark),
#'  \item the LEFT mouse button (secondary) is used to rotate mesh, 
#'  \item the mouse SCROLLER (third/middle) is used to zoom in and out.
#' }
#' NOTE: Digitizing functions on MACINTOSH computers using a standard 3-button mice works as specified. 
#' Macs using platform specific single button mice: 
#' \enumerate{
#'  \item press button to rotate 3D mesh,
#'  \item press button while pressing CONTROL key to select points to be digitized (click-drag a box around a vertex to select as landmark),
#'  \item press button while pressing OPTION key to to zoom in and out.}
#'
#' Mac mouse settings need adjustment to acquire zooming functions of the "middle/third" mouse button/scroller.
#' Some instructions can be found at \url{http://blog.boastr.net/}. In addition, freeware providing "Middle Click" 
#' functionality is available at \url{http://magicprefs.com/} for "magicmice" now standard 
#' on many Macintosh machines.
#' 
#' NOTE: there is no pan (translate) functionality in rgl library for all platforms at this time.
#' }
#'
#' @param spec Name of surface file, as either an object of class shape3d/mesh3d, or matrix of three-dimensional vertex coordinates.
#' @param fixed numeric: The number landmarks (fixed, and curve sliders if desired)
#' @param index logical: Whether selected landmark addresses should be returned
#' @param ptsize numeric: Size to plot the mesh points (vertices), e.g. 0.1 for dense meshes, 3 for sparse meshes 
#' @return Function writes to the working directory an NTS file with the name of the specimen and .nts suffix containing the landmark coordinates.
#' If index=FALSE function returns to the console an n x 3 matrix containing the x,y,z coordinates of the digitized landmarks.  
#' If index=TRUE, function returns a list:
#'  \item{selected}{a matrix containing the x,y,z coordinates of the digitized landmarks} 
#'  \item{fix}{a matrix of addresses for landmarks that are "fixed" (for internal use)} 
#' @export
#' @keywords digit.fixed
#' @author Erik Otarola-Castillo & Emma Sherratt
digit.fixed <- function(spec, fixed, index=FALSE, ptsize = 1)    {
  spec.name<-deparse(substitute(spec))
  mesh <- NULL
  if (inherits(spec, "shape3d") == TRUE || inherits(spec, "mesh3d") == TRUE){
    specimen <- scale(as.matrix(t(spec$vb)[,-4]), scale = FALSE)
    spec$vb <- rbind(t(specimen), 1)
    mesh <- spec 
    if (is.null(mesh$material)) { mesh$material <- "gray" }
  } else if (inherits(spec, "matrix") == FALSE) {
    stop ("File is not a shape3d/mesh3d object or xyz matrix")
  } else if (inherits(spec, "matrix") == TRUE && dim(spec)[2]==3) {
    specimen <- scale(spec, scale = FALSE)
  } else { stop ("File is not matrix in form: vertices by xyz")} 
  clear3d();plot3d(specimen[,1],specimen[,2],specimen[,3],size=ptsize,aspect=FALSE)
  if (!is.null(mesh)) { shade3d(mesh, add=TRUE) }
  selected<-matrix(NA,nrow=fixed,ncol=3);fix<-NULL    
  for (i in 1:fixed)      {
    f<-keep<-ans<-NULL
    f<-select3d(button="right")
    keep<-f(specimen)
    selected[i,]<-as.numeric(specimen[which(keep==TRUE)[1],])
    fix<-c(fix,which(keep==TRUE)[1])   
    points3d(selected[i,1],selected[i,2],selected[i,3],size=10,color="red",add=TRUE)
    cat(paste("Keep Landmark ",i,"(y/n)?"),"\n")
    ans<-readLines(n=1)
    if(ans=="y" & length(fix)!=fixed) {
      cat("Select Landmark ",i+1,"\n")
      rgl.bringtotop(stay = FALSE)
    } 
    if(ans=="n" ) {
      cat(paste("Select Landmark ",i," Again"),"\n")
    }
    while(ans=="n") {
      selected[i,]<-NA
      fix<-fix[-i] 
      rgl.bringtotop(stay = FALSE)
      clear3d();plot3d(specimen[,1],specimen[,2],specimen[,3],size=ptsize,aspect=FALSE)
      if (!is.null(mesh)) { shade3d(mesh, add=TRUE) }
      if(sum(1-is.na(selected))>0){
        points3d(selected[,1],selected[,2],selected[,3],size=10,color="red",add=TRUE)
      }      
      f<-keep<-NULL
      f<-select3d(button="right")
      keep<-f(specimen)
      selected[i,]<-as.numeric(specimen[which(keep==TRUE)[1],])
      fix<-c(fix,which(keep==TRUE)[1])   
      points3d(selected[i,1],selected[i,2],selected[i,3],size=10,color="red",add=TRUE)
      cat(paste("Keep Landmark ",i,"(y/n)?"),"\n")
      ans<-readLines(n=1)
      if(ans=="y" & length(fix)!=fixed) {
        cat("Select Landmark ",i+1,"\n")
      } 
      if(ans=="n") {
        cat(paste("Select Landmark ",i," Again"),"\n")
      }
    } 
  } 
  if(index==FALSE){
    writeland.nts(selected, spec.name, comment=NULL)
    rownames(selected)<-seq(from=1,to=nrow(selected))
    colnames(selected)<-c("xpts","ypts","zpts")
    return(selected) 
  }
  if(index==TRUE){
    return(list(selected=selected,fix=fix))
  }
}


#' Digitize 3D fixed landmarks and surface semilandmarks.
#' 
#' An interactive function to digitize three-dimensional (3D) landmarks on a surface lacking known landmarks.
#' Input for the function is either a matrix of vertex coordinates defining a 3D surface object
#' or a mesh3d object as obtained from \code{\link{read.ply}}. 
#' 
#' Function for digitizing fixed 3D landmarks and placing "surface sliders", semilandmarks that slide over a surface.
#'  Following selection of fixed points (see digitizing below), function finds surface semilandmarks following 
#'  algorithm outlined in Gunz et al. (2005:90-92) and Mitteroecker and Gunz (2009:242). digitsurface finds the same number of surface 
#'  semilandmarks as the template (created by \code{\link{buildtemplate}}) by downsampling scanned mesh, registering template with 
#'  current specimen via GPA. A nearest neighbor algorithm is used to match template surface landmarks to current specimen's. 
#'  To use function digitsurface, the template must be constructed first, and 'template.txt' be in the working directory. 
#' 
#' Some of the "fixed" landmarks digitized with digitsurface can be later designated as "curve sliders" using function \code{\link{define.sliders.3d}} 
#'  if required (see details in \code{\link{digit.fixed}}).
#'  NOTE: Function centers the mesh before digitizing (required for further operations).
#' 
#' \subsection{Digitizing}{
#' Digitizing is interactive between landmark selection using a mouse (see below for instructions), 
#' and the R console. Once a point is selected, the user is asked if the system should keep or discard the 
#' selection #'(y/n). If "y", the user is asked to continue to select the next landmark. If "n" the removes the last chosen
#' landmark, and the user is asked to select it again. This can be repeated until the user is comfortable with the landmark
#' chosen. 
#' 
#' To digitize with a standard 3-button (PC):
#' \enumerate{
#'  \item the RIGHT mouse button (primary) to select points to be digitized (click-drag a box around a vertex to select as landmark),
#'  \item the LEFT mouse button (secondary) is used to rotate mesh, 
#'  \item the mouse SCROLLER (third/middle) is used to zoom in and out.
#' }
#' NOTE: Digitizing functions on MACINTOSH computers using a standard 3-button mice works as specified. 
#' Macs using platform specific single button mice: 
#' \enumerate{
#'  \item press button to rotate 3D mesh,
#'  \item press button while pressing CONTROL key to select points to be digitized (click-drag a box around a vertex to select as landmark), 
#'  \item press button while pressing OPTION key to adjust mesh perspective.
#' }
#' Mac mouse settings need adjustment to acquire zooming functions of the "middle/third" mouse button/scroller.
#' Some instructions can be found at \url{http://blog.boastr.net/}. In addition, freeware providing "Middle Click" 
#' functionality is available at \url{http://magicprefs.com/} for "magicmice" now standard 
#' on many Macintosh machines.
#' 
#' NOTE: there is no pan (translate) functionality in rgl library for all platforms at this time.
#' }
#'
#' @param spec Name of surface file, as either an object of class shape3d/mesh3d, or matrix of three-dimensional vertex coordinates.
#' @seealso \code{\link{buildtemplate}}
#' @param fixed numeric: The number of fixed template landmarks 
#' @param ptsize numeric: Size to plot the mesh points (vertices), e.g. 0.1 for dense meshes, 3 for sparse meshes
#' @return Function writes to the working directory an NTS file with the name of the specimen and .nts suffix containing the landmark coordinates.   
#' Function also returns a list:
#'  \item{FIX.LANDMARKS}{an n x 3 matrix containing the x,y,z coordinates of the digitized landmarks} 
#'  \item{SURFSLIDERSxyz}{a matrix of addresses for landmarks that are "surface sliders"} 
#' @references Gunz P, Mitteroecker P, & Bookstein FJ (2005) Semilandmarks in Three Dimensions. Modern Morphometrics in Physical Anthropology, ed Slice DE (Springer-Verlag, New York), pp 73-98.
#' @references Mitteroecker P & Gunz P (2009) Advances in Geometric Morphometrics. Evolutionary Biology 36(2):235-247. 
#' @export
#' @keywords digitsurface
#' @author Erik Otarola-Castillo & Emma Sherratt
digitsurface<-function(spec, fixed, ptsize = 1)    {
  spec.name<-deparse(substitute(spec))
  mesh <- NULL
  if (inherits(spec, "shape3d") == TRUE || inherits(spec, "mesh3d") == TRUE){
    specimen <- scale(as.matrix(t(spec$vb)[,-4]), scale = FALSE)
    spec$vb <- rbind(t(specimen), 1)
    mesh <- spec 
    if (is.null(mesh$material)) { mesh$material <- "gray" }
  } else if (inherits(spec, "matrix") == FALSE) {
    stop ("File is not a shape3d/mesh3d object or xyz matrix")
  } else if (inherits(spec, "matrix") == TRUE && dim(spec)[2]==3) {
    specimen <- scale(spec, scale = FALSE)
  } else { stop ("File is not matrix in form: vertices by xyz")} 
  selected<-digit.fixed(spec, fixed,index=T,ptsize)
  template<-as.matrix(read.table("template.txt",header=TRUE))
  specimen<-trans(as.matrix(specimen))
  template<-trans(template)*(csize(specimen[selected$fix,])[[1]]/csize(template[(1:fixed),])[[1]])  
  template<-template%*%(pPsup(template[(1:fixed),],specimen[selected$fix,]))[[3]] 
  template.tps<-tps2d3d(template[-(1:fixed),],template[(1:fixed),],specimen[selected$fix,])             
  spec.surfs<-specimen[-selected$fix,]
  nei<-numeric(dim(template.tps)[1])
  sliders<-matrix(NA,nrow=dim(template.tps)[1],ncol=3)
  for (i in 1:dim(template.tps)[1])     {
    nei[i]<-which.min(sqrt((template.tps[i,1]-spec.surfs[,1])^2+(template.tps[i,2]-spec.surfs[,2])^2+(template.tps[i,3]-spec.surfs[,3])^2))[1] #3D NN
    sliders[i,]<-spec.surfs[nei[i],]
    spec.surfs<-spec.surfs[-nei[i],]  
  }
  clear3d();plot3d(specimen[,1],specimen[,2],specimen[,3],size=ptsize,aspect=FALSE)
  if (!is.null(mesh)) { 
    mesh$vb <- rbind(t(specimen), 1)
    shade3d(mesh, add=TRUE) 
  }
  points3d(specimen[selected$fix,],col="red",size=10)
  points3d(template.tps,col="blue",size=10)
  points3d(sliders[,1:3],size=10,col="green")
  selected.out <- rbind(specimen[selected$fix,],sliders)
  writeland.nts(selected.out, spec.name, comment=NULL)
  return(list(FIX.LANDMARKS=selected$selected,SURFSLIDERSxyz=sliders))
}

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
#' NOTE: Digitizing functions on MACINTOSH computers using a standard 3-button mice works as specified. 
#' Macs using platform specific single button mice: 
#' \enumerate{
#'  \item press button to rotate 3D mesh,
#'  \item press button while pressing CONTROL key to choose points to be deleted (click-drag a box around landmark),
#'  \item press button while pressing OPTION key to zoom in and out.
#'  }
#' Mac mouse settings need adjustment to acquire zooming functions of the "middle/third" mouse button/scroller.
#' Some instructions can be found at \url{http://blog.boastr.net/}. In addition, freeware providing "Middle Click" 
#' functionality is available at \url{http://magicprefs.com/} for "magicmice" now standard 
#' on many Macintosh machines.
#' }
#' @param template Matrix of template 3D coordinates.
#' @param fixed Number of "fixed" landmark points (non surface sliding points)
#' @param n Number of points to be removed 
#' @return Function returns a matrix containing the x,y,z coordinates of the new template landmarks.
#' Function also writes to working directory 'template.txt' containing the x,y,z coordinates of the template 
#' @export
#' @keywords editTemplate
#' @author Erik Otarola-Castillo & Emma Sherratt
editTemplate<-function(template, fixed, n){
  if (is.null(dim(template))) stop ("File is not a matrix of 3D coordinates.")
  if (dim(template)[2]!=3) stop ("File is not a matrix of 3D coordinates.")
  spec.name<-deparse(substitute(template))
  clear3d();plot3d(template[-(1:fixed[1]),],size=7,col="blue",aspect=FALSE)
  points3d(template[(1:fixed),],size=10,color="red",add=TRUE)  
  cat("Remove Template Points","\n")
  selected2<-NULL
  for (i in 1:n)      	{
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


#' Digitize 2D landmarks.
#'
#' An interactive function to digitize two-dimensional(2D) landmarks from .jpg files.
#'
#' Function for digitizing 2D landmarks on specimen images (.jpg). "nlandmarks" 
#'  are digitized by the user (see below). The landmarks are "fixed" (traditional landmarks). 
#'  They can be later designated as "curve sliders" (semilandmarks, that will "slide" along curves 
#'  lacking known landmarks (see Bookstein 1991:376-382, 1997 for algorithm details) 
#'  if required. These landmarks are defined as "sliders" using function \code{\link{define.sliders.2d}} 
#'  or with similar format matrix made outside R. Digitizing sequence should proceed 
#'  by selecting all "fixed" landmark points first, followed by selection of "curve sliders". 
#' 
#' \subsection{Digitizing}{ 
#' Digitizing landmarks involves landmark selection using a mouse in the plot window, 
#' using the LEFT mouse button (or regular button for Mac users):
#' \enumerate{
#'  \item Digitize the scale bar by selecting the two end points (single click for start and end),
#'  NOTE: no dot is printed on the image on the first click; only when both ends have been digitized.
#'  \item Digitize each landmark with single click and the landmark is shown in red,
#'  When selection of n landmarks is completed, an ".nts" file is created in working directory using the specimen name, adding ".nts" as a suffix.
#' }
#' }
#' @param spec Name of jpeg file to be digitized. File names can be 
#' written in manually, including paths, or obtained using directory/file manipulation functions 
#' e.g., \code{\link{list.files}}
#' @seealso \code{\link{list.files}}, 
#' @param nlandmarks Number of landmarks to be digitized.
#' @param scale Length of scale placed in image.
#' @return Function writes to the working directory an NTS file with the name of the specimen and .nts suffix containing the landmark coordinates.   
#' Function also returns to console an n x 2 matrix containing the x,y coordinates of the digitized landmarks. 
#' @export 
#' @keywords digitize2d
#' @author Erik Otarola-Castillo & Emma Sherratt
digitize2d<-function(spec, nlandmarks, scale){
  spec.name<-unlist(strsplit(basename(file), "\\."))[1]
  specimen<-readJPEG(file, native = T)  
  plot(seq(0,dim(specimen)[2],length.out=10),seq(0,dim(specimen)[1],length.out=10), type='n',xlab="x",ylab="y",asp=1,tck=0,xaxt="n",yaxt="n")
  rasterImage(specimen, 1, 1,dim(specimen)[2],  dim(specimen)[1])
  cat("Set scale =",scale,"\n")
  scalebar<-picscale(scale)
  selected <- matrix(NA, nrow = nlandmarks, ncol = 2)
  fix <- NULL
  cat("Select landmarks 1:",nlandmarks,"\n",sep="")
    for (i in 1:nlandmarks) {
      fix<- locator(n=1, type ="p",col="black",cex=4,pch=21,bg="red")
      selected[i,1] <- fix$x
      selected[i,2] <- fix$y  
    }
  output<-selected/scalebar 
  writeland.nts(output, spec.name, comment=NULL)  
  output<-matrix(output,dimnames=list(paste("LM",seq(1,nlandmarks)),c("x","y")),ncol=2)
  return(LANDMARKS=output)
}

# picscale is called by digitize2d
picscale<- function(scale){
  digscale<-NULL
  digscale<-locator(2,type="o",lwd=2,col="red",lty="11")
  sqrt(sum(diff(digscale$x)^2+diff(digscale$y)^2))*scale
}

#' Defunct functions in geomorph
#'
#' The following functions are no longer supported in geomorph
#'
#' This function is now defunct. Below shows the original function and the replacement function.
#'
#' curves2d:  use define.sliders.2d
#'
#' digit.curves: use define.sliders.3d
#' @export
curves2d<-function(){
  .Defunct(define.sliders.2d)
}


#' Select points to "slide" along two-dimensional curves.
#'
#' An interactive function to define which landmarks will "slide" along two-dimensional (2D) curves.
#'
#' Function takes a matrix of digitized landmark coordinates, such as made by \code{\link{digitize2d}},
#'  and helps user choose which landmarks will be treated as "curve sliders" in Generalized Procrustes analysis
#'  \code{\link{gpagen}}. This type of semilandmark "slides" along curves lacking known landmarks 
#'  (see Bookstein 1991:376-382, 1997 for algorithm details). 
#'  Each sliding semilandmark ("sliders") will slide between two designated points, along a line 
#'  tangent to the specified curvature.
#' 
#' \subsection{Selection}{ 
#' Choosing which landmarks will be sliders involves landmark selection using a mouse in the plot window. 
#' To define the sliders, for each sliding landmark along the curve in the format 'before-slider-after',
#' using the LEFT mouse button (or regular button for Mac users), click on the hollow circle to choose the landmark
#' in the following order:
#' \enumerate{
#'  \item Click to choose the first landmark between which semi-landmark will "slide",
#'  \item Click to choose sliding landmark,
#'  \item Click to choose the last landmark between which semi-landmark will "slide",
#'  Selected landmarks will be filled in and lines are drawn connecting the three landmarks, 
#'  and will highlight the sliding semilandmark in red and the flanking landmarks in blue. 
#' }
#' This procedure is overlapping, so for example, a curve defined by a sequence of semilandmarks,
#' the user must select the 2nd point of the first three to be the 1st for the next
#' e.g. 1 2 3 then 2 3 4, then 3 4 5 etc.
#' }
#' @param spec Name of specimen, as an object matrix containing 2D landmark coordinates
#' @param nsliders Number of landmarks to be semilandmarks that slide along curves
#' @return Function returns a 'nsliders-x-3' matrix containing the landmark address of the curve sliders, indicating the 
#' landmarks between which the slider landmarks will "slide". The matrix is also written to working directory
#' as "curveslide.csv". Matrix (or "curveslide.csv") is designed for use by \code{\link{gpagen}} during GPA.
#' @export
#' @keywords define.sliders.2d
#' @seealso  \code{\link{digitize2d}}, \code{\link{gpagen}}
#' @author Dean Adams,
#' @author Erik Otarola-Castillo
#' @author Emma Sherratt
#' @references  Bookstein, F. J. 1991  Morphometric Tools for Landmark Data: Geometry and Biology. 
#' Cambridge University Press, New York.
#' @references Bookstein, F. J. 1997 Landmark Methods for Forms without Landmarks: Morphometrics of 
#' Group Differences in Outline Shape. Medical Image Analysis 1(3):225-243.
define.sliders.2d<-function(spec, nsliders){
  spec.name <- deparse(substitute(spec))
  checkmat <- is.matrix(spec)
  if (checkmat==FALSE) { stop("Input must be a p-x-k matrix of landmark coordinates")}
  checkdim <- dim(spec)[2]
  if (checkdim==3) {stop("Input must be a p-x-k matrix of two-dimensional landmark coordinates") }
  plot(spec[,1],spec[,2],cex=1,pch=21,bg="white",xlim=range(spec[,1]),ylim=range(spec[,2]),asp=1)
  text(spec[,1],spec[,2],label=paste(1:dim(spec)[1]),adj=.5,pos=4)
  selected<-matrix(NA,ncol=3,nrow=nsliders)
  select<-NULL
  for(i in 1:nsliders){
    for(j in 1:3){
      select<-identify(spec,n=1,plot=FALSE,cex=5,pch=25)
      selected[i,j]<-select
      if(j==2){
        points(spec[select,][1],spec[select,][2],cex=1.5,pch=19,col="red")
        arrows(spec[selected[i,j],][1],spec[selected[i,j],][2],spec[selected[i,j-1],][1],spec[selected[i,j-1],][2]
               ,col="red",lwd=2,length=.15)
      } else {
        points(spec[select,][1],spec[select,][2],cex=1.1,pch=19,col="blue")
      }
      if(j==3){
        arrows(spec[selected[i,j],][1],spec[selected[i,j],][2],spec[selected[i,j-1],][1],spec[selected[i,j-1],][2]
               ,col="red",lwd=2,length=.15,code=1)
      } else { NA
      } 
    }
  }
  write.table(selected,file="curveslide.csv",row.names=FALSE,col.names=c("before","slide","after"),sep=",")
  return(curveslide=selected)
}

#' Plot 3D specimen, fixed landmarks and surface semilandmarks.
#'
#' A function to plot three dimensional (3D) specimen along with its landmarks.
#' 
#' Function to plot 3D specimens along with their digitized "fixed" landmarks and semilandmarks
#' "surface sliders" and "curve sliders".
#' If specimen is a 3D surface (ply file), mesh is plotted.
#'
#' @param spec Name of surface file, as either an object of class shape3d/mesh3d, or matrix of three-dimensional vertex coordinates.
#' @param digitspec Name of data matrix containing three-dimensional fixed and/or surface sliding coordinates.
#' @param fixed numeric: The number of fixed template landmarks  
#' @param ptsize numeric: Size to plot the mesh points (vertices), e.g. 0.1 for dense meshes, 3 for sparse meshes                                       
#' @export
#' @keywords plotspec
#' @examples
#' data(scallopPLY)
#' ply <- scallopPLY$ply
#' digitdat <- scallopPLY$coords
#' plotspec(spec=ply,digitspec=digitdat,fixed=16)
#' @author Erik Otarola-Castillo & Emma Sherratt
plotspec <- function (spec, digitspec, fixed, ptsize = 1) 
{
  mesh <- NULL
  if (inherits(spec, "shape3d") == TRUE || inherits(spec, "mesh3d") == TRUE){
    specimen <- scale(as.matrix(t(spec$vb)[,-4]), scale = FALSE)
    spec$vb <- rbind(t(specimen), 1)
    mesh <- spec 
    if (is.null(mesh$material)) { mesh$material <- "gray" }
  } else if (inherits(spec, "matrix") == FALSE) {
    stop ("File is not a shape3d/mesh3d object or xyz matrix")
  } else if (inherits(spec, "matrix") == TRUE && dim(spec)[2]==3) {
    specimen <- scale(spec, scale = FALSE)
  } else { stop ("File is not matrix in form: vertices by xyz")} 
  if (is.null(dim(digitspec)) || dim(digitspec)[2] != 3) stop("Digitized file is not xyz matrix in form: p x k")
  clear3d()
  plot3d(specimen[, 1], specimen[, 2], specimen[, 3], size = ptsize, 
         aspect = FALSE)
  if (!is.null(mesh)) { shade3d(mesh, add=TRUE) }
  points3d(digitspec[1:fixed, ], aspect = FALSE, size = 10, 
           col = "red")
  points3d(digitspec[(fixed + 1):nrow(digitspec), ], aspect = F, 
           size = 10, col = "green")
}

#' Read landmark data from vrml files
#'
#' Read vrml files (Virtual Reality Modeling Language) to obtain landmark coordinates and triangulations
#'
#' This function reads three-dimensional surface data from simple vrml files
#' (Virtual Reality Modeling Language). The landmarks of this surface may then be 
#'  used to digitize three-dimensional points, and semilandmarks on curves and surfaces. 
#' .vrml files are stored either as centralized data, where 3D scanned object information 
#' (i.e., coordinates, triangles, color, and surface) are stored within single data blocks for
#' the complete object. This version of read.vrml will import .wrl files written in "utf8" and "ascii" format. Mesh triangle facets will be imported if present as coordinate connections. Argument plotspec allows users to plot file. This is helpful to help inspect object for potential
#' errors. Argument plottri plots triangle facets if present. write.nts provides users with the option of writing a .nts coordinate file.
#'
#' NOTE: At present, only simple vrml files containing vertices and faces is supported.
#'
#' @param file A .vrml file with coordinates in a "centralized" format, or in "stitched" format. File names 
#' can be written in manually, including paths, or obtained using directory/file manipulation functions e.g., 
#' \code{\link{list.files}}
#' @seealso \code{\link{list.files}} 
#' @param plotspec Logical should object be plotted. Defaults to TRUE.
#' @param plottri Logical should triangles be plotted. Defaults to TRUE.
#' @param write.nts Logical should .nts file be created. Defaults to FALSE.
#' @export
#' @keywords read.vrml
#' @author Erik Otarola-Castillo
#' @return Function returns a list with the following components:
#'   \item{coords}{The x,y,z coordinates of the .vrml surface}
#'   \item{triangles}{Triangle facet connections between coordinate points (if available)}
read.vrml<-function(file,plotspec=TRUE,plottri=TRUE,write.nts=FALSE) {
  spec<-basename(file)
  spec.name<-unlist(strsplit(spec, "\\."))[1]
  b<-scan(spec,what="character",skip=0,sep="",quiet=T)
  if (sum(b[1:20]=="utf8")==0 & sum(b[1:20]=="ascii")==0) {stop('Not a recognized VRML format. Use only "ascii" or "utf8"',call.=FALSE)}
  if (sum(b[1:20]=="utf8")==1){
    cos<-which(b=="Coordinate")
  }
  if (sum(b[1:20]=="ascii")==1){
    cos<-which(b=="Coordinate3")
  }
  tri<-which(b=="coordIndex")
  if (length(cos)>1) {
    zz<-triang<-NULL
    for (i in 1:length(cos)){ 
      zerocoord<-NULL
      co<-cos[i]+4 
      if (length(tri)>0){
        co.tri<-tri[i]+2
        newrange.tri<-b[co.tri:length(b)]
        en.tri<-which(b[co.tri:length(b)]=="]")[1]-1
        tria<-as.matrix(newrange.tri[1:en.tri])
        triang.temp<-matrix(as.numeric(unlist(strsplit(tria,split=","))),ncol=4,byrow=T)[,-4]
      }
      newrange<-b[co:length(b)]
      en<-(which(b[co:length(b)]=="]")[1]-1)
      z<-as.matrix(newrange[1:en])
      zz.temp<-matrix(as.numeric(unlist(strsplit(z,split=","))),ncol=3,byrow=T)
      zerocoord<-which(as.character(zz.temp[,1])==as.character(0) & as.character(zz.temp[,2])==as.character(0) & as.character(zz.temp[,3])==as.character(0))
      if (length(zerocoord)>0){
        zz.temp<-zz.temp[-zerocoord,]
      }
      if (length(tri)>0){
        if(i>1){triang.temp<-triang.temp + dim(zz)[1]}
        triang<-rbind(triang,triang.temp)
      }
      zz<-rbind(zz,zz.temp)      
    }     
    colnames(zz)<-c("x","y","z")
    if(plotspec==TRUE){
      plot3d(zz[,1],zz[,2],zz[,3],xlab=paste(spec.name,"_x",sep=""),ylab=paste(spec.name,"_y",sep=""),zlab=paste(spec.name,"_z",sep=""),size=1,col="black",aspect=FALSE)
    }     
    if(plottri==TRUE & plotspec==TRUE & length(tri)>0){
      rgl.triangles(zz[t(triang),])
    }    
    if(plottri==TRUE & plotspec==FALSE){
      plot3d(zz[,1],zz[,2],zz[,3],type="n",xlab=paste(spec.name,"_x",sep=""),ylab=paste(spec.name,"_y",sep=""),zlab=paste(spec.name,"_z",sep=""),size=1,col="black",aspect=FALSE)
      if(plottri==TRUE & plotspec==FALSE & length(tri)>0){
        rgl.triangles(zz[t(triang),])
      }
    }
  }    
  if (length(cos)==1) {
    co<-cos[1]+4
    coor<-length(co)
    co.tri<-tri[1]+2
    if (length(tri)>0){
      co.tri<-tri[1]+2
      newrange.tri<-b[co.tri:length(b)]
      en.tri<-which(b[co.tri:length(b)]=="]")[1]-1
      tria<-as.matrix(newrange.tri[1:en.tri])
      triang<-matrix(as.numeric(unlist(strsplit(tria,split=","))),ncol=4,byrow=T)[,-4]
    }
    newrange<-b[co:length(b)]
    en<-(which(b[co:length(b)]=="]")[1]-1)
    z<-as.matrix(newrange[1:en])
    zz<-matrix(as.numeric(unlist(strsplit(z,split=","))),ncol=3,byrow=T)
    zerocoord<-which(as.character(zz[,1])==as.character(0) & as.character(zz[,2])==as.character(0) & as.character(zz[,3])==as.character(0))
    if (length(zerocoord)>0){
      zz<-zz[-zerocoord,]
    }
    colnames(zz)<-c("x","y","z")
    if(plotspec==TRUE){
      plot3d(zz[,1],zz[,2],zz[,3],xlab=paste(spec.name,"_x",sep=""),ylab=paste(spec.name,"_y",sep=""),zlab=paste(spec.name,"_z",sep=""),size=1,col="black",aspect=FALSE)
    } else {NA}    
    if(plottri==TRUE & plotspec==TRUE & length(tri)>0){
      rgl.triangles(zz[t(triang),])
    }    
    if(plottri==TRUE & plotspec==FALSE){
      plot3d(zz[,1],zz[,2],zz[,3],type="n",xlab=paste(spec.name,"_x",sep=""),ylab=paste(spec.name,"_y",sep=""),zlab=paste(spec.name,"_z",sep=""),size=1,col="black",aspect=FALSE)
      if(plottri==TRUE & plotspec==FALSE & length(tri)>0){
        rgl.triangles(zz[t(triang),])
      }
    }
  }  
  if(write.nts==TRUE){
    output<-zz
    cat(paste('"',spec.name,sep=""),file=paste(spec.name,"coords.nts",sep=""),sep="\n")
    cat(paste(1,dim(output)[1],3,0, "dim=3"),file=paste(spec.name,"coords.nts",sep=""),sep="\n",append=TRUE)
    write.table(output,file=paste(spec.name,"coords.nts",sep=""),col.names = FALSE, row.names = FALSE,sep="  ",append=TRUE)
  }
  if (length(tri)>0){
    return(list(coords=zz,triangles=triang))
  } else {
    return(zz)
  }
}
