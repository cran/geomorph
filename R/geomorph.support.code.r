#' @name geomorph-package
#' @docType package
#' @aliases geomorph
#' @title Geometric morphometric analyses for 2D/3D data
#' @author Dean C. Adams, Michael Collyer, Erik Otarola-Castillo, & Emma Sherratt
#'
#' Functions in this package allow one to read, manipulate, and digitize landmark data; generate shape
#'  variables via Procrustes analysis for points, curves and surface data, perform statistical analyses
#'  of shape variation and covariation, and provide graphical depictions of shapes and patterns of
#'  shape variation.
#' 
#' @import rgl 
#' @importFrom geiger sim.char
#' @importFrom jpeg readJPEG
#' @importFrom phytools fastAnc
#' @importFrom vegan adonis
#' @importFrom ape is.binary.tree compute.brlen vcv.phylo summary.phylo Ntip

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
#' @keywords datasets
NULL

#' Head shape and food use data from Plethodon salamanders
#'
#' @name plethShapeFood
#' @docType data
#' @author Dean Adams
#' @references Adams, D. C., and F. J. Rohlf. 2000. Ecological character 
#' displacement in Plethodon: biomechanical differences found from a geometric 
#' morphometric study. Proceedings of the National Academy of Sciences, 
#' U.S.A. 97:4106-4111
#' @keywords datasets
NULL

#' Landmark data from dataset rat
#'
#' @name ratland
#' @docType data
#' @author Dean Adams
#' @references Bookstein, F. L. 1991. Morphometric tools for landmark data: Geometry and Biology. 
#'  Cambridge Univ. Press, New York.
#' @keywords datasets
NULL

#' Landmark data from hummingbird bills (includes sliding semilandmarks on curves)
#'
#' @name hummingbirds
#' @docType data
#' @author Chelsea Berns and Dean Adams
#' @references Berns, C.M., and Adams, D.C. 2010. Bill shape and sexual shape dimorphism between two species 
#' of temperate hummingbirds: Archilochus alexandri (black-chinned hummingbirds) and Archilochus colubris 
#' (ruby-throated hummingbirds). The Auk. 127:626-635.
#' @keywords datasets
NULL

#' Average head shape and phylogenetic relationships for several Plethodon salamander species
#'
#' @name plethspecies
#' @docType data
#' @author Dean Adams
#' @references Phylogeny pruned from: Wiens et al. (2006). Evol.
#' @references Data from: Adams and Rohlf (2000); Adams et al. (2007); Arif et al. (2007) Myers and Adams (2008)
#' @keywords datasets
NULL

#' Landmark data from scallop shells
#'
#' @name scallops
#' @docType data
#' @author Dean Adams and Erik Otarola-Castillo
#' @references Serb et al. (2011). "Morphological convergence of shell shape in distantly related
#' scallop species (Mollusca: Pectinidae)." Zoological Journal of the Linnean Society 163: 571-584.
#' @keywords datasets
NULL

#' 3D scan of a scallop shell from a .ply file in mesh3d format
#'
#' @name scallopPLY
#' @docType data
#' @author Emma Sherratt
#' @references Serb et al. (2011). "Morphological convergence of shell shape in distantly related
#' scallop species (Mollusca: Pectinidae)." Zoological Journal of the Linnean Society 163: 571-584.
#' @keywords datasets
NULL

#' Simulated motion paths
#'
#' @name motionpaths
#' @docType data
#' @author Dean Adams
#' @references Adams, D. C., and M. L. Collyer. 2009. A general framework for the analysis of phenotypic 
#'   trajectories in evolutionary studies. Evolution 63:1143-1154.
#' @keywords datasets
NULL

#' landmarks on mosquito wings
#'
#' @name mosquito
#' @docType data
#' @author Dean Adams
#' @keywords datasets
NULL

#### TPS and GPA routines (DCA and J Claude code) 

tps<-function(matr, matt, n,sz=1.5){		#DCA: altered from J. Claude: 2D only	
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
  points(matt,pch=21,bg="black",cex=sz)
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
#' @keywords utilities
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

TransRot<-function(M1,M2){  			#DCA: 4, 2014 translate, rotate one specimen to another
  k<-ncol(M1)
  Z1<-trans(M1)
  Z2<-trans(M2)
  sv<-svd(t(Z2)%*%Z1)
  U<-sv$v; V<-sv$u; Delt<-sv$d
  sig<-sign(det(t(Z2)%*%Z1))
  Delt[k]<-sig*abs(Delt[k]); V[,k]<-sig*V[,k]
  Gam<-U%*%t(V)
  list(Mp1=Z1%*%Gam,Mp2=Z2,rotation=Gam)}    


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

# picscale is called by digitize2d
picscale<- function(scale){
  digscale<-NULL
  digscale<-locator(2,type="o",lwd=2,col="red",lty="11")
  cat(paste("Keep scale (y/n)?"), "\n")
  ans <- readLines(n = 1)
  if (ans == "n") {
    cat(paste("Set scale again"), "\n")
  }
  while (ans == "n") {
    digscale<-NULL
    digscale<-locator(2,type="o",lwd=2,col="red",lty="11")
    cat(paste("Keep scale (y/n)?"), "\n")
    ans <- readLines(n = 1)
    if (ans == "y") { 
    }
    if (ans == "n") {
      cat(paste("Set scale again"), "\n")
    }
  }
  scale/sqrt(sum(diff(digscale$x)^2+diff(digscale$y)^2))      
}

# Function written by person who wrote identify() - called by define.modules
identifyPch <- function(x, y = NULL, n = length(x), pch = 19, col="red", ...)
{
  xy <- xy.coords(x, y); x <- xy$x; y <- xy$y
  sel <- rep(FALSE, length(x)); res <- integer(0)
  while(sum(sel) < n) {
    ans <- identify(x[!sel], y[!sel], n = 1, plot = FALSE, ...)
    if(!length(ans)) break
    ans <- which(!sel)[ans]
    points(x[ans], y[ans], pch = pch, col=col)
    sel[ans] <- TRUE
    res <- c(res, ans)
  }
  res
} 

# submodel prediction  (for procD.lm RRPP method)
submod.pred = function(X,Y){
	fit = lm(Y~X-1)
	B = as.matrix(coef(fit))
	if(nrow(B) == nrow(na.omit(B))) Yhat = X%*%B
	if(nrow(B) > nrow(na.omit(B))){
			B = na.omit(B)
			Yhat = (X)[,-attr(B, "na.action")]%*%B}
	Yhat
}

# residual randomization  (for procD.lm RRPP method)
RRP.submodels = function(Xs, Y){
	p = ncol(Y)
	n = nrow(Y)
	k = dim(Xs)[3]
	E = Yh = array(0,c(n,p,k))
	for(i in 1:k){
		yhat = submod.pred(Xs[,,i],Y)
		Yh[,,i] = yhat
		E[,,i] = Y-yhat
	}
	Er = E[sample(nrow(E)),,]
	Yr = Yh + Er
}

# P-values  (for procD.lm RRPP method)
pval = function(s){# s = sampling distribution
	p = length(s)
	r = rank(s)[1]-1
	pv = 1-r/p
	pv
}

#P value matrix  (for procD.lm RRPP method)
Pval.matrix = function(M){
	P = matrix(0,dim(M)[1],dim(M)[2])
	for(i in 1:dim(M)[1]){
		for(j in 1:dim(M)[2]){
			y = M[i,j,]
			p = pval(y)
			P[i,j]=p
		}
	}
	if(dim(M)[1] > 1 && dim(M)[2] >1) diag(P)=1
	rownames(P) = dimnames(M)[[1]]
	colnames(P) = dimnames(M)[[2]]
	P
}

# SSE via squared distances  (for procD.lm RRPP method)
SSE = function(L){# L is a linear model
	r = as.matrix(resid(L))
	S = r%*%t(r)
	sse = sum(diag(S))
	sse
}

##Support functions for pairwise.slope.test
vec.cor <- function(v1, v2){
  v1 = matrix(v1)
  v2 = matrix(v2)
  d1 = sqrt(t(v1)%*%v1)
  d2 = sqrt(t(v2)%*%v2)
  r =t(v1%*%solve(d1))%*%v2%*%solve(d2)
  r
}


vec.cor.matrix <- function(M) {
  M= as.matrix(M)
  w = solve(diag(sqrt(diag(M%*%t(M)))))
                 w = solve(diag(w))
                 z = w%*%as.matrix(M)
                 vc = z%*%t(z)
                 vc
}

vec.ang.matrix <- function(M){
  M= as.matrix(M)
  w = sqrt(diag(M%*%t(M)))
  w = solve(diag(w))
  z = w%*%as.matrix(M)
  vc = z%*%t(z)
  options(warn = -1)
  vc = acos(vc)
  diag(vc) = 0
  vc
}

# PLS calculations for two.b.pls analysis

pls = function(x,y){ # x and y must be vectors or matrices
  px <- ncol(x)
  py <- ncol(y)
  XY.vcv <- var(cbind(x, y))
  S12 <- XY.vcv[1:px, (px + 1):(px + py)]
  pls <- svd(S12)
  U <- pls$u
  V <- pls$v
  if(px && py == 1) {
    XScores <- x 
    YScores <- y 
  }
  if(px > 1 && py > 1) {
    XScores <- x %*% U
    YScores <- y %*% V
  }
  if(px == 1 && py > 1){
    XScores <- x %*% V
    YScores <- y %*% U
  } 
  if(px > 1 && py == 1) {
    XScores <- x %*% U
    YScores <- y %*% V
  }
  
  r <- cor(XScores[, 1], YScores[, 1])
  list(r=r, XScores = matrix(XScores[,1]), YScores = matrix(YScores[,1]))
}


