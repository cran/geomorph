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
#'   Landmark groups can be defined using \code{\link{define.modules}}. 
#'
#' @param A An array (p x k x n) containing GPA-aligned coordinates for all specimens
#' @param landgroups A list of which landmarks belong in which partition (e.g. A,A,A,B,B,B,C,C,C)
#' @param iter Number of iterations for significance testing
#' @export
#' @keywords analysis
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
    stop("Data matrix not a 3D array (see 'arrayspecs').")  }
  if(any(is.na(A))==T){
    stop("Data matrix contains missing values. Estimate these first (see 'estimate.missing').")}
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