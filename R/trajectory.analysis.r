#'  Quantify and compare shape change trajectories
#'
#'  Function estimates attributes of shape change trajectories or motion trajectories for a set of 
#'  Procrustes-aligned specimens and compares them statistically
#'
#'  The function quantifies phenotypic shape change trajectories from a set of specimens, and assesses variation 
#'  in these parameters via permutation. A shape change trajectory is defined by a sequence 
#'  of shapes in tangent space. These trajectories can be quantified various attributes (their size, orientation, 
#'  and shape), and comparisons of these attribute enables the statistical comparison of shape change 
#'  trajectories (see Collyer and Adams 2013; Collyer and Adams 2007; Adams and Collyer 2007; Adams and Collyer 2009). 
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
#' @keywords analysis
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
  if(any(is.na(y))==T){
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
    traj.specs.obs<- aperm(array(t(lsmeans.obs), c(k1,p1,n1)), c(2,1,3)) 
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
      traj.specs.r<- aperm(array(t(lsmeans.r), c(k1,p1,n1)), c(2,1,3)) 
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
    traj.specs.obs<-arrayspecs(y,p1,k1) 
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