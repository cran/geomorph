#' Phylogenetic ANOVA/regression for shape data
#'
#' Function performs Procrustes ANOVA in a phylogenetic framework and uses  permutation procedures to assess 
#' statistical hypotheses describing patterns of shape variation and covariation for a set of Procrustes-aligned coordinates
#'
#' The function performs ANOVA and regression models in a phylogenetic context under a Brownian motion model of evolution, 
#' in a manner that can accommodate 
#' high-dimensional datasets. The approach is derived from the statistical equivalency between parametric methods 
#' utilizing covariance matrices and methods based on distance matrices (Adams 2014). Data input is specified by 
#' a formula (e.g., y~X), where 'y' specifies the response variables (shape data), and 'X' contains one or more 
#' independent variables (discrete or continuous). The response matrix 'y' must be in the form of a two-dimensional data 
#'   matrix of dimension (n x [p x k]), rather than a 3D array.  It is assumed that the landmarks have previously 
#'   been aligned using Generalized Procrustes Analysis (GPA) [e.g., with \code{\link{gpagen}}].
#'   The user must also specify a phylogeny describing the evolutionary relationships among species (of class phylo).
#'   Note that the specimen labels for both x and y must match the labels on the tips of the phylogeny.
#'
#'   From the phylogeny, a phylogenetic transformation matrix is obtained under a Brownian motion model, and used to 
#'   transform the x and y variables. Next, the Gower-centered distance matrix is obtained from predicted values from the
#'   model (y~x), from which sums-of-squares, F-ratios, and R^2 are estimated for each factor in the model (see Adams, 2014). 
#'   Data are then permuted across the tips of the phylogeny, and estimates of statistical values are obtained for the permuted data,
#'   which are  compared to the observed value to assess significance. 
#'
#' @param f1 A formula for the linear model (e.g., y~x1+x2)
#' @param phy A phylogenetic tree of {class phylo} - see \code{\link[ape]{read.tree}} in library ape
#' @param iter Number of iterations for significance testing
#' @keywords analysis
#' @export
#' @author Dean Adams
#' @return Function returns an ANOVA table of statistical results for all factors: df (for each factor), SS, MS,
#' F ratio, Prand, and Rsquare.
#' @references Adams, D.C. 2014. A method for assessing phylogenetic least squares models for shape and other high-dimensional 
#' multivariate data. Evolution. 68. DOI:10.1111/evo.12463. 
#' @examples
#' ### Example of D-PGLS for high-dimensional data 
#' data(plethspecies)
#' Y.gpa<-gpagen(plethspecies$land)    #GPA-alignment
#' 
#' procD.pgls(two.d.array(Y.gpa$coords)~Y.gpa$Csize,plethspecies$phy,iter=49)
procD.pgls<-function(f1,phy,iter=999){
    data=NULL
    form.in<-formula(f1)
    Terms<-terms(form.in,keep.order=TRUE)
    Y<-as.matrix(eval(form.in[[2]],parent.frame()))
    N<-length(phy$tip.label)
    p<-ncol(Y)
    if(is.null(rownames(Y))){
      stop("No species names with Y-data.")  }
    if(length(match(rownames(Y), phy$tip.label))!=N) 
      stop("Data matrix missing some taxa present on the tree.")
    if(length(match(phy$tip.label,rownames(Y)))!=N) 
      stop("Tree missing some taxa in the data matrix.")
    if (any(is.na(match(sort(phy$tip.label), sort(rownames(Y)))) == T)) {
      stop("Names do not match between tree and data matrix.")}

    C<-vcv.phylo(phy); C<-C[rownames(Y),rownames(Y)]  
    eigC <- eigen(C)
    D.mat<-solve(eigC$vectors %*% diag(sqrt(eigC$values)) %*% t(eigC$vectors)) 
    Y.new<-D.mat %*% (Y)    
    ones.new<-D.mat%*%(array(1,N))
    pred.1<- predict(lm(Y.new~ones.new-1)) 
    dat<-model.frame(form.in,data)
    df<-df.tmp<-SS.tmp<-SS.obs<-F<-array()
    for (i in 1:ncol(attr(Terms, "factors"))){
      mod.mat<-model.matrix(Terms[1:i],data=dat)
      x.new<-D.mat%*%mod.mat
      pred.y<-predict(lm(Y.new~x.new-1))
      G<-(pred.y-pred.1)%*%t(pred.y-pred.1)
      SS.tmp[i]<-sum(diag(G))   
      ifelse(i==1, SS.obs[i]<-SS.tmp[i], SS.obs[i]<-(SS.tmp[i]-SS.tmp[i-1]))
      df.tmp[i]<-ifelse(ncol(mod.mat)==1,1,(ncol(mod.mat)-1))
      ifelse(i==1, df[i]<-df.tmp[i], df[i]<-(df.tmp[i]-df.tmp[i-1]))
    }
    MS<-SS.obs/df
    mod.mat<-model.matrix(Terms)
    x.new<-D.mat%*%mod.mat
    y.res<-residuals(lm(Y.new~x.new-1))
    SS.res<-sum(diag(y.res%*%t(y.res)))  
    df.res<-nrow(Y)-1-sum(df)
    MS.res<-SS.res/df.res
    Rsq<-SS.obs/(sum(SS.obs)+SS.res)
    F<-MS/MS.res
    F.r<-P.val<-array(1,dim=length(SS.obs))
    for(i in 1:iter){
      SS.tmp<-SS.r<-array()
      Y.r<-as.matrix(Y[sample(nrow(Y)),])  
      row.names(Y.r)<-row.names(Y)
      Y.r.new<-D.mat %*% (Y.r)    
      pred.1.r<- predict(lm(Y.r.new~ones.new-1)) 
      for (ii in 1:ncol(attr(Terms, "factors"))){
        mod.mat<-model.matrix(Terms[1:ii])
        x.new<-D.mat%*%mod.mat
        pred.y.r<-predict(lm(Y.r.new~x.new-1))
        G.r<-(pred.y.r-pred.1.r)%*%t(pred.y.r-pred.1.r)
        SS.tmp[ii]<-sum(diag(G.r))   
        ifelse(ii==1, SS.r[ii]<-SS.tmp[ii], SS.r[ii]<-(SS.tmp[ii]-SS.tmp[ii-1]))
      }
      MS.r<-SS.r/df
      mod.mat<-model.matrix(Terms)
      x.new<-D.mat%*%mod.mat
      y.res.r<-residuals(lm(Y.r.new~x.new-1))
      SS.r.res<-sum(diag(y.res.r%*%t(y.res.r)))  
      MS.r.res<-SS.r.res/df.res
      F.r<-MS.r/MS.r.res
      P.val<-ifelse(F.r>=F, P.val+1,P.val) 
    }
    P.val<-P.val/(iter+1)
    
    anova.tab <- data.frame(df = c(df,df.res), 
    SS = c(SS.obs, SS.res), 
    MS = c(MS, MS.res),
    Rsq = c(Rsq, NA),
    F = c(F, NA),
    P.val = c(P.val, NA))
    rownames(anova.tab) <- c(attr(Terms, "term.labels"), "Residuals")
    class(anova.tab) <- c("anova", class(anova.tab))
    return(anova.tab)
  }
