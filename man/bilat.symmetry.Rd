% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bilat.symmetry.r
\name{bilat.symmetry}
\alias{bilat.symmetry}
\title{Analysis of bilateral symmetry}
\usage{
bilat.symmetry(
  A,
  ind = NULL,
  side = NULL,
  replicate = NULL,
  object.sym = FALSE,
  land.pairs = NULL,
  data = NULL,
  iter = 999,
  seed = NULL,
  RRPP = TRUE,
  SS.type = c("I", "II", "III"),
  turbo = TRUE,
  Parallel = FALSE,
  print.progress = TRUE,
  ...
)
}
\arguments{
\item{A}{One of either A 3D array (p x k x n) containing raw landmarks (requiring GPA to be 
performed) or a gpagen object (if GPA has been previously performed) or a geomorphShapes object 
(requiring GPA to be performed).  Any gpagen argument should work within bilat.symmetry.}

\item{ind}{A vector containing labels for each individual. For matching symmetry, the matched 
pairs receive the same label (replicates also receive the same label).}

\item{side}{An optional vector (for matching symmetry) designating which object belongs to which
'side-group'}

\item{replicate}{An optional vector designating which objects belong to which group of replicates.
Alternatively, this can be a character value to indicate the name of the variable in the data frame to use.}

\item{object.sym}{A logical value specifying whether the analysis should proceed based on object 
symmetry = TRUE or matching symmetry = FALSE}

\item{land.pairs}{An optional matrix (for object symmetry) containing numbers for matched pairs 
of landmarks across the line of symmetry}

\item{data}{A data frame for the function environment, see \code{\link{geomorph.data.frame}}. It 
is imperative that the variables "ind", "side", and "replicate" in the data frame match these 
names exactly (as shown in examples below).}

\item{iter}{Number of iterations for significance testing.}

\item{seed}{An optional argument for setting the seed for random permutations of the resampling 
procedure. If left NULL (the default), the exact same P-values will be found for repeated runs 
of the analysis (with the same number of iterations). If seed = "random", a random seed will be 
used, and P-values will vary.  One can also specify an integer for specific seed values,
which might be of interest for advanced users.}

\item{RRPP}{A logical value indicating whether residual randomization should be used for 
significance testing.}

\item{SS.type}{A choice between type I (sequential), type II (hierarchical), or type III (marginal).}

\item{turbo}{A logical value that if TRUE, suppresses coefficient estimation in every random permutation, 
in order to speed up computation time.}

\item{Parallel}{Either a logical value to indicate whether parallel processing 
should be used or a numeric value to indicate the number of cores to use in 
parallel processing via the \code{parallel} library. 
If TRUE, this argument invokes forking of all processor cores, except one.  If
FALSE, only one core is used. A numeric value directs the number of cores to use,
but one core will always be spared.}

\item{print.progress}{A logical value to indicate whether a progress bar should be printed to 
the screen.  
This is helpful for long-running analyses.}

\item{...}{Arguments to pass onto gpagen}
}
\value{
An object of class "bilat.symmetry" returns a list of the following
\item{shape.anova}{An analysis of variance table for the shape data.}
\item{size.anova}{An analysis of variance table for the shape data (when object.sym = FALSE).}
\item{symm.shape}{The symmetric component of shape variation.}
\item{asymm.shape}{The asymmetric component of shape variation.}
\item{DA.component}{The directional asymmetry component, found as the mean shape for each side.}
\item{FA.component}{The fluctuating asymmetry component for each specimen, 
found as the specimen specific side deviation adjusted for the mean 
directional asymmetry in the dataset.}
\item{unsigned.AI}{Individual unsigned asymmetry index, as per Klingenberg and McIntyre, 1998;
Lazić et al 2015.}
\item{data.type}{A value indicating whether the analysis was performed as Object or Matching 
symmetry.}
\item{permutations}{The number of random permutations used.}
\item{random.shape.F}{A matrix of random F-values from the Shape analysis.}
\item{random.size.F}{A matrix of random F-values from the Centroid Size analysis (when 
object.sym = FALSE).}
\item{perm.method}{A value indicating whether "Raw" values were shuffled or "RRPP" performed.}
\item{procD.lm.shape}{A list of typical output from an object of class procD.lm, for shape}
\item{procD.lm.size}{If applicable, a list of typical output from an object of class procD.lm, 
for size (when object.sym = FALSE).}
\item{call}{The matched call.}
}
\description{
Function performs an analysis of directional and fluctuating asymmetry for bilaterally symmetric objects
}
\details{
The function quantifies components of shape variation for a set of specimens as described by 
their patterns of symmetry and asymmetry. Here, shape variation is decomposed into variation 
among individuals, variation among sides (directional asymmetry), and variation due to an 
individual x side interaction (fluctuating symmetry). These components are then statistically 
evaluated using Procrustes ANOVA. Statistical assessment of model effects for shape variation 
is accomplished using permutation procedures. Methods for both matching symmetry and object 
symmetry can be implemented. Matching symmetry is when each object contains mirrored pairs of 
structures (e.g., right and left hands) while object symmetry is when a single object is 
symmetric about a midline (e.g., right and left sides of human faces). Details on general 
approaches for the study of symmetry in geometric morphometrics may be found in: Mardia et 
al. 2000; Klingenberg et al. 2002. 
 
As input, the function receives either A 3D array (p x k x n) containing raw landmarks (requiring 
GPA to be performed) or a gpagen object (if GPA has been previously performed) or a geomorphShapes 
object. If one wishes to incorporate semilandmarks, GPA can either be performed first using gpagen,
or within bilat.symmetry by passing adequate GPA arguments (i.e. curves, surfaces, ProcD etc, 
see \code{\link{gpagen}}. If a geomorphShapes object is provided, semilandmarks are automatically 
identified and slid during GPA. For object.sym = FALSE, landmarks should be of dimension (p x k 
x 2n), as each specimen is represented by both left and right configurations.
   
Analyses of symmetry for matched pairs of objects is implemented when object.sym = FALSE. Here, 
a 3D array [p x k x 2n] contains the landmark coordinates for all pairs of structures (2 
structures for each of n specimens). Because the two sets of structures are on opposite sides,
they represent mirror images, and one set must be reflected prior to the analysis to allow 
landmark correspondence. IT IS ASSUMED THAT THE USER HAS DONE THIS PRIOR TO PERFORMING THE
SYMMETRY ANALYSIS. Reflecting a set of specimens may be accomplished by multiplying one coordinate 
dimension by '-1' for these structures (either the x-, the y-, or the z-dimension). A vector 
containing information on individuals and sides must also be supplied. Replicates of each 
specimen may also be included in the dataset, and when specified will be used as measurement 
error (see Klingenberg and McIntyre 1998). 

Analyses of object symmetry is implemented when object.sym = TRUE. Here, a 3D array [p x k x n] 
contains the landmark coordinates for all n specimens. To obtain information about asymmetry, 
the function generates a second set of objects by reflecting them about one of their coordinate 
axes. The landmarks across the line of symmetry are then relabeled to obtain landmark 
correspondence. The user must supply a list of landmark pairs. A vector containing information 
on individuals must also be supplied. Replicates of each specimen may also be included in the 
dataset, and when specified will be used as measurement error. 

The function also provides individual measures of unsigned asymmetry, calculated as the
Procrustes distance between the right and left element (for paired structures, as detailed in 
Klingenberg and McIntyre 1998) or side of the structure (for object symmetry, following Lazić 
et al 2015). The computational difference between the two approaches consists in that, for object
symmetry, only paired landmarks are considered, excluding the landmarks of the midline.
 
\subsection{Notes for geomorph 3.0}{ 
Compared to older versions of geomorph, some results can be expected to be slightly different. 
Starting with geomorph 3.0, results use only type I sums of squares (SS) with either full 
randomization of raw shape values or RRPP (preferred with nested terms) for analysis of variance
(ANOVA).  Older versions used a combination of parametric and non-parametric results, as well as 
a combination of type I and type III SS.  While analytical conclusions should be consistent 
(i.e., "significance" of effects is the same), these updates maintain consistency in analytical 
philosophy.  This change will require longer computation time for large datasets, but the 
trade-off allows users to have more flexibility and eliminates combining disparate analytical 
philosophies. 
 
Note also that significance of terms in the model are found by comparing F-values for each term 
to those obtained via permutation.  F-ratios and df are not strictly necessary (a ratio of SS 
would suffice), but they are reported as is standard for anova tables. Additionally, users will 
notice that the df reported are based on the number of observations rather than a combination 
of objects * coordinates * dimensions, as is sometimes found in morphometric studies of symmetry. 
However, this change has no effect on hypothesis testing, as only SS vary among permutations (df, 
coordinates, and dimensions are constants). 
}
 
The generic functions, \code{\link{print}}, \code{\link{summary}}, and \code{\link{plot}} all work 
with \code{\link{bilat.symmetry}}.
}
\examples{
\dontrun{

#Example of matching symmetry
data(mosquito)
gdf <- geomorph.data.frame(wingshape = mosquito$wingshape, 
ind = mosquito$ind, 
side = mosquito$side,
replicate = mosquito$replicate)
mosquito.sym <- bilat.symmetry(A = wingshape, ind = ind, side = side,
replicate = replicate, object.sym = FALSE, RRPP = TRUE, 
data = gdf)
summary(mosquito.sym)
plot(mosquito.sym, warpgrids = TRUE)
mosquito.sym$shape.anova # extract just the anova table on shape

# Previous example, performing GPA first
Y.gpa <- gpagen(mosquito$wingshape)
mosquito.sym2 <- bilat.symmetry(A = Y.gpa, ind = ind, side = side,
replicate = replicate, object.sym = FALSE, RRPP = TRUE, 
data = gdf)
summary(mosquito.sym2)
summary(mosquito.sym) # same results

#Example of object symmetry

data(lizards)
gdf <- geomorph.data.frame(shape = lizards$coords, 
ind = lizards$ind, 
replicate = lizards$rep)
liz.sym <- bilat.symmetry(A = shape, ind = ind, rep = rep, 
object.sym = TRUE, 
land.pairs = lizards$lm.pairs, data = gdf, RRPP = TRUE)
summary(liz.sym)

# Example of object symmetry in 3D and including semilandmarks

data(scallops)
gdf <- geomorph.data.frame(shape = scallops$coorddata, 
ind = scallops$ind)
scallop.sym <- bilat.symmetry(A = shape, ind = ind, 
object.sym = TRUE, 
curves= scallops$curvslide, surfaces = scallops$surfslide,
land.pairs=scallops$land.pairs, data = gdf, RRPP = TRUE)
summary(scallop.sym)
# NOTE one can also: plot(scallop.sym, warpgrids = TRUE, mesh = NULL)
# NOTE one can also: scallop.sym$data.type # recall the symmetry type
}
}
\references{
Klingenberg, C.P. and G.S. McIntyre. 1998. Quantitative genetics of geometric shape 
in the mouse mandible. Evolution. 55:2342-2352.

Mardia, K.V., F.L. Bookstein, and I.J. Moreton. 2000. Statistical assessment of 
bilateral symmetry of shapes. Biometrika. 87:285-300.

Klingenberg, C.P., M. Barluenga, and A. Meyer. 2002. Shape analysis of symmetric 
structures: quantifying variation among individuals and asymmetry. Evolution. 56:1909-1920.

Lazić, M. M., M. A. Carretero, J. Crnobrnja-Isailović, and A. Kaliontzopoulou. 2015.
Effects of environmental disturbance on phenotypic variation: an integrated assessment of 
canalization, developmental stability, modularity, and allometry in lizard head shape. The American
Naturalist 185:44–58.
}
\author{
Dean Adams, Michael Collyer and Antigoni Kaliontzopoulou
}
\keyword{analysis}
