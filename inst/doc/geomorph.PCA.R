## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ----eval=TRUE----------------------------------------------------------------
library(geomorph)
data(plethspecies) 
Y.gpa <- gpagen(plethspecies$land, print.progress = F)

## ----eval=TRUE, fig.width = 6, fig.height = 6---------------------------------
PCA <- gm.prcomp(Y.gpa$coords)
summary(PCA)
plot(PCA, main = "PCA")

## ----eval=TRUE, fig.width = 6, fig.height = 6---------------------------------
plot(PCA, main = "PCA", pch = 22, bg = "green", cex = 1.5, cex.lab = 1.5, font.lab = 2)

## ----eval=TRUE, fig.width = 6, fig.height = 3---------------------------------
msh <- mshape(Y.gpa$coords)
plotRefToTarget(PCA$shapes$shapes.comp1$min, msh)
plotRefToTarget(msh, PCA$shapes$shapes.comp1$max)

## ----eval=TRUE, fig.width = 6, fig.height = 3---------------------------------
plotRefToTarget(PCA$shapes$shapes.comp1$min, PCA$shapes$shapes.comp1$max, method = "vector", mag = 2)

## ----eval=FALSE---------------------------------------------------------------
#  picknplot.shape(plot(PCA), method = "vector")

## ----eval=TRUE, fig.width = 6, fig.height = 6---------------------------------
PCA.w.phylo <- gm.prcomp(Y.gpa$coords, phy = plethspecies$phy)
summary(PCA.w.phylo)
plot(PCA.w.phylo, phylo = TRUE, main = "PCA.w.phylo")

## ----eval=TRUE, fig.width = 6, fig.height = 6---------------------------------

# Phylo PCA without projecting untransformed residuals
phylo.PCA <- gm.prcomp(Y.gpa$coords, phy = plethspecies$phy, GLS = TRUE)
summary(phylo.PCA)
plot(phylo.PCA, phylo = TRUE, main = "phylo PCA")


## ----eval=TRUE, fig.width = 6, fig.height = 6---------------------------------

# Phylo PCA without projecting transformed residuals
phylo.tPCA <- gm.prcomp(Y.gpa$coords, phy = plethspecies$phy, GLS = TRUE, transform = TRUE)
summary(phylo.tPCA)
plot(phylo.tPCA, phylo = TRUE, main = "phylo PCA with transformed projection")
par(mfrow = c(1, 1))


## ----eval=TRUE, fig.width = 6, fig.height = 6---------------------------------
PaCA <- gm.prcomp(Y.gpa$coords, phy = plethspecies$phy, align.to.phy = TRUE)
summary(PaCA)
plot(PaCA, phylo = TRUE, main = "PaCA")

## ----eval=FALSE---------------------------------------------------------------
#  plot(PCA.w.phylo, time.plot = TRUE, pch = 22, bg = c(rep("red", 5), rep("green", 4)), cex = 2,
#       phylo.par = list(edge.color = "grey60", edge.width = 1.5, tip.txt.cex = 0.75,
#                        node.labels = F, anc.states = F))
#  

## ---- echo = FALSE, out.width="40%"-------------------------------------------
knitr::include_graphics("figs/PCAplot2D.png")  
knitr::include_graphics("figs/tree.plot.png")  

