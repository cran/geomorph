## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ---- include=TRUE, fig.height=4, fig.width=5----------------------------
library(geomorph)
data("larvalMorph")

Y.gpa <- gpagen(larvalMorph$tailcoords, curves = larvalMorph$tail.sliders,
                ProcD = FALSE, print.progress = FALSE)
plot(Y.gpa)

## ---- include=TRUE-------------------------------------------------------
gdf <- geomorph.data.frame(Y.gpa, treatment = larvalMorph$treatment, 
                           family = larvalMorph$family)

fit.size <- procD.lm(coords ~ log(Csize), 
                     data = gdf, print.progress = FALSE) # simple allometry model
fit.family<- procD.lm(coords ~ log(Csize) * family, 
                     data = gdf, print.progress = FALSE) # unique family allometries
fit.treatment<- procD.lm(coords ~ log(Csize) * treatment/family, 
                     data = gdf, print.progress = FALSE) # unique treatment: family allometries

fit.size
fit.family
fit.treatment

## ---- include=TRUE-------------------------------------------------------
anova(fit.size)
anova(fit.family)
anova(fit.treatment)

## ---- include=TRUE-------------------------------------------------------
anova(fit.size, fit.family, fit.treatment, print.progress = FALSE)

## ---- include=TRUE-------------------------------------------------------
summary(fit.size)

## ---- include=TRUE, fig.height=4, fig.width=5----------------------------
plot(fit.size, type = "regression", reg.type = "PredLine", predictor = log(gdf$Csize))
plot(fit.size, type = "regression", reg.type = "RegScore", predictor = log(gdf$Csize))

## ---- include=TRUE, fig.height=4, fig.width=5----------------------------
plotAllometry(fit.size, size = gdf$Csize, logsz = TRUE, method = "PredLine")
plotAllometry(fit.size, size = gdf$Csize, logsz = TRUE, method = "RegScore")

## ---- include=TRUE, fig.height=4, fig.width=5----------------------------
PLS <- two.b.pls(log(gdf$Csize), gdf$coords, print.progress = FALSE)
PLS
plot(PLS)

## ---- include=TRUE, fig.height=4, fig.width=6----------------------------
plotAllometry(fit.size, size = gdf$Csize, logsz = TRUE, method = "CAC")

## ---- include=TRUE, fig.height=4, fig.width=5----------------------------
plotAllometry(fit.size, size = gdf$Csize, logsz = TRUE, method = "size.shape")

## ---- include = TRUE-----------------------------------------------------

fit.unique <- procD.lm(coords ~ log(Csize) * treatment/family, 
                     data = gdf, print.progress = FALSE) # unique allometries
fit.common <- procD.lm(coords ~ log(Csize) + treatment/family, 
                     data = gdf, print.progress = FALSE) # common allometry
anova(fit.common, fit.unique, print.progress = FALSE)


## ---- include=TRUE,  fig.height=4, fig.width=5---------------------------
plotAllometry(fit.common, size = gdf$Csize, logsz = TRUE, method = "PredLine",
              pch = 19, col = as.numeric(gdf$treatment))
plotAllometry(fit.common, size = gdf$Csize, logsz = TRUE, method = "RegScore",
              pch = 19, col = as.numeric(gdf$treatment))

## ---- include=TRUE-------------------------------------------------------
anova(fit.common)

## ---- include=TRUE-------------------------------------------------------
anova(fit.common, error = c("Residuals", "treatment:family", "Residuals"))

## ---- include=TRUE-------------------------------------------------------
reveal.model.designs(fit.common)

## ---- include=TRUE-------------------------------------------------------
fit.null <- procD.lm(coords ~ log(Csize) + family, data = gdf, print.progress = FALSE)
PW <- pairwise(fit.common, fit.null, groups = gdf$treatment, print.progress = FALSE)
PW

## ---- include=TRUE-------------------------------------------------------
summary(PW, test.type = "dist", confidence = 0.95)

## ---- include=TRUE-------------------------------------------------------
anova(fit.null, fit.common, print.progress = FALSE)
summary(PW, test.type = "dist", confidence = 0.95, stat.table = FALSE)

## ---- include=TRUE-------------------------------------------------------
summary(PW, test.type = "var", confidence = 0.95)

## ---- include=TRUE-------------------------------------------------------
morphol.disparity(fit.common, groups = gdf$treatment, print.progress = FALSE)

## ---- include=TRUE, fig.height=4, fig.width=5----------------------------
TA <- trajectory.analysis(fit.common, 
                          groups = gdf$treatment, traj.pts = gdf$family,
                          pca = TRUE, print.progress = FALSE)
summary(TA, attribute = "MD")
TP <- plot(TA, pch = 19, cex = 0.7, col = as.numeric(gdf$treatment))
add.trajectories(TP, traj.bg = 1:nlevels(gdf$treatment), 
                 start.bg = 1:nlevels(gdf$treatment),
                 end.bg = 1:nlevels(gdf$treatment))

