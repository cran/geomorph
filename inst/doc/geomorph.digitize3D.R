## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ----eval=TRUE----------------------------------------------------------------
library(geomorph)
data("scallopPLY")
my.ply <- scallopPLY$ply

## ----eval=FALSE---------------------------------------------------------------
# fixed.lms1 <- digit.fixed(spec = my.ply, fixed = 5)

## ----echo = FALSE, out.width="80%"--------------------------------------------
knitr::include_graphics("figs/fixed.3D.png")  

## ----eval=FALSE---------------------------------------------------------------
# my.ply.2 <- scallopPLY$ply
# fixed.lms2 <- digit.fixed(my.ply.2, 5)

## ----eval=FALSE---------------------------------------------------------------
# surf.pts1 <- buildtemplate(spec = my.ply, fixed = fixed.lms1, surface.sliders = 100)

## ----echo = FALSE, out.width="80%"--------------------------------------------
knitr::include_graphics("figs/example_template.png")  

## ----eval=FALSE---------------------------------------------------------------
# surf.pts2 <- digitsurface(spec = my.ply.2, fixed = fixed.lms2)

## ----echo = FALSE, out.width="80%"--------------------------------------------
knitr::include_graphics("figs/example_spec2.png")  

