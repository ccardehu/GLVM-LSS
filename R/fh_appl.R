rm(list = ls())
set.seed(1234)
# rm(.Random.seed, envir=globalenv())

library(MASS)
library(abind)
library(splines)
library(mvtnorm)
library(numDeriv)
library(gamlss)
library(haven)

# `%!in%` <- Negate(`%in%`)

source("fa_main.R")
source("fb_auxf.R")
source("fc_mfit.R")
source("fd_mghq.R")
source("ff_msim.R")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# Example No.1 : European Social Survey - ESS #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

trPOL <- read_dta("C:/Users/carde/Dropbox/Camilo and Irini/Research/Data/ESS/Immigration/trPOL.dta")
trPOL <- as.data.frame(trPOL[complete.cases(trPOL),-c(1,ncol(trPOL))])
summary(trPOL)

trPOL <- trPOL[,1:7]; #colnames(trPOL) <- paste0("h",1:7)

p. <- ncol(trPOL)
famd <- c(rep("normal", p.))#, rep("binomial", 3))
l01. <- l11. <- l21. <- l31. <- NULL
formd0 <- NULL; formd0$mu <- "~ Z1"; formd0$sigma <- "~ 1"
formd1 <- NULL; formd1$mu <- "~ Z1"; formd1$sigma <- "~ Z1"
formd2 <- NULL; formd2$mu <- "~ Z1 + I(Z1^2)"; formd2$sigma <- "~ 1"
formd3 <- NULL; formd3$mu <- "~ Z1 + I(Z1^2) "; formd3$sigma <- "~ Z1"
l01.$mu <- matrix(1,ncol = 2, nrow = p.)
l01.$sigma <- matrix(1, ncol = 1, nrow = p.)
l11.$mu <- matrix(1,ncol = 2, nrow = p.)
l11.$sigma <- matrix(1, ncol = 2, nrow = p.)
l21.$mu <- matrix(1,ncol = 3, nrow = p.)
l21.$sigma <- matrix(1, ncol = 1, nrow = p.)
l31.$mu <- matrix(1,ncol = 3, nrow = p.)
l31.$sigma <- matrix(1, ncol = 2, nrow = p.)

# Linear + Homoscedastic Factor Model
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

test0 <- splvm.fit(trPOL,famd,formd0,
          control = list(method = "EM", constraint = l01., full.hess = F,
          ghQqp = 30, iter.lim = 1e3, tol = sqrt(.Machine$double.eps), silent = F,information = "Fisher"))

# test00 <- NULL; test00$b$mu <- test0$b$mu[1:7,]; test00$b$sigma <- test0$b$sigma[1:7,,drop=F]
# syntax = 'Z1 =~ h1 + h2 + h3 + h4 + h5 + h6 + h7'
# testmod <- lavaan::cfa(syntax, data = trPOL, orthogonal = T, meanstructure = T, std.lv = TRUE, estimator="ML",
#                information = "observed"); summary(testmod)
# testcoef <- mb2lb(matrix(coef(testmod),nrow = ncol(trPOL),byrow = F)[,c(3,1,2)],test00$b)
# testcoef$sigma <- log(sqrt(testcoef$sigma))

GBIC(test0)
test0$loglik
test0$b

# Linear + Heteroscedastic Factor Model
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

test1 <- splvm.fit(trPOL,famd,formd1,
          control = list(method = "EM", constraint = l11., full.hess = F,
          ghQqp = 30, iter.lim = 1e3, tol = sqrt(.Machine$double.eps), silent = F,information = "Fisher"))
GBIC(test1)
test1$loglik
test1$b

test1p <- splvm.fit(trPOL,famd,formd1,
          control = list(method = "PEM", constraint = l11., full.hess = F, start.val = test1$b,
          ghQqp = 30, iter.lim = 1e3, tol = sqrt(.Machine$double.eps), silent = F,information = "Fisher",
          pml.control = list(type = "mcp", w.alasso = test1$b, lambda = 0.1, pen.load = T)))
GBIC(test1p)
test1p$loglik
test1p$b

# Quadratic + Homoscedastic Factor Model
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

test2 <- splvm.fit(trPOL,famd,formd2,
          control = list(method = "EM", constraint = l21., full.hess = F,
          ghQqp = 30, iter.lim = 1e3, tol = sqrt(.Machine$double.eps), silent = F,information = "Fisher"))
GBIC(test2)
test2$loglik
test2$b

test2p <- splvm.fit(trPOL,famd,formd2,
          control = list(method = "PEM", constraint = l21., full.hess = F, start.val = test2$b,
          ghQqp = 30, iter.lim = 1e3, tol = sqrt(.Machine$double.eps), silent = F,information = "Fisher",
          pml.control = list(type = "mcp", w.alasso = test2$b, lambda = 0.2, pen.load = T)))
GBIC(test2p)
test2p$loglik
test2p$b

# Quadratic + Heteroscedastic Factor Model
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

test3 <- splvm.fit(trPOL,famd,formd3,
          control = list(method = "EM", constraint = l31., full.hess = F,
          ghQqp = 30, iter.lim = 1e3, tol = sqrt(.Machine$double.eps), silent = F,information = "Fisher"))
GBIC(test3)
test3$loglik
test3$b

test3p <- splvm.fit(trPOL,famd,formd3,
          control = list(method = "PEM", constraint = l31., full.hess = F, start.val = test3$b,
          ghQqp = 30, iter.lim = 1e3, tol = sqrt(.Machine$double.eps), silent = F,information = "Fisher",
          pml.control = list(type = "mcp", w.alasso = test3$b, lambda = 0.2, pen.load = T)))
GBIC(test3p)
test3p$loglik
test3p$b

round(c(GBIC(test0),GBIC(test1),GBIC(test1p),GBIC(test2),GBIC(test2p),GBIC(test3),GBIC(test3p)),3)
which.min(c(GBIC(test0),GBIC(test1),GBIC(test1p),GBIC(test2),GBIC(test2p),GBIC(test3),GBIC(test3p)))
round(c(test0$loglik,test1$loglik,test1p$loglik,test2$loglik,test2p$loglik,test3$loglik,test3p$loglik),3)
round(cbind(test0$b$mu,test1$b$mu,test1p$b$mu,test2$b$mu,test2p$b$mu,test3$b$mu,test3p$b$mu),3)
round(cbind(test0$b$sigma,test1$b$sigma,test1p$b$sigma,test2$b$sigma,test2p$b$sigma,test3$b$sigma,test3p$b$sigma),3)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# Example No.1 : European Social Survey - ESS #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
 
