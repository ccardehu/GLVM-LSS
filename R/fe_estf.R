
rm(list = ls())
set.seed(1234)
# rm(.Random.seed, envir=globalenv())

library(mvtnorm)
library(MASS)
library(numDeriv)
library(abind)

`%!in%` <- Negate(`%in%`)

source("fa_main.R")
source("fb_auxf.R")
source("fc_mfit.R")
source("SimFA.R")
source("GHFun.R")
# source("GLVMfit.R")

n = 1000     # Number of individuals
p = 10       # Number of items
nsim = 1000  # Number of simulations
form <- list("mu" = "~ Z1 + I(Z1^2) + I(cos(Z1))", "sigma" = "~ Z1")
fam <- rep("normal",p)

l1 <- NULL
l1$mu <- matrix(1,ncol = 4, nrow = p)
l1$sigma <- matrix(1, ncol = 2, nrow = p)

lc <- NULL
lc$mu <- matrix(runif(length(l1$mu), min = -0.5, max = 0.5),nrow = p)
lc$mu[,1] <- runif(p,1,2)
lc$mu[,2] <- runif(p,1.0,3.0)
lc$mu[,3] <- runif(p,-1.0,-0.5)
# lc$mu[,4] <- runif(p,-3.0,-1.0)
lc$sigma <- matrix(runif(length(l1$sigma), min = 0.1, max = 0.3), nrow = p)

simR <- simGLVM(n = n, p = p, form = form, dist = fam, loadmt = l1, coefs = lc)
Y <- simR$Y
Z <- simR$Z
borg <- simR$borg

testa <- splvm.fit(Y,fam,form,
          control = list(method = "EM", constraint = l1, full.hess = F, start.val = borg, 
          ghQqp = 25, iter.lim = 250, tol = sqrt(.Machine$double.eps), silent = F))

testb <- splvm.fit(Y,fam,form,
          control = list(method = "EM", constraint = l1, full.hess = T, start.val = borg, 
          ghQqp = 25, iter.lim = 250, tol = sqrt(.Machine$double.eps), silent = F))

testc <- splvm.fit(Y,fam,form,
          control = list(method = "ML", constraint = l1, start.val = borg,
          ghQqp = 25, iter.lim = 250, tol = sqrt(.Machine$double.eps), silent = F))

testd <- splvm.fit(Y,fam,form,
          control = list(method = "hybrid", constraint = l1, full.hess = F, start.val = borg,
          ghQqp = 10, iter.lim = 250, tol = sqrt(.Machine$double.eps), silent = F))

cbind(testa$b$mu, testb$b$mu, testc$b$mu, testd$b$mu, borg$mu)

rbenchmark::benchmark(
  "new" = {testa <- splvm.fit(Y,fam,form,
          control = list(method = "EM", constraint = l1, start.val = borg, full.hess = F,
          ghQqp = 20, iter.lim = 350, tol = sqrt(.Machine$double.eps), silent = F))},
  
  "old" = {test2 <- GLVM.fit3(Y = Y, fam = fam, form = form , silent = F, ghp = 20, iter.lim = 1000,
                  tol = sqrt(.Machine$double.eps), loadmt = l1, icoefs = lc, useoptim = F, skipEM = F)},
  replications = 1,
  columns = c("test", "replications", "elapsed",
              "relative", "user.self", "sys.self")
)
