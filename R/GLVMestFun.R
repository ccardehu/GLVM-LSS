rm(list = ls())
#rm(list= ls()[!(ls() %in% c("ex1", "simR"))])
set.seed(1234)
# rm(.Random.seed, envir=globalenv())

library(mvtnorm)
library(MASS)
#library(profvis)
#library(ltm)
#library(gamlss.dist)
library(numDeriv)

`%!in%` <- Negate(`%in%`)

source("ddFun.R")
source("SimFA.R")
source("GHFun.R")
source("scFun.R")
source("GLVMfit.R")
source("graphFun.R")
source("EMFun.R")

n = 1000     # Number of individuals
p = 5       # Number of items
nsim = 1000  # Number of simulations
form <- list("mu" = "~ Z1")
#form1 <- list("mu" = "~ Z1 + I(Z1^2)", "sigma" = "~ 1") #  , "sigma" = "~ Z1"
#form2 <- list("mu" = "~ Z1", "sigma" = "~ 1") #  , "sigma" = "~ Z1"
fam <- rep("poisson",p)
# fam <- c(rep("normal", p/2), rep("poisson", p/2)) # rep("poisson",p)#
# sample(c("normal","poisson","binom"),size = 10,replace = T)

l1 <- NULL
l1$mu <- matrix(1,ncol = 2, nrow = p)
l1$sigma <- matrix(1, ncol = 1, nrow = p)

# Restrictions
# ____________
#
#l1$mu[6:10,2] <- 0
#l1$mu[1:5,3] <- 0
#l1$mu[sample(length(l1$mu), 20)] <- 0
#l1$sigma[sample(length(l1$sigma), 10)] <- 0

# l0 <- NULL
# l0$mu <- l1$mu * 0
# l0$sigma <- l1$sigma * 0

# l1 when model misspecification
# ______________________________
#
#l11 <- l1
#l11$sigma <- l11$sigma[,-c(2:ncol(l11$sigma)), drop = F]
#l12 <- l11
#l12$mu <- l12$mu[,-ncol(l12$mu), drop = F]

lc <- NULL
lc$mu <- matrix(runif(length(l1$mu), min = 0.1, max = 0.5),nrow = p)
#lc$mu[,1] <- runif(p,1,2)
#lc$mu[,2] <- runif(p,-0.5,0.5)
#lc$mu[,3] <- runif(p,-0.5,0.5)
#lc$mu[,4] <- runif(p,-0.5,-0.1)
lc$sigma <- matrix(runif(length(l1$sigma), min = 0.1, max = 0.5), nrow = p)
#lc$sigma[,1] <- runif(p,-1,1)
#lc$sigma[,2] <- runif(p,2,4)
#lc$sigma[,3] <- runif(p,-0.5,-0.1)


# lc when model misspecification
# ______________________________
#
#lc1 <- lc
#lc1$sigma <- lc1$sigma[,-c(2:ncol(lc1$sigma)), drop = F]
#lc2 <- lc1
#lc2$mu <- lc2$mu[,-ncol(lc2$mu), drop = F]
#lc2$sigma <- lc2$sigma + 1

simR <- simGLVM(n = n, p = p, form = form, dist = fam, loadmt = l1, coefs = lc)
Y <- simR$Y
Z <- simR$Z
borg <- simR$borg

#profvis({
ex1 <- GLVM.fit(Y = Y, fam = fam, form = form , silent = F, ghp = 50, iter.lim = 700, tol = 1e-7, loadmt = l1, icoefs = lc, useoptim = F, skipEM = F)
#ex2 <- GLVM.fit(Y = Y, fam = fam, form = form1, silent = F, ghp = 50, iter.lim = 700, tol = 1e-7, loadmt = l11, icoefs = lc1, useoptim = F, skipEM = F)
#ex2 <- GLVM.fit(Y = Y, fam = fam2, form = form1, silent = F, ghp = 50, iter.lim = 700, tol = 1e-7, loadmt = l11, icoefs = lc1, useoptim = F)
#ex3 <- GLVM.fit(Y = Y, fam = fam, form = form2, silent = F, ghp = 50, iter.lim = 700, tol = 1e-7, loadmt = l12, icoefs = lc2, useoptim = F)
#})

# sourceCpp("C:/Users/carde/Dropbox/Camilo and Irini/Research/GitHub/SPLVM/C/ddFun.cpp")
# ex2 <- NRb(as.matrix(Y),form,fam,1,50,maxit = 700)

round(ex1$b$mu - borg$mu,4)
round(ex1$b$sigma - borg$sigma,4)


plotGLVM(item = 1, mod = ex1, morg = simR, plot.org = F,
         plot.mean = F, plot.sd = F, quant = c(0.025,0.25,0.75,0.975),
         sep.plots = F, plot.3D = T, plot.dist = T, plot.addpoints = T)

plot.score(ex1)

ex1$b$mu; borg$mu
ex1$b$sigma; borg$sigma
#plot(ex1$cvgRes[,sample(1:ncol(ex1$cvgRes),2,replace = F)])

testB1 <- rep(1,length(unlist(ex1$b)))
testB1 <- unname(unlist(ex1$b))

num.score <- grad(func = loglikFun, x = testB1, method = "Richardson",
                  method.args=list(r = 10), Y = ex1$Y, ghQ = ex1$gr, fam = ex1$fam, beta = ex1$b)

ana.score <- ScoreFun(testB1, ex1$Y, ex1$b, ex1$gr, ex1$fam)

plot.ts(ana.score)
lines(num.score[1:length(ana.score)], col = "red")

# Think of a) change loglikFun to have as input the parameters that we actually use. or 
#          b) Change ScoreFun to return scores for non existent parameters (e.g. sigma in Poisson RV)

