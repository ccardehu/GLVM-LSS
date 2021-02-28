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
p = 10       # Number of items
nsim = 1000  # Number of simulations
form <- list("mu" = "~ Z1", "sigma" = "~ Z1")
#form1 <- list("mu" = "~ Z1 + I(Z1^2)", "sigma" = "~ 1") #  , "sigma" = "~ Z1"
#form2 <- list("mu" = "~ Z1", "sigma" = "~ 1") #  , "sigma" = "~ Z1"
fam <- rep("normal",p)
# fam <- c(rep("normal", p/2), rep("poisson", p/2)) # rep("poisson",p)#
# sample(c("normal","poisson","binom"),size = 10,replace = T)

l1 <- NULL
l1$mu <- matrix(1,ncol = 2, nrow = p)
l1$sigma <- matrix(1, ncol = 2, nrow = p)

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
lc$mu <- matrix(runif(length(l1$mu), min = -0.5, max = 0.5),nrow = p)
#lc$mu[,1] <- runif(p,1,2)
#lc$mu[,2] <- runif(p,0.5,1.5)
#lc$mu[,3] <- runif(p,-0.5,0.5)
#lc$mu[,c(4,7)] <- runif(p,-0.5,-0.1)
lc$sigma <- matrix(runif(length(l1$sigma), min = 0.1, max = 0.3), nrow = p)
#lc$sigma[,1] <- runif(p,-1,1)
#lc$sigma[,2] <- runif(p,2,4)
#lc$sigma[,3] <- runif(p,-0.5,-0.1)
#lc$sigma[,4] <- runif(p,-0.5,-0.1)
# 
# lc. <- NULL
# lc.$mu <- l1$mu * -4
# lc.$sigma <- l1$sigma * -2

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
ex1 <- GLVM.fit(Y = Y, fam = fam, form = form , silent = F, ghp = 30, iter.lim = 1000, tol = .Machine$double.eps, loadmt = l1, icoefs = lc, useoptim = F, skipEM = F)
#ex10 <- GLVM.fit(Y = Y, fam = fam, form = form , silent = F, ghp = 30, iter.lim = 1000, tol = .Machine$double.eps, loadmt = l1, icoefs = lc, useoptim = T, skipEM = F)
#ex11 <- GLVM.fit(Y = Y, fam = fam, form = form , silent = F, ghp = 10, iter.lim = 700, tol = 1e-7, loadmt = l1, icoefs = lc, useoptim = F, skipEM = F)
#ex10 <- GLVM.fit(Y = Y, fam = fam, form = form , silent = F, ghp = 50, iter.lim = 700, tol = 1e-7, loadmt = l1, icoefs = lc, useoptim = F, skipEM = F)
#ex2 <- GLVM.fit(Y = Y, fam = fam, form = form1, silent = F, ghp = 50, iter.lim = 700, tol = 1e-7, loadmt = l11, icoefs = lc1, useoptim = F, skipEM = F)
#ex2 <- GLVM.fit(Y = Y, fam = fam2, form = form1, silent = F, ghp = 50, iter.lim = 700, tol = 1e-7, loadmt = l11, icoefs = lc1, useoptim = F)
#ex3 <- GLVM.fit(Y = Y, fam = fam, form = form2, silent = F, ghp = 50, iter.lim = 700, tol = 1e-7, loadmt = l12, icoefs = lc2, useoptim = F)
#})

# sourceCpp("C:/Users/carde/Dropbox/Camilo and Irini/Research/GitHub/SPLVM/C/ddFun.cpp")
# ex2 <- NRb(as.matrix(Y),form,fam,1,50,maxit = 700)

round(ex1$b$mu - borg$mu,4)
round(ex1$b$sigma - borg$sigma,4)


plotGLVM(item = 5, mod = ex1, morg = simR, plot.org = F,
         plot.mean = T, plot.sd = T, quant = c(0.025,0.25,0.75,0.975),
         sep.plots = F, plot.3D = T, plot.dist = T, plot.addpoints = T)

# plot.score(ex1)
pairs(Y)
round(Matrix::tril(cor(Y)),2)

ex1$b$mu; borg$mu
ex1$b$sigma; borg$sigma
#plot(ex1$cvgRes[,sample(1:ncol(ex1$cvgRes),2,replace = F)])


modt <- ex1

testB1 <- unname(unlist(modt$b))
testB2 <- unname(unlist(borg))
testB3 <- unname(rep(1,length(unlist(borg))))
#testB4 <- unname(rep(15,length(unlist(borg))))

num.score1 <- grad(func = loglikFun, x = testB1, method = "Richardson",
                   method.args=list(r = 10, v = 10), Y = modt$Y, ghQ = modt$gr, fam = modt$fam, beta = modt$b)
num.score2 <- grad(func = loglikFun, x = testB2, method = "Richardson",
                   method.args=list(r = 10, v = 10), Y = modt$Y, ghQ = modt$gr, fam = modt$fam, beta = modt$b)
num.score3 <- grad(func = loglikFun, x = testB3, method = "Richardson",
                   method.args=list(r = 10, v = 10), Y = modt$Y, ghQ = modt$gr, fam = modt$fam, beta = modt$b)
# num.score4 <- grad(func = loglikFun, x = testB4, method = "Richardson",
                   # method.args=list(r = 10), Y = modt$Y, ghQ = modt$gr, fam = modt$fam, beta = modt$b)

ana.score1 <- ScoreFun(testB1, modt$Y, modt$b, modt$gr, modt$fam)
ana.score2 <- ScoreFun(testB2, modt$Y, modt$b, modt$gr, modt$fam)
ana.score3 <- ScoreFun(testB3, modt$Y, modt$b, modt$gr, modt$fam)
# ana.score4 <- ScoreFun(testB4, modt$Y, modt$b, modt$gr, modt$fam)

plot(ana.score1, main = "Score (first deriv. of log-likelihood) @ MLE", xlab = "Parameter index", ylab = "value", pch = 16, col = "gray50", ylim = c(1.5*min(ana.score1), 1.5*max(ana.score1)))
abline(h = 0, col = "forestgreen", lty = 2); abline(v = 20.5, lty = 2); abline(v = 10.5, lty = 2); abline(v = 30.5, lty = 2)
points(num.score1[1:length(ana.score1)], col = "red", pch = 10)
segments(x0 = seq_along(ana.score1), y0 = unlist(ana.score1), y1 = num.score1, col = "blue", lty = 3, lwd = 1)
legend("bottomright",legend=c("Analytical", "Numerical"), col=c("gray50", "red"), pch = c(16,10), cex=1, inset = 0.02, bty = "n")
#legend(19,-0.0025,legend=c("Sigma Params."), cex=1, inset = 0.02, bty = "n")
#points(unlist(ex1$Score),col = "green", pch = 16)

plot(ana.score2, main = "Score (first deriv. of log-likelihood) @ Original Parameters", xlab = "Parameter index", ylab = "value", pch = 16, col = "gray50")
abline(h = 0, col = "forestgreen", lty = 2); abline(v = 20.5, lty = 2); abline(v = 10.5, lty = 2); abline(v = 30.5, lty = 2)
points(num.score2[1:length(ana.score2)], col = "red", pch = 10)
segments(x0 = seq_along(ana.score2), y0 = unlist(ana.score2), y1 = num.score2, col = "blue", lty = 3, lwd = 1)
legend("bottomright",legend=c("Analytical", "Numerical"), col=c("gray50", "red"), pch = c(16,10), cex=1, inset = 0.02, bty = "n")
#legend(19,-20,legend=c("Sigma Params."), cex=1, inset = 0.02, bty = "n")

plot(ana.score3, main = "Score (first deriv. of log-likelihood) @ Vector of 1s", xlab = "Parameter index", ylab = "value", pch = 16, col = "gray50")
abline(h = 0, col = "forestgreen", lty = 2); abline(v = 20.5, lty = 2); abline(v = 10.5, lty = 2); abline(v = 30.5, lty = 2)
points(num.score3[1:length(ana.score3)], col = "red", pch = 10)
segments(x0 = seq_along(ana.score3), y0 = unlist(ana.score3), y1 = num.score3, col = "blue", lty = 3, lwd = 1)
legend("bottomleft",legend=c("Analytical", "Numerical"), col=c("gray50", "red"), pch = c(16,10), cex=1, inset = 0.02, bty = "n")
#legend(19,-300,legend=c("Sigma Params."), cex=1, inset = 0.02, bty = "n")


# Think of a) Change loglikFun to have as input the parameters that we actually use. or 
#          b) Change ScoreFun to return scores for non existent parameters (e.g. sigma in Poisson RV)


hist(Z$mu$Z1, breaks = 100, freq = F, xlim = c(-5,5))
points(ex1$gr$points, ex1$gr$weights, col = "red", pch = 16)
lines(seq(-5,5,length.out = 100), dnorm(seq(-5,5,length.out = 100)), col = "blue", lwd = 2)

data <- (psych::bfi)
data <- data[complete.cases(data),1:5]
p. <- ncol(data)
famd <- rep("normal", p.)
l1. <- l2. <- NULL
l1.$mu <- matrix(1,ncol = 3, nrow = p.)
l2.$mu <- matrix(0.5,ncol = 3, nrow = p.)
l1.$sigma <- matrix(1, ncol = 1, nrow = p.)
l2.$sigma <- matrix(0.5, ncol = 1, nrow = p.)
formd <- NULL; formd$mu <- "~ Z1 + I(Z1^2)"; formd$sigma <- "~ 1"
exBin <- GLVM.fit(Y = data, fam = famd, form = formd , silent = F, ghp = 25, iter.lim = 800, tol = 1e-7, loadmt = l1., useoptim = F, skipEM = F, icoefs = l2.)
exBin$b
exBin$loglik

summary(ltm(WIRS ~ z1*z2, constraint = matrix(c(1:6,rep(4,6),rep("A1",3), rep(1,3)),ncol = 3)))

ltm(Y ~ z1, IRT.param = F, control = list(GHk = 30))
round(ex1$b$mu,3)
ex1$loglik
