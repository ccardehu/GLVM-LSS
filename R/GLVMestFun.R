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
source("GLVMfitv2.R")
source("GLVMfitv3.R")
source("graphFun.R")
source("EMFun.R")

n = 100     # Number of individuals
p = 10       # Number of items
nsim = 1000  # Number of simulations
form <- list("mu" = "~ Z1", "sigma" = "~ 1")
#form1 <- list("mu" = "~ Z1 + I(Z1^2)", "sigma" = "~ 1") #  , "sigma" = "~ Z1"
#form2 <- list("mu" = "~ Z1", "sigma" = "~ 1") #  , "sigma" = "~ Z1"
fam <- rep("normal",p)
# fam <- c(rep("normal", p/2), rep("ZIpoisson", p/2)) # rep("poisson",p)#
# sample(c("normal","poisson","binom"),size = 10,replace = T)

l1 <- NULL
l1$mu <- matrix(1,ncol = 2, nrow = p)
l1$sigma <- matrix(1, ncol = 1, nrow = p)

# Restrictions
# ____________
#
#l1$mu[1,c(4)] <- 0
#l1$mu[1:p/2,c(3)] <- 0
#l1$mu[(p/2+1):p,c(2)] <- 0
#l1$mu[2,c(2,4)] <- 0
#l1$mu[3,c(2,3)] <- 0
#l1$sigma[c(2,4,9),2] <- 0

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
lc$mu[,1] <- runif(p,1,2)
lc$mu[,c(2)] <- runif(p,1.0,3.0)
# lc$mu[,3] <- runif(p,-0.4,-0.2)
# lc$mu[,4] <- runif(p,-0.7,-0.2)
lc$sigma <- matrix(runif(length(l1$sigma), min = 0.1, max = 0.3), nrow = p)
#lc$sigma[,1] <- runif(p,-1,1)
#lc$sigma[,2] <- runif(p,2,4)
# lc$sigma[,3] <- runif(p,-0.7,0.1)
#lc$sigma[,4] <- runif(p,-0.5,-0.1)
# 
# lc. <- NULL
# lc.$mu <- l1$mu * 0# + runif(length(l1$mu),-1,1)
# lc.$sigma <- l1$sigma * 0# + runif(length(l1$sigma),-1,1)

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

# profvis::profvis({
# rbenchmark::benchmark(
# "ex1" = {
ex1 <- GLVM.fit(Y = Y, fam = fam, form = form , silent = F, ghp = 25, iter.lim = 1000, tol = sqrt(.Machine$double.eps), loadmt = l1, icoefs = lc, useoptim = F, skipEM = F)
# },
# "ex2" = {
# exint <- GLVM.fit3(Y = Y, fam = fam, form = form , silent = T, ghp = 10, iter.lim = 10, tol = sqrt(.Machine$double.eps), loadmt = l1, icoefs = lc, useoptim = F, skipEM = F)
# ex2 <- GLVM.fit2(Y = Y, fam = fam, form = form , silent = T, ghp = 10, iter.lim = 1000, tol = sqrt(.Machine$double.eps), loadmt = l1, icoefs = lc, useoptim = F, skipEM = F)
# },
# "ex3" = {
ex2 <- GLVM.fit3(Y = Y, fam = fam, form = form , silent = F, ghp = 25, iter.lim = 1000, tol = sqrt(.Machine$double.eps), loadmt = l1, icoefs = lc, useoptim = F, skipEM = F)
# },
# replications = 1, columns = c("test", "replications", "elapsed","relative", "user.self", "sys.self") )
# })
# 
rmtest <- ltm::ltm(Y ~ z1,IRT.param = F, control = list(GHk = 25, iter.em = 30))
rmtestfac <- factor.scores(rmtest, method = "EAP")
round(cbind(coef(rmtest), ex1$b$mu, borg$mu),5)

# rmtest <- ltm::ltm(Y ~ z1*z2,IRT.param = F,constraint = matrix(c(1,1,2,2,3,4,2,4,0,0,0,0),nrow = 4),control = list(GHk = 30, iter.em = 50))
# round(cbind(coef(rmtest), ex1$b$mu, ex2$b$mu,borg$mu),2)

# library(lavaan)
# CFA.model <- ' Z1 =~ 0*Y1+Y2+Y3+Y4+Y5+Y6+Y7+Y8+Y9+Y10 
#                Z2 =~ NA*Y1+0*Y2+Y3+Y4+Y5+Y6+Y7+Y8+Y9+Y10 '
# rmtest <- cfa(CFA.model, data = Y, orthogonal = T, meanstructure = TRUE, std.lv = TRUE)
# coef(rmtest)
# summary(rmtest)
# rmtest <- psych::fa(Y,2, fm="mle", rotate = "none")$loadings

# sourceCpp("C:/Users/carde/Dropbox/Camilo and Irini/Research/GitHub/SPLVM/C/ddFun.cpp")
# ex2 <- NRb(as.matrix(Y),form,fam,1,50,maxit = 700)

round(ex1$b$mu - borg$mu,4)
round(ex1$b$sigma - borg$sigma,4)
round(ex2$b$mu - borg$mu,4)
round(ex2$b$sigma - borg$sigma,4)


plotGLVM(item = 10, mod = ex1, morg = simR, plot.org = T,
         plot.mean = T, plot.sd = T,# quant = c(0.025,0.25,0.75,0.975),
         sep.plots = F, plot.3D = T, plot.dist = T, plot.addpoints = T)

# plot.score(ex1)
#pairs(Y)
#round(Matrix::tril(var(Y)),2)

round(cbind(ex1$b$mu, ex2$b$mu, borg$mu),3)
round(cbind(ex1$b$sigma, ex2$b$sigma, borg$sigma),3)
#plot(ex1$cvgRes[,sample(1:ncol(ex1$cvgRes),2,replace = F)])



# Think of a) Change loglikFun to have as input the parameters that we actually use. or 
#          b) Change ScoreFun to return scores for non existent parameters (e.g. sigma in Poisson RV)


hist(Z$mu$Z1, breaks = 100, freq = F, xlim = c(-5,5), border = "gray", ylim = c(0,0.7))
# lines(density(rmtestfac$score.dat$z1), col = 5, lwd = 2)
lines(density(Z$mu$Z1), col = 1, lwd = 2)
lines(density(zsc(ex1)$Z1), col = 2, lwd = 2)
# lines(density((rowSums(Y)-mean(rowSums(Y)))/sd(rowSums(Y))), col = 3, lwd = 2)
# hist(rmtestfac$score.dat$z1, breaks = 100, freq = F, xlim = c(-5,5), border = "gray")
lines(density(rep(rmtestfac$score.dat$z1,rmtestfac$score.dat$Obs)), col = "orange", lwd = 2)
hist(zsc(ex1)$Z1, breaks = 100, freq = F, xlim = c(-5,5), border = "gray")
lines(density(zsc(ex1)$Z1), col = 2, lwd = 2)
plot(Z$mu$Z1, zsc(ex1)$Z1)
# cor(Z$mu$Z1, rmtestfac$score.dat$z1,method = "kendall")
cor(Z$mu$Z1, zsc(ex1)$Z1,method = "kendall")
cor(Z$mu$Z1, rowSums(Y),method = "kendall")
cor(Z$mu$Z1, zsc(ex1)$Z1,method = "pearson")
cor(Z$mu$Z1, rowSums(Y),method = "pearson")
# plot(Z$mu$Z2, zsc(ex1)$Z2)
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


library(mvabund)
data <- cbind(antTraits$abund)
data <- as.matrix(unname(data))
p. <- ncol(data)
famd <- rep("poisson",p.)
l1. <- NULL; l1.$mu <- matrix(1,ncol = 2,nrow = p.)
formd <- NULL; formd$mu <- "~ Z1"
exPois <- GLVM.fit(Y = data, fam = famd, form = formd , silent = F, ghp = 40, iter.lim = 800, tol = sqrt(.Machine$double.eps), loadmt = l1., useoptim = F, skipEM = F)
faPois <- zsc(exPois)
plot(faPois$Z1,data[,41])

fitp <- gllvm::gllvm(data,family = poisson(),num.lv = 1)
summary(fitp)$Coefficients

