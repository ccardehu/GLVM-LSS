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
source("graphFun.R")
source("EMFun.R")

n = 1000     # Number of individuals
p = 10       # Number of items
nsim = 1000  # Number of simulations
form <- list("mu" = "~ Z1")# + I(Z1^2)", "sigma" = "~ Z1")
#form1 <- list("mu" = "~ Z1 + I(Z1^2)", "sigma" = "~ 1") #  , "sigma" = "~ Z1"
#form2 <- list("mu" = "~ Z1", "sigma" = "~ 1") #  , "sigma" = "~ Z1"
fam <- rep("binomial",p)
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
#lc$mu[,3] <- runif(p,-0.7,-0.2)
#lc$mu[,c(4)] <- runif(p,-0.5,-0.1)
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
ex1 <- GLVM.fit(Y = Y, fam = fam, form = form , silent = F, ghp = 10, iter.lim = 1000, tol = sqrt(.Machine$double.eps), loadmt = l1, icoefs = lc, useoptim = F, skipEM = F)
ex2 <- GLVM.fit2(Y = Y, fam = fam, form = form , silent = F, ghp = 10, iter.lim = 1000, tol = sqrt(.Machine$double.eps), loadmt = l1, icoefs = lc, useoptim = F, skipEM = F)
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
round(ex2$b$mu - borg$mu,4)
round(ex2$b$sigma - borg$sigma,4)


plotGLVM(item = 6, mod = ex1, morg = simR, plot.org = F,
         plot.mean = T, plot.sd = T, quant = c(0.025,0.25,0.75,0.975),
         sep.plots = F, plot.3D = T, plot.dist = T, plot.addpoints = T)

# plot.score(ex1)
#pairs(Y)
#round(Matrix::tril(var(Y)),2)

cbind(ex1$b$mu, ex2$b$mu, borg$mu)
cbind(ex1$b$sigma, ex2$b$sigma, borg$sigma)
#plot(ex1$cvgRes[,sample(1:ncol(ex1$cvgRes),2,replace = F)])


modt <- ex2

xlabn <- NULL
for(i in 1:ncol(Y)){
  xlabn <- append(xlabn,paste0("mu[",i,",",1:ncol(modt$b$mu),"]"))
  if("sigma" %in% pFun(modt$fam[i])) xlabn <- append(xlabn, paste0("sigma[",i,",",1:ncol(modt$b$sigma),"]"))
  if("tau" %in% pFun(modt$fam[i])) xlabn <- append(xlabn, paste0("tau[",i,",",1:ncol(modt$b$tau),"]"))
  if("nu" %in% pFun(modt$fam[i])) xlabn <- append(xlabn, paste0("nu[",i,",",1:ncol(modt$b$nu),"]"))
}

# testB1 <- unname(unlist(modt$b))
testH1 <- decoefmod(modt$b)
# testB2 <- unname(unlist(borg))
testH2 <- decoefmod(borg)
testH3 <- unname(rep(1,length(unlist(borg))))
#testB4 <- unname(rep(15,length(unlist(borg))))

num.score1 <- grad(func = loglikFun, x = testH1, method = "Richardson", method.args=list(r = 6, v = 2),
                   Y = modt$Y, ghQ = modt$gr, fam = modt$fam, beta = modt$b)
num.score2 <- grad(func = loglikFun, x = testH2, method = "Richardson", method.args=list(r = 6, v = 2),
                   Y = modt$Y, ghQ = modt$gr, fam = modt$fam, beta = modt$b)
num.score3 <- grad(func = loglikFun, x = testH3, method = "Richardson", method.args=list(r = 6, v = 2),
                   Y = modt$Y, ghQ = modt$gr, fam = modt$fam, beta = modt$b)
# num.score4 <- grad(func = loglikFun, x = testB4, method = "Richardson",
                   # method.args=list(r = 10), Y = modt$Y, ghQ = modt$gr, fam = modt$fam, beta = modt$b)

num.hess1 <- hessian(func = loglikFun, x = testH1, method = "Richardson", method.args=list(r = 6, v = 2),
                     Y = modt$Y, ghQ = modt$gr, fam = modt$fam, beta = modt$b);
num.hess2 <- hessian(func = loglikFun, x = testH2, method = "Richardson", method.args=list(r = 6, v = 2),
                     Y = modt$Y, ghQ = modt$gr, fam = modt$fam, beta = modt$b);
num.hess3 <- hessian(func = loglikFun, x = testH3, method = "Richardson", method.args=list(r = 6, v = 2),
                     Y = modt$Y, ghQ = modt$gr, fam = modt$fam, beta = modt$b);
#num.hess1 <- diag(num.hess1);  num.hess2 <- diag(num.hess2);  num.hess3 <- diag(num.hess3)

plot.jac = F
if(plot.jac == T){
num.jac1 <- jacobian(func = ScoreFun, x = testH1, method = "Richardson", method.args=list(r = 6, v = 2),
                     Y = modt$Y, ghQ = modt$gr, fam = modt$fam, beta = modt$b)
num.jac2 <- jacobian(func = ScoreFun, x = testH2, method = "Richardson", method.args=list(r = 6, v = 2),
                     Y = modt$Y, ghQ = modt$gr, fam = modt$fam, beta = modt$b)
num.jac3 <- jacobian(func = ScoreFun, x = testH3, method = "Richardson", method.args=list(r = 6, v = 2),
                     Y = modt$Y, ghQ = modt$gr, fam = modt$fam, beta = modt$b)
}

ana1 <- ScHesFun(testH1, modt$Y, modt$b, modt$gr, modt$fam)
ana2 <- ScHesFun(testH2, modt$Y, modt$b, modt$gr, modt$fam)
ana3 <- ScHesFun(testH3, modt$Y, modt$b, modt$gr, modt$fam)

profvis::profvis({t1 <- ScHesFun(testH1, modt$Y, modt$b, modt$gr, modt$fam)})

ana.score1 <- ana1$score
ana.score2 <- ana2$score
ana.score3 <- ana3$score
# ana.score4 <- ScoreFun(testB4, modt$Y, modt$b, modt$gr, modt$fam)

ana.hess1 <- ana1$hess
ana.hess2 <- ana2$hess
ana.hess3 <- ana3$hess


par("mar" = c(6, 4, 4, 2) + 0.1)
plot(ana.score1, main = "Score (first deriv. of log-likelihood) @ MLE", xlab = "", ylab = "Value", pch = 16,
     col = "gray50", xaxt = "n")
mtext(side = 1, text = "Parameter index", line = 4)
for(i in seq_along(xlabn)){ axis(1, at=i, labels = xlabn[i] , las = 2, cex.axis = 0.7) }
abline(h = 0, col = "forestgreen", lty = 2)
abline(v = seq(from = ((sum(length(modt$b$mu), length(modt$b$sigma)) %/% ncol(Y))),length.out = ncol(Y)-1,
               by = ((sum(length(modt$b$mu), length(modt$b$sigma)) %/% ncol(Y)))) + 0.5, lty = 3, col = "gray80")
points(num.score1[1:length(ana.score1)], col = "red", pch = 10)
segments(x0 = seq_along(ana.score1), y0 = unlist(ana.score1), y1 = num.score1, col = "blue", lty = 3, lwd = 1)
legend("topleft",legend=c("Analytical", "Numerical"), col=c("gray50", "red"),
       pch = c(16,10), cex=1, inset = 0.02, box.col = "white")

# plot(ana.score1 - num.score1,main = "Differences in Score values: Analytical vs. Numerical (@ MLE)",xlab ="", ylab = "Value", pch = 16, col = "gray50", xaxt = "n")
# mtext(side = 1, text = "Parameter index", line = 4)
# for(i in seq_along(xlabn)){ axis(1, at=i, labels = xlabn[i] , las = 2, cex.axis = 0.7) }
# abline(h = 0, col = "forestgreen", lty = 2)
# abline(v = seq(from = ((sum(length(modt$b$mu), length(modt$b$sigma)) %/% ncol(Y))),length.out = ncol(Y)-1,
#                by = ((sum(length(modt$b$mu), length(modt$b$sigma)) %/% ncol(Y)))) + 0.5, lty = 3, col = "gray80")

plot(ana.score2, main = "Score (first deriv. of log-likelihood) @ Original Parameters", xlab = "", ylab = "value", pch = 16,
     col = "gray50", xaxt = "n")
mtext(side = 1, text = "Parameter index", line = 4)
for(i in seq_along(xlabn)){ axis(1, at=i, labels = xlabn[i] , las = 2, cex.axis = 0.7) }
abline(h = 0, col = "forestgreen", lty = 2)
abline(v = seq(from = ((sum(length(modt$b$mu), length(modt$b$sigma)) %/% ncol(Y))),length.out = ncol(Y)-1,
               by = ((sum(length(modt$b$mu), length(modt$b$sigma)) %/% ncol(Y)))) + 0.5, lty = 3, col = "gray80")
points(num.score2[1:length(ana.score2)], col = "red", pch = 10)
segments(x0 = seq_along(ana.score2), y0 = unlist(ana.score2), y1 = num.score2, col = "blue", lty = 3, lwd = 1)
legend("bottomleft",legend=c("Analytical", "Numerical"), col=c("gray50", "red"), pch = c(16,10), cex=1, inset = 0.02, box.col = "white")

# plot(ana.score2 - num.score2,main = "Differences in Score values: Analytical vs. Numerical (@ Original)",xlab ="", ylab = "Value", pch = 16, col = "gray50", xaxt = "n")
# mtext(side = 1, text = "Parameter index", line = 4)
# for(i in seq_along(xlabn)){ axis(1, at=i, labels = xlabn[i] , las = 2, cex.axis = 0.7) }
# abline(h = 0, col = "forestgreen", lty = 2)
# abline(v = seq(from = ((sum(length(modt$b$mu), length(modt$b$sigma)) %/% ncol(Y))),length.out = ncol(Y)-1,
#                by = ((sum(length(modt$b$mu), length(modt$b$sigma)) %/% ncol(Y)))) + 0.5, lty = 3, col = "gray80")

plot(ana.score3, main = "Score (first deriv. of log-likelihood) @ Vector of 1s", xlab = "", ylab = "value", pch = 16,
     col = "gray50", xaxt = "n")
mtext(side = 1, text = "Parameter index", line = 4)
for(i in seq_along(xlabn)){ axis(1, at=i, labels = xlabn[i] , las = 2, cex.axis = 0.7) }
abline(h = 0, col = "forestgreen", lty = 2)
abline(v = seq(from = ((sum(length(modt$b$mu), length(modt$b$sigma)) %/% ncol(Y))),length.out = ncol(Y)-1,
               by = ((sum(length(modt$b$mu), length(modt$b$sigma)) %/% ncol(Y)))) + 0.5, lty = 3, col = "gray80")
points(num.score3[1:length(ana.score3)], col = "red", pch = 10)
segments(x0 = seq_along(ana.score3), y0 = unlist(ana.score3), y1 = num.score3, col = "blue", lty = 3, lwd = 1)
legend("topleft",legend=c("Analytical", "Numerical"), col=c("gray50", "red"), pch = c(16,10), cex=1, inset = 0.02, box.col = "white")

# plot(ana.score3 - num.score3,main = "Differences in Score values: Analytical vs. Numerical (@ 1s)",xlab ="", ylab = "Value", pch = 16, col = "gray50", xaxt = "n")
# mtext(side = 1, text = "Parameter index", line = 4)
# for(i in seq_along(xlabn)){ axis(1, at=i, labels = xlabn[i] , las = 2, cex.axis = 0.7) }
# abline(h = 0, col = "forestgreen", lty = 2)
# abline(v = seq(from = ((sum(length(modt$b$mu), length(modt$b$sigma)) %/% ncol(Y))),length.out = ncol(Y)-1,
#                by = ((sum(length(modt$b$mu), length(modt$b$sigma)) %/% ncol(Y)))) + 0.5, lty = 3, col = "gray80")

int <- 1:sum(length(modt$b$mu), length(modt$b$sigma))^2
# int <- 1:((sum(length(modt$b$mu), length(modt$b$sigma))^2 %/% ncol(Y)) * 2)
plot(c(ana.hess1)[int], main = "Hessian (2nd deriv. of log-likelihood) @ MLE",xlab ="Hessian entry Index", ylab = "Value",
     cex = 0.5, pch = 16, col = "gray50", ylim = c(min(ana.hess1,num.hess1),max(ana.hess1,num.hess1)))#, xaxt = "n")
abline(h = 0, col = "forestgreen", lty = 2);
abline(v = seq(from = sum(length(modt$b$mu), length(modt$b$sigma))^2 %/% ncol(Y),length.out = ncol(Y)-1,
               by = sum(length(modt$b$mu), length(modt$b$sigma))^2 %/% ncol(Y)) + 0.5, lty = 3, col = "gray80")
points(c(num.hess1)[int], col = "red", pch = 16, cex = 0.5)
segments(x0 = seq_along(ana.hess1[int]), y0 = unlist(ana.hess1)[int], y1 = num.hess1[int], col = "blue", lty = 3, lwd = 1)
legend(0,-2500,legend=c("Analytical", "Numerical (Hessian)"), col=c("gray50", "red",0.6),
       pch = c(16,16), cex=0.8, inset = 0.02, box.col = "white")
if(plot.jac == T){ points(c(num.jac1)[int], col = scales::alpha("blue",0.6), pch = 16, cex = 0.5) 
  legend("bottomright",legend=c("Analytical", "Numerical (Hessian)", "Numerical (Jacobian)"),
         col=c("gray50",scales::alpha("red",0.6), scales::alpha("blue",0.6)), pch = c(16,16,16),
         cex=1, inset = 0.02, box.col = "white")}

# plot(c(num.hess1 - ana.hess1), main = "Differences in Hessian Analytical vs. Numerical (@ MLE)",xlab ="Hessian entry Index", ylab = "Value", pch = 16, col = "gray50")
# abline(h = 0, col = "forestgreen", lty = 2)
# abline(v = seq(from = sum(length(modt$b$mu), length(modt$b$sigma))^2 %/% ncol(Y),length.out = ncol(Y)-1,
#               by = sum(length(modt$b$mu), length(modt$b$sigma))^2 %/% ncol(Y)) + 0.5, lty = 3, col = "gray80")

plot(c(ana.hess2)[int], main = "Hessian (2nd deriv. of log-likelihood) @ Original Parameters",
     xlab = "Hessian entry Index", ylab = "value", cex = 0.5, pch = 16, col = "gray50")
abline(h = 0, col = "forestgreen", lty = 2);
abline(v = seq(from = sum(length(modt$b$mu), length(modt$b$sigma))^2 %/% ncol(Y),length.out = ncol(Y)-1,
               by = sum(length(modt$b$mu), length(modt$b$sigma))^2 %/% ncol(Y)) + 0.5, lty = 3, col = "gray80")
points(c(num.hess2)[int], col = "red", pch = 16, cex = 0.5)
segments(x0 = seq_along(ana.hess2[int]), y0 = unlist(ana.hess2)[int], y1 = num.hess2[int], col = "blue", lty = 3, lwd = 1)
legend(0,-3000,legend=c("Analytical", "Numerical (Hessian)"), col=c("gray50", "red"),
       pch = c(16,16), cex=1, inset = 0.02, box.col = "white")
if(plot.jac == T){ points(c(num.jac2)[int], col = scales::alpha("blue",0.6), pch = 16, cex = 0.5) 
  legend("bottomright",legend=c("Analytical", "Numerical (Hessian)", "Numerical (Jacobian)"),
         col=c("gray50",scales::alpha("red",0.6), scales::alpha("blue",0.6)), pch = c(16,16,16),
         cex=1, inset = 0.02, box.col = "white")}

# plot(c(num.hess2 - ana.hess2),main = "Differences in Hessian Analytical vs. Numerical (@ Original)",xlab ="Hessian entry Index", ylab = "Value", pch = 16, col = "gray50")
# abline(h = 0, col = "forestgreen", lty = 2)
# abline(v = seq(from = sum(length(modt$b$mu), length(modt$b$sigma))^2 %/% ncol(Y),length.out = ncol(Y)-1,
#               by = sum(length(modt$b$mu), length(modt$b$sigma))^2 %/% ncol(Y)) + 0.5, lty = 3, col = "gray80")

plot(c(ana.hess3)[int], main = "Hessian (2nd deriv. of log-likelihood) @ Vector of 1s",
     xlab = "Hessian entry Index", ylab = "value", cex = 0.5, pch = 16, col = "gray50", ylim = c(min(ana.hess3,num.hess3),max(ana.hess3,num.hess3)))
abline(h = 0, col = "forestgreen", lty = 2);
abline(v = seq(from = sum(length(modt$b$mu), length(modt$b$sigma))^2 %/% ncol(Y),length.out = ncol(Y)-1,
               by = sum(length(modt$b$mu), length(modt$b$sigma))^2 %/% ncol(Y)) + 0.5, lty = 3, col = "gray80")
points(c(num.hess3)[int], col = "red", pch = 16, cex = 0.5)
segments(x0 = seq_along(ana.hess3[int]), y0 = unlist(ana.hess3)[int], y1 = num.hess3[int], col = "blue", lty = 3, lwd = 1)
legend(2000,-850,legend=c("Analytical", "Numerical (Hessian)"), col=c("gray50", "red",0.6),
       pch = c(16,16), cex=0.8, inset = 0.02, box.col = "white")
if(plot.jac == T){ points(c(num.jac3)[int], col = scales::alpha("blue",0.6), pch = 16, cex = 0.5) 
  legend("bottomleft",legend=c("Analytical", "Numerical (Hessian)", "Numerical (Jacobian)"),
         col=c("gray50",scales::alpha("red",0.6), scales::alpha("blue",0.6)), pch = c(16,16,16),
         cex=1, inset = 0.02, box.col = "white")}

# plot(c(num.hess3 - ana.hess3),main = "Differences in Analytical Hessian vs. Numerical Hessian (@ 1s)",xlab ="Hessian entry Index", ylab = "Value", pch = 16, col = "gray50")
# abline(h = 0, col = "forestgreen", lty = 2)
# abline(v = seq(from = sum(length(modt$b$mu), length(modt$b$sigma))^2 %/% ncol(Y),length.out = ncol(Y)-1,
#               by = sum(length(modt$b$mu), length(modt$b$sigma))^2 %/% ncol(Y)) + 0.5, lty = 3, col = "gray80")

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
