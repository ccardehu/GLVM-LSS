
rm(list = ls())
set.seed(1234)
# rm(.Random.seed, envir=globalenv())

library(MASS)
library(abind)
library(splines)
library(mvtnorm)
library(numDeriv)

# `%!in%` <- Negate(`%in%`)

source("fa_main.R")
source("fb_auxf.R")
source("fc_mfit.R")
source("fd_mghq.R")
source("ff_msim.R")

n = 1000     # Number of individuals
p = 10       # Number of items
nsim = 1000  # Number of simulations
s.form <- list("mu" = "~ Z1 + I(Z1^2)")#, "sigma" = "~ 1")
fam <- rep("binomial",p)

l1 <- NULL
l1$mu <- matrix(1,ncol = 3, nrow = p)
# l1$sigma <- matrix(1, ncol = 1, nrow = p)

l1. <- l1
l1.$mu[sample(1:p, 5, replace = F),3] <- 0

lc <- NULL
lc$mu <- matrix(rnorm(length(l1$mu)),nrow = p)
lc$mu[,1] <- runif(p,-0.5,0.5)
lc$mu[,2] <- runif(p,1,2)
# lc$mu[,3] <- runif(p,1.0,2.0)
# lc$mu[,4] <- runif(p,-2.0,-1.0)
# lc$sigma <- matrix(runif(length(l1$sigma), min = 0.1, max = 0.3), nrow = p)

simR <- splvm.sim(n = n, form = s.form, fam = fam, constraints = l1, coefs = lc)
Y <- simR$Y
Z <- simR$Z
borg <- simR$b
e.form <- s.form

# lc. <- NULL; a <- 1
# lc.$mu <- lc$mu + runif(length(lc$mu),-a,a)# *0 + l1$mu
# lc.$sigma <- lc$sigma + runif(length(lc$sigma),-a,a); rm(a)

# l1. <- l1
# l1.$mu <- l1$mu*1
# l1.$sigma <- l1$sigma*1

# profvis::profvis({
testa1 <- splvm.fit(Y,fam,e.form,
          control = list(method = "EM", constraint = l1, full.hess = F, start.val = lc,
          ghQqp = 20, iter.lim = 200, tol = sqrt(.Machine$double.eps), silent = F))

testa2 <- splvm.fit(Y,fam,e.form,
          control = list(method = "PEM", constraint = l1, full.hess = F, start.val = lc,
          ghQqp = 20, iter.lim = 200, tol = sqrt(.Machine$double.eps), silent = F,
          pml.control = list(type = "lasso", lambda = 1, w.alasso = "par", gamma = 1, a = 3.7)))
# })

# library(ltm)
# rmltma <- ltm::ltm(Y ~ z1, IRT.param = F, start.val = lc$mu, control = list(GHk = 20, iter.em = 200))
# coef(rmltma)*l1$mu
# rmltma$log.Lik
# # 
# library(lavaan)
# CFA.model <- ' Z1 =~ NA*Y1+Y2+Y3+Y4+Y5+Y6+Y7+Y8+Y9+Y10 '
# rmtest <- cfa(CFA.model, data = Y, orthogonal = T, meanstructure = TRUE, std.lv = TRUE)
# coef(rmtest)
# summary(rmtest)

testb1 <- splvm.fit(Y,fam,e.form,
          control = list(method = "EM", constraint = l1, full.hess = T, start.val = lc, 
          ghQqp = 20, iter.lim = 250, tol = sqrt(.Machine$double.eps), silent = F))

testb2 <- splvm.fit(Y,fam,e.form,
          control = list(method = "PEM", constraint = l1, full.hess = T, start.val = lc, 
          ghQqp = 20, iter.lim = 250, tol = sqrt(.Machine$double.eps), silent = F,
          pml.control = list(type = "lasso", lambda = 1, w.alasso = NULL,gamma = 1, a = 3.7)))

testc1 <- splvm.fit(Y,fam,e.form,
          control = list(method = "ML", constraint = l1, start.val = lc,
          ghQqp = 20, iter.lim = 250, tol = sqrt(.Machine$double.eps), silent = F))

testc2 <- splvm.fit(Y,fam,e.form,
          control = list(method = "PML", constraint = l1, start.val = lc,
          ghQqp = 20, iter.lim = 250, tol = sqrt(.Machine$double.eps), silent = F,
          pml.control = list(type = "alasso", lambda = 1, w.alasso = 1, gamma = 1, a = 3.7)))

testd1 <- splvm.fit(Y,fam,e.form,
          control = list(method = "hybrid", constraint = l1, full.hess = F, start.val = lc,
          EM.iter.lim = 10, ghQqp = 20, iter.lim = 250, tol = sqrt(.Machine$double.eps), silent = F))

testd2 <- splvm.fit(Y,fam,e.form,
          control = list(method = "P-hybrid", constraint = l1, full.hess = F, start.val = lc,
          EM.iter.lim = 10, ghQqp = 20, iter.lim = 250, tol = sqrt(.Machine$double.eps), silent = F,
          pml.control = list(type = "lasso", lambda = 1, w.alasso = NULL,gamma = 1, a = 3.7)))

round(cbind(testa1$b$mu, testb1$b$mu, testc1$b$mu, testd1$b$mu, borg$mu),3)
round(cbind(testa2$b$mu, testb2$b$mu, testc2$b$mu, testd2$b$mu, borg$mu),3)
round(cbind(testa1$b$sigma, testb1$b$sigma, testc1$b$sigma, testd1$b$sigma, borg$sigma),3)
round(cbind(testa2$b$sigma, testb2$b$sigma, testc2$b$sigma, testd2$b$sigma, borg$sigma),3)

zsco <- fscore(testc1) # check function

b1 <- lb2cb(testa$b)
b2 <- lb2cb(lc)
b3 <- rep(1,length(b1))*lb2cb(l1)
b4 <- lb2cb(lc.)

modt <- testa
num.score <- grad(func = llkf, x = b4, method = "Richardson", method.args=list(r = 6, v = 2),
                   Y = modt$Y, ghQ = modt$ghQ, fam = modt$fam, bg = modt$b)
num.hess <- hessian(func = llkf, x = b4, method = "Richardson", method.args=list(r = 6, v = 2),
                   Y = modt$Y, ghQ = modt$ghQ, fam = modt$fam, bg = modt$b)
ana.res <- sche.test(b4,modt)
ana.score <- ana.res$gradient
ana.hess <- ana.res$hessian

eigen(-ana.hess,symmetric = T,only.values = T)$values

xlabn <- NULL
for(i in 1:ncol(modt$Y)){
  xlabn <- append(xlabn,paste0("mu[",i,",",1:ncol(modt$b$mu),"]"))
  if("sigma" %in% pFun(modt$fam[i])) xlabn <- append(xlabn, paste0("sigma[",i,",",1:ncol(modt$b$sigma),"]"))
  if("tau" %in% pFun(modt$fam[i])) xlabn <- append(xlabn, paste0("tau[",i,",",1:ncol(modt$b$tau),"]"))
  if("nu" %in% pFun(modt$fam[i])) xlabn <- append(xlabn, paste0("nu[",i,",",1:ncol(modt$b$nu),"]"))
}

par("mar" = c(6, 4, 4, 2) + 0.1)
plot(ana.score, main = "Score (first deriv. of log-likelihood) @ MLE", xlab = "", ylab = "Value", pch = 16,
     col = "gray50", xaxt = "n")
mtext(side = 1, text = "Parameter index", line = 4)
for(i in seq_along(xlabn)){ axis(1, at=i, labels = xlabn[i] , las = 2, cex.axis = 0.7) }
abline(h = 0, col = "forestgreen", lty = 2)
abline(v = seq(from = ((sum(length(modt$b$mu), length(modt$b$sigma)) %/% ncol(Y))),length.out = ncol(Y)-1,
               by = ((sum(length(modt$b$mu), length(modt$b$sigma)) %/% ncol(Y)))) + 0.5, lty = 3, col = "gray80")
points(num.score[1:length(ana.score)], col = "red", pch = 10)
segments(x0 = seq_along(ana.score), y0 = unlist(ana.score), y1 = num.score, col = "blue", lty = 3, lwd = 1)
legend("topleft",legend=c("Analytical", "Numerical"), col=c("gray50", "red"),
       pch = c(16,10), cex=1, inset = 0.02, box.col = "white")

int <- 1:sum(length(modt$b$mu), length(modt$b$sigma))^2
# int <- 1:((sum(length(modt$b$mu), length(modt$b$sigma))^2 %/% ncol(Y)) * 2)
plot(c(ana.hess)[int], main = "Hessian (2nd deriv. of log-likelihood) @ MLE",xlab ="Hessian entry Index", ylab = "Value",
     cex = 0.5, pch = 16, col = "gray50", ylim = c(min(ana.hess,num.hess),max(ana.hess,num.hess)))#, xaxt = "n")
abline(h = 0, col = "forestgreen", lty = 2);
abline(v = seq(from = sum(length(modt$b$mu), length(modt$b$sigma))^2 %/% ncol(Y),length.out = ncol(Y)-1,
               by = sum(length(modt$b$mu), length(modt$b$sigma))^2 %/% ncol(Y)) + 0.5, lty = 3, col = "gray80")
points(c(num.hess)[int], col = "red", pch = 16, cex = 0.5)
segments(x0 = seq_along(ana.hess[int]), y0 = unlist(ana.hess)[int], y1 = num.hess[int], col = "blue", lty = 3, lwd = 1)
legend(0,-2500,legend=c("Analytical", "Numerical (Hessian)"), col=c("gray50", "red",0.6),
       pch = c(16,16), cex=0.8, inset = 0.02, box.col = "white")
rm(modt)


rbenchmark::benchmark(
  "test20" = {testa <- splvm.fit(Y,fam,e.form,
          control = list(method = "EM", constraint = l1, full.hess = F, start.val = lc,
          ghQqp = 20, iter.lim = 200, tol = sqrt(.Machine$double.eps), silent = F))},
  "test40" = {testa <- splvm.fit(Y,fam,e.form,
          control = list(method = "EM", constraint = l1, full.hess = F, start.val = lc,
          ghQqp = 40, iter.lim = 200, tol = sqrt(.Machine$double.eps), silent = F))},
  replications = 3, columns = c("test", "replications", "elapsed", "relative", "user.self", "sys.self")
)

