
rm(list = ls())
set.seed(1234)
# rm(.Random.seed, envir=globalenv())

library(MASS)
library(abind)
library(splines)
library(mvtnorm)
library(numDeriv)
library(gamlss)

# `%!in%` <- Negate(`%in%`)

source("fa_main.R")
source("fb_auxf.R")
source("fc_mfit.R")
source("fd_mghq.R")
source("ff_msim.R")

n = 1000     # Number of individuals
p = 10       # Number of items
nsim = 1000  # Number of simulations
s.form <- list("mu" = "~ Z1 + Z2", "sigma" = "~ 1")
fam <- rep("normal",p)

l1 <- NULL
l1$mu <- matrix(1,ncol = 3, nrow = p)
l1$sigma <- matrix(1, ncol = 1, nrow = p)

l1. <- l1
# l1.$mu[1,3] <- 0
# l1.$mu[6,2] <- 0
# l1.$mu[sample(1:p, p/2, replace = F),3] <- 0
# l1.$mu[sample(1:p, p/2, replace = F),4] <- 0
# l1.$sigma[sample(1:p, p/2, replace = F),2] <- 0
# l1.$sigma[,3] <- 0
# l1.$sigma[sample(which(l1.$sigma[,2] == 1), sum(l1.$sigma[,2] == 1)/2, replace = F),3] <- 1

lc <- NULL
lc$mu <- matrix(runif(length(l1$mu),-1,1),nrow = p)
lc$mu[,1] <- runif(p,1,2) 
lc$mu[,2] <- runif(p,-1,1)
if(lc$mu[1,2] < 0) lc$mu[1,2] <- -lc$mu[1,2]
# lc$mu[,3] <- runif(p,-1,-0.5)
# lc$mu[,4] <- runif(p,-0.3,-0.1)
lc$sigma <- matrix(runif(length(l1$sigma), min = 0.1, max = 0.3), nrow = p)
if(ncol(lc$sigma) >= 2 && lc$sigma[1,2] < 0) lc$sigma[1,2] <- -lc$sigma[1,2]
# lc$sigma[,1] <- runif(p,-0.5,0.5)
# lc$sigma[,2] <- runif(p,-0.5,1)
# # lc$sigma[,3] <- runif(p,-0.3,-0.1)

lc. <- NULL
lc.$mu <- lc$mu * l1.$mu
lc.$sigma <- lc$sigma * l1.$sigma

simR <- splvm.sim(n = n, form = s.form, fam = fam, constraints = l1., coefs = lc.)
Y <- simR$Y
Z <- simR$Z
borg <- simR$b
e.form <- s.form
# e.form <- list("mu" = "~ bs(Z1, degree = 3)", "sigma" = "~ bs(Z1, degree = 3)")

# profvis::profvis({
testa0 <- splvm.fit(Y,fam,e.form,
          control = list(method = "EM", constraint = l1., full.hess = F, start.val = lc.,
          ghQqp = 10, iter.lim = 1e3, tol = sqrt(.Machine$double.eps), silent = F,information = "Fisher"))

testa1 <- splvm.fit(Y,fam,e.form,
          control = list(method = "EM", constraint = l1., full.hess = F, 
          ghQqp = 10, iter.lim = 1e3, tol = sqrt(.Machine$double.eps), silent = F,information = "Fisher"))

testa2 <- splvm.fit(Y,fam,e.form,
          control = list(method = "PEM", constraint = l1., full.hess = F, start.val = testa1$b,
          ghQqp = 10, iter.lim = 1e3, tol = sqrt(.Machine$double.eps), silent = F,information = "Fisher",
          pml.control = list(type = "scad", lambda = 0.01, w.alasso = testa1$b, gamma = 1, a = 3.7,pen.load = T)))

testa3 <- splvm.fit(Y,fam,e.form,
          control = list(method = "PEM", constraint = l1., full.hess = F, start.val = testa1$b,
          ghQqp = 10, iter.lim = 1e3, tol = sqrt(.Machine$double.eps), silent = F,information = "Fisher",
          pml.control = list(type = "mcp", lambda = 0.03, w.alasso = testa1$b, gamma = 1, a = 3.7,pen.load = T)))

round(cbind(testa1$b$mu,testa2$b$mu,testa3$b$mu,borg$mu),4)
round(cbind(testa1$b$sigma,testa2$b$sigma,testa3$b$sigma, borg$sigma),4)
# })

GBIC(testa1); GBIC(testa2); GBIC(testa3)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~

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

testb0 <- splvm.fit(Y,fam,e.form,
          control = list(method = "EM", constraint = l1, full.hess = T, start.val = lc, 
          ghQqp = 20, iter.lim = 1e3, tol = sqrt(.Machine$double.eps), silent = F, information = "Fisher"))

testb1 <- splvm.fit(Y,fam,e.form,
          control = list(method = "EM", constraint = l1, full.hess = T, #start.val = borg, 
          ghQqp = 20, iter.lim = 1e3, tol = sqrt(.Machine$double.eps), silent = F, information = "Fisher"))

testb2 <- splvm.fit(Y,fam,e.form,
          control = list(method = "PEM", constraint = l1, full.hess = T, start.val = testb0$b, 
          ghQqp = 20, iter.lim = 1e3, tol = sqrt(.Machine$double.eps), silent = F, information = "Fisher",
          pml.control = list(type = "mcp", lambda = 0.1, w.alasso = testb1$b,gamma = 1, a = 3.7,pen.load = F)))

testb3 <- splvm.fit(Y,fam,e.form,
          control = list(method = "PEM", constraint = l1, full.hess = T,
          ghQqp = 20, iter.lim = 1e3, tol = sqrt(.Machine$double.eps), silent = F, information = "Fisher",
          pml.control = list(type = "alasso", lambda = 0.01, w.alasso = testb1$b,gamma = 1, a = 3.7,pen.load = F)))

round(cbind(testb1$b$mu,testb2$b$mu,testb3$b$mu,borg$mu),4)
round(cbind(testb1$b$sigma,testb2$b$sigma,testb3$b$sigma, borg$sigma),4)

GBIC(testb1);GBIC(testb2); GBIC(testb3)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~

testc0 <- splvm.fit(Y,fam,e.form,
          control = list(method = "ML", constraint = l1, start.val = lc,
          ghQqp = 20, iter.lim = 1e3, tol = sqrt(.Machine$double.eps), silent = F, information = "Fisher"))

testc1 <- splvm.fit(Y,fam,e.form,
          control = list(method = "ML", constraint = l1, #start.val = borg, 
          ghQqp = 20, iter.lim = 1e3, tol = sqrt(.Machine$double.eps), silent = F, information = "Fisher"))

testc2 <- splvm.fit(Y,fam,e.form,
          control = list(method = "PML", constraint = l1, start.val = testc1$b, 
          ghQqp = 20, iter.lim = 1e3, tol = sqrt(.Machine$double.eps), silent = F, information = "Fisher",
          pml.control = list(type = "lasso", lambda = 0.2, w.alasso = testc1$b, gamma = 1, a = 3.7,pen.load = F)))

testc3 <- splvm.fit(Y,fam,e.form,
          control = list(method = "PML", constraint = l1, start.val = testc1$b,
          ghQqp = 20, iter.lim = 1e3, tol = sqrt(.Machine$double.eps), silent = F, information = "Fisher",
          pml.control = list(type = "alasso", lambda = 0.01, w.alasso = testc1$b, gamma = 1, a = 3.7,pen.load = F)))

round(cbind(testc1$b$mu,testc2$b$mu,testc3$b$mu,borg$mu),3)
round(cbind(testc1$b$sigma,testc2$b$sigma,testc3$b$sigma, borg$sigma),4)
testc1$loglik; testc2$loglik; testc3$loglik
GBIC(testc1);GBIC(testc2); GBIC(testc3)
 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~

testd0 <- splvm.fit(Y,fam,e.form,
          control = list(method = "hybrid", constraint = l1, full.hess = F, start.val = lc, 
          EM.iter.lim = 10, ghQqp = 20, iter.lim = 1e3, tol = sqrt(.Machine$double.eps), silent = F, information = "Fisher"))

testd1 <- splvm.fit(Y,fam,e.form,
          control = list(method = "hybrid", constraint = l1, full.hess = F,
          EM.iter.lim = 10, ghQqp = 20, iter.lim = 1e3, tol = sqrt(.Machine$double.eps), silent = F, information = "Fisher"))

testd2 <- splvm.fit(Y,fam,e.form,
          control = list(method = "P-hybrid", constraint = l1, full.hess = F, start.val = testd1$b,
          EM.iter.lim = 10, ghQqp = 20, iter.lim = 1e3, tol = sqrt(.Machine$double.eps), silent = F, information = "Fisher",
          pml.control = list(type = "lasso", lambda = 0.01, w.alasso = testd1$b, gamma = 1, a = 3.7,pen.load = F)))

testd3 <- splvm.fit(Y,fam,e.form,
          control = list(method = "P-hybrid", constraint = l1, full.hess = F,start.val = testd1$b,
          EM.iter.lim = 10, ghQqp = 20, iter.lim = 1e3, tol = sqrt(.Machine$double.eps), silent = F, information = "Fisher",
          pml.control = list(type = "alasso", lambda = 0.01, w.alasso = testd1$b, gamma = 1, a = 3.7,pen.load = F)))

round(cbind(testd1$b$mu,testd2$b$mu,testd3$b$mu,borg$mu),4)
round(cbind(testd1$b$sigma,testd2$b$sigma,testd3$b$sigma, borg$sigma),4)

GBIC(testd1); GBIC(testd2); GBIC(testd3)
 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~

round(cbind(testa0$b$mu, testb0$b$mu, testc0$b$mu, testd0$b$mu, borg$mu),3)
round(cbind(testa1$b$mu, testb1$b$mu, testc1$b$mu, testd1$b$mu, borg$mu),3)
round(cbind(testa2$b$mu, testb2$b$mu, testc2$b$mu, testd2$b$mu, borg$mu),3)
round(cbind(testa3$b$mu, testb3$b$mu, testc3$b$mu, testd3$b$mu, borg$mu),3)
round(cbind(testa1$b$sigma, testb1$b$sigma, testc1$b$sigma, testd1$b$sigma, borg$sigma),3)
round(cbind(testa2$b$sigma, testb2$b$sigma, testc2$b$sigma, testd2$b$sigma, borg$sigma),3)
round(cbind(testa3$b$sigma, testb3$b$sigma, testc3$b$sigma, testd3$b$sigma, borg$sigma),3)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Test for different values of lambda in PML
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

lambda <- seq(0.02,0.005,length.out = 7)

test <- NULL
for(i in seq_along(lambda)){
 test[[i]] <- splvm.fit(Y,fam,e.form,
           control = list(method = "PEM", constraint = l1,full.hess = T,
           ghQqp = 20, iter.lim = 1e3, tol = sqrt(.Machine$double.eps), silent = T,information = "Fisher",
           pml.control = list(type = "lasso", lambda = lambda[i]),pen.load = T))
}

round(cbind(test[[which(lambda == 0.0150)]]$b$mu,
            test[[which(lambda == 0.0125)]]$b$mu,
            test[[which(lambda == 0.0100)]]$b$mu,
            test[[which(lambda == 0.0075)]]$b$mu,
            test[[which(lambda == 0.0050)]]$b$mu),3)

c(GBIC(test[[which(lambda == 0.0150)]]),
GBIC(test[[which(lambda == 0.0125)]]),
GBIC(test[[which(lambda == 0.0100)]]),
GBIC(test[[which(lambda == 0.0075)]]), 
GBIC(test[[which(lambda == 0.0050)]]))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Test: Comparison vs penfa & lavaan
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 
library(penfa)
data(ccdata)

ccdata.subset <- ccdata[ccdata$country == "LEB",2:8]
### Single-group analysis (no mean-structure, unit factor variances)
syntax = 'Z1 =~ h1 + h2 + h3 + h4 + h5 + h6 + h7 '
 # h1 + h2 + h3 + h4 + h5 + h6 + h7 ~ 1'
alasso_fit <- penfa(model = syntax, data = ccdata.subset, std.lv = TRUE,
pen.shrink = "alasso", information = "hessian",
eta = list(shrink = c("lambda" = 0.1), diff = c("none" = 0)),
strategy = "auto", gamma = 3.7, meanstructure = F)
summary(alasso_fit)
penfaParEstim(alasso_fit)
penmat(alasso_fit)

library(lavaan)
CFA.model <- ' Z1 =~ h1 + h2 + h3 + h4 + h5 + h6 + h7 '
rmtest <- cfa(CFA.model, data = ccdata.subset, orthogonal = T, meanstructure = F, std.lv = TRUE, estimator="ML",
              information = "observed"); # coef(rmtest)
summary(rmtest)
# rm(CFA.model,rmtest)
 
fam.ex <- rep("normal",ncol(ccdata.subset))
ex.form <- list("mu" = "~ Z1-1", "sigma" = "~ 1")
test <- t(t(ccdata.subset)-colMeans(ccdata.subset)) # test <- ccdata.subset
ex.test.unp <- splvm.fit(test,fam.ex,ex.form,
           control = list(method = "EM", full.hess = T,#start.val = ex.ulb,
           ghQqp = 100, iter.lim = 1e3, tol = .Machine$double.eps^(1/2), silent = F,information = "Hessian"))

# rmtes2 <- ERP::emfa(as.matrix(ccdata.subset), 1, min.err = .Machine$double.eps, verbose = FALSE, svd.method = c("irlba")) #"irlba" "fast.svd"
# ex.EMfa <- NULL
# ex.EMfa$mu <- matrix(rmtes2$B, dimnames = dimnames(ex.test.unp$b$mu))
# ex.EMfa$sigma <- matrix(log(sqrt(rmtes2$Psi)),dimnames = dimnames(ex.test.unp$b$sigma))
# # ~~~~~~~~~~~~~
ex.ulb <- mb2lb(matrix(coef(rmtest),nrow = ncol(test),byrow = F),ex.test.unp$b) #[,c(3,1,2)]
ex.ulb$sigma <- log(sqrt(ex.ulb$sigma))
c(logLik(rmtest)); ex.test.unp$loglik
round(cbind(ex.ulb$mu,ex.test.unp$b$mu),3)
round(cbind(ex.ulb$sigma,ex.test.unp$b$sigma),3)
 
ex.lb <- mb2lb(matrix(coef(alasso_fit),nrow = ncol(test),byrow = F),ex.test.unp$b) # [,c(2,1,3)]
ex.lb$sigma <- log(sqrt(ex.lb$sigma))

ex.test.pen <- splvm.fit(test,fam.ex,ex.form,
           control = list(method = "PEM", full.hess = T, start.val = ex.lb,
           ghQqp = 100, iter.lim = 1e3, tol = .Machine$double.eps^(1/2), silent = F,information = "Hessian",
           pml.control = list(type = "alasso", w.alasso = ex.lb, lambda = alasso_fit@Options$eta$shrink, pen.load = T)))

mod <- NULL
mod$b <- ex.lb
mod$Y <- ex.test.pen$Y
mod$ghQ <- ex.test.pen$ghQ
mod$fam <- ex.test.pen$fam
mod$pml.control <- ex.test.pen$pml.control
mod$loglik <- llkf(lb2cb(mod$b),mod$Y,mod$ghQ,mod$b,mod$fam) # Difference with penfa; c(alasso_fit@Optim$logl.unpen)
mod$info = "Hessian"
mod$hessian <- sche.test(lb2cb(mod$b),mod)$hessian
# mod$hessian <- unname(-alasso_fit@Optim$hessian.pen[a,a])

round(cbind(ex.lb$mu,ex.test.pen$b$mu),3)
round(cbind(ex.lb$sigma,ex.test.pen$b$sigma),3)

round(cbind(ex.test.unp$b$mu,ex.test.pen$b$mu,ex.lb$mu),3)
round(cbind(exp(ex.test.unp$b$sigma),exp(ex.test.pen$b$sigma),exp(ex.lb$sigma))^2,3)
# coef(alasso_fit)
c(alasso_fit@Optim$logl.unpen);ex.test.unp$loglik
c(alasso_fit@Optim$logl.pen); ex.test.pen$loglik; llkf(lb2cb(ex.test.pen$b),ex.test.unp$Y,ex.test.unp$ghQ,ex.test.unp$b,ex.test.unp$fam)

pen.idx <- pidx(ex.test.pen$b,T)
penmat(alasso_fit)[1:7,1:7]
nrow(ex.test.pen$Y)*penM(lb2mb(ex.test.pen$b)[pen.idx],type = ex.test.pen$pml.control$type, lambda = ex.test.pen$pml.control$lambda,
                  w.alasso = lb2mb(ex.test.pen$pml.control$w.alasso)[pen.idx], gamma = ex.test.pen$pml.control$gamma,
                  a = ex.test.pen$pml.control$a)
nrow(ex.test.pen$Y)*penM(lb2mb(ex.lb)[pen.idx],type = ex.test.pen$pml.control$type, lambda = ex.test.pen$pml.control$lambda,
                  w.alasso = lb2mb(ex.test.pen$pml.control$w.alasso)[pen.idx], gamma = ex.test.pen$pml.control$gamma,
                  a = ex.test.pen$pml.control$a)
# slotNames(alasso_fit)

rmMat <- mod$hessian
rmMat2 <- ex.test.pen$hessian
rmMat3 <- ex.test.unp$hessian
# a <- NULL; for(i in 1:ncol(test)){a <- c(a,i+ncol(test),i,i+2*ncol(test))}
a <- NULL; for(i in 1:ncol(test)){a <- c(a,i,i+ncol(test))}
plot(c(rmMat2), col = "black", pch = 16, cex = 0.5) # mine
points(c(rmMat), col = "blue", pch = 16, cex = 0.5) # Implied Elena's
points(c(-alasso_fit@Optim$hessian.pen[a,a]),col = "red", pch = 16, cex = 0.5) # Elena's
rm(rmMat,rmMat2,a)

GBIC(ex.test.pen); GIC(ex.test.pen)
GBIC(ex.test.unp); GIC(ex.test.unp)
GBIC(mod); GIC(mod)
summary(alasso_fit,estimates = F,nd = 3L,extra = TRUE)

# ~~~~~~~~~~~~~~~~~~~~~~~
# Test: Comparison vs ltm
# ~~~~~~~~~~~~~~~~~~~~~~~

library(ltm)
test <- WIRS
test <- Mobility
# data(test)
rmltm <- ltm::ltm(test ~ z1*z2, IRT.param = F,control = list(GHk = 15, iter.em = 350))
# rmltm <- ltm(test ~ z1 + I(z1^2), IRT.param = F,control = list(GHk = 50, iter.em = 350))
# coef(rmltm)
# rmltm$log.Lik

ex.lb <- mb2lb(matrix(coef(rmltm),nrow = ncol(test),byrow = F),ex.test.unp$b)

fam.ex <- rep("binomial",ncol(test))
ex.form <- list("mu" = "~ Z1*Z2")
# ex.form <- list("mu" = "~ Z1 + I(Z1^2)")
ex.test.unp <- splvm.fit(test,fam.ex,ex.form,
           control = list(method = "EM", full.hess = F, start.val = ex.lb,
           ghQqp = 15, iter.lim = 1e3, tol = .Machine$double.eps^0.5, silent = F,information = "Fisher"))

round(cbind(ex.test.unp$b$mu,ex.lb$mu),3)
rmltm$log.Lik; ex.test.unp$loglik

ex.test.pen <- splvm.fit(test,fam.ex,ex.form,
           control = list(method = "PEM", full.hess = F,start.val = ex.test.unp$b,
           ghQqp = 15, iter.lim = 1e3, tol = .Machine$double.eps^0.5, silent = F,information = "Fisher",
           pml.control = list(type = "scad", w.alasso = ex.test.unp$b, lambda = 0.1, pen.load = T)))

round(cbind(ex.test.pen$b$mu,ex.test.unp$b$mu,ex.lb$mu),3)
rmltm$log.Lik; ex.test.unp$loglik; ex.test.pen$loglik

GBIC(ex.test.pen)
GBIC(ex.test.unp)

mod <- NULL
mod$b <- ex.lb
mod$Y <- ex.test.unp$Y
mod$ghQ <- ex.test.unp$ghQ
mod$fam <- ex.test.unp$fam
mod$loglik <- llkf(lb2cb(mod$b),mod$Y,mod$ghQ,mod$b,mod$fam) # rmltm$log.Lik


# # ~~~~~~~~~~~~~~~~~~~~~~~
# # ~~~~~~~~~~~~~~~~~~~~~~~
# # ~~~~~~~~~~~~~~~~~~~~~~~
# # ~~~~~~~~~~~~~~~~~~~~~~~
# # ~~~~~~~~~~~~~~~~~~~~~~~
# # ~~~~~~~~~~~~~~~~~~~~~~~

zsco <- fscore(testc1) # check function

modt <- ex.test.pen
b1 <- lb2cb(modt$b)
b2 <- lb2cb(lc)
b3 <- rep(1,length(b1))*lb2cb(modt$loadmt)
b4 <- lb2cb(lc.)
num.score <- grad(func = llkf, x = b1, method = "Richardson", method.args=list(r = 6, v = 2),
                   Y = modt$Y, ghQ = modt$ghQ, fam = modt$fam, bg = modt$b)
num.hess <- hessian(func = llkf, x = b1, method = "Richardson", method.args=list(r = 6, v = 2),
                   Y = modt$Y, ghQ = modt$ghQ, fam = modt$fam, bg = modt$b)
ana.res <- sche.test(b1,modt)
ana.score <- ana.res$gradient
ana.hess <- ana.res$hessian

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


hist(Z$Z1, breaks = 100, freq = F, xlim = c(-5,5), border = "gray", ylim = c(0,0.7))
hist(Z$Z2, breaks = 100, freq = F, xlim = c(-5,5), border = "gray", ylim = c(0,0.7))
# lines(density(rmtestfac$score.dat$z1), col = 5, lwd = 2)
lines(density(Z$Z1), col = 1, lwd = 2)
lines(density(Z$Z2), col = 1, lwd = 2)
lines(density(fscore(testa0)$Z1), col = 2, lwd = 2)
lines(density(fscore(testa0)$Z2), col = 2, lwd = 2)
plot(Z$Z1, fscore(testa0)$Z1); abline(0,1, col = "red", lty = 3)
plot(Z$Z2, fscore(testa0)$Z2); abline(0,1, col = "red", lty = 3)
# plot(Z$Z2, fscore(testa0)$Z2)
cor(Z$Z1, fscore(testa0)$Z1,method = "kendall")
cor(Z$Z1, fscore(testa0)$Z1,method = "pearson")
cor(Z$Z2, fscore(testa0)$Z2,method = "pearson")
points(testa0$ghQ$points, testa0$ghQ$weights, col = "red", pch = 16)
lines(seq(-5,5,length.out = 100), dnorm(seq(-5,5,length.out = 100)), col = "blue", lwd = 2)
