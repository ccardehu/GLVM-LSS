
rm(list = ls())
set.seed(1234)
# # rm(.Random.seed, envir=globalenv())
# 
# library(MASS)
# library(abind)
# library(splines)
# library(mvtnorm)
# library(numDeriv)
# library(gamlss)
# library(foreach)
# library(doSNOW)
# library(parallel)
# library(itertools)
# library(xtable)
# 
# # `%!in%` <- Negate(`%in%`)
# 
# source("fa_main.R")
# source("fb_auxf.R")
# source("fc_mfit.R")
# source("fd_mghq.R")
# source("ff_msim.R")
# source("fg_grph.R")
# source("fi_esim.R")

source("f0_prep.R")

n = 1000     # Number of individuals
nsim = 1000  # Number of simulations

# Simulation 1: Heteroscedastic Normal model:
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# s.form <- list("mu" = "~ Z1 + Z2", "sigma" = "~ Z1 + Z2"); p = 10
# fam <- rep("normal",p)
# l1 <- lc <- NULL
# l1$mu <- lc$mu <- matrix(1,ncol = 3, nrow = p)
# l1$sigma <- lc$sigma <- matrix(1, ncol = 3, nrow = p)
# lc$mu[,1] <- runif(p,1,2)
# lc$mu[,c(2,3)] <- runif(length(l1$mu[,c(2,3)]),-1.5,1.5)# * sample(c(-1,1),size = 2*p,replace = T)
# if(lc$mu[1,2] < 0) lc$mu[1,2] <- -lc$mu[1,2]; if(lc$mu[2,3] < 0) lc$mu[2,3] <- -lc$mu[2,3]
# lc$sigma <- matrix(runif(length(l1$sigma),0.1,0.4), nrow = p)
# lc$sigma[,c(2,3)] <- lc$sigma[,c(2,3)]# * sample(c(-1,1),size = 2*p,replace = T)
# if(ncol(lc$sigma) >= 2 && lc$sigma[1,2] < 0) lc$sigma[1,2] <- -lc$sigma[1,2]
# if(ncol(lc$sigma) >= 2 && lc$sigma[2,3] < 0) lc$sigma[2,3] <- -lc$sigma[2,3]


# Simulation 2: IRT model:
# ~~~~~~~~~~~~~~~~~~~~~~~~
s.form <- list("mu" = "~ Z1+Z2"); p = 10
fam <- rep("binomial",p)
l1 <- lc <- NULL
l1$mu <- lc$mu <- matrix(1,ncol = 3, nrow = p)
lc$mu[,1] <- runif(length(l1$mu[,1]),-1,1)
lc$mu[,c(2,3)] <- runif(length(l1$mu[,c(2,3)]),1.5,2.5) * sample(c(-1,1),size = 2*p,replace = T)
if(lc$mu[1,2] < 0) lc$mu[1,2] <- -lc$mu[1,2]; if(lc$mu[2,3] < 0) lc$mu[2,3] <- -lc$mu[2,3]


# Simulation 3: ZI-Poisson:
# ~~~~~~~~~~~~~~~~~~~~~~~~~
# s.form <- list("mu" = "~ Z1 + Z2", "sigma" = "~ Z1 + Z2"); p = 15
# fam <- rep("ZIpoisson", p)
# l1 <- lc <- NULL
# l1$mu <- matrix(1,ncol = 3, nrow = p)
# l1$sigma <- matrix(1, ncol = 3, nrow = p)
# lc$mu <- matrix(runif(length(l1$mu),0.5,1),nrow = p)
# lc$sigma <- matrix(runif(length(l1$sigma),1.5,2.5), nrow = p)
# lc$mu[,1] <- runif(p,0.5,1); if(lc$mu[1,2] < 0) lc$mu[1,2] <- -lc$mu[1,2] ; if(lc$mu[1,3] < 0) lc$mu[1,3] <- -lc$mu[1,3]
# lc$sigma[,1] <- runif(p,-1,1)

# ~~~~~~~~~~~
# Starts here:
# ~~~~~~~~~~~

l1. <- l1

# Parameter restrictions:
# Simulation 1
# ~~~~~~~~~~~~
# l1.$mu[1,c(3)] <- l1.$mu[2,c(2)] <- 0
# l1.$mu[sample(3:p, floor((p-2)/2), replace = F),2] <- 0
# for(i in 3:nrow(l1.$mu)){ if(l1.$mu[i,2] != 0) l1.$mu[i,3] <- rbinom(1,1,0.5) }
# l1.$sigma[1,3] <- l1.$sigma[2,2] <- 0
# l1.$sigma[sample(3:p, floor((p-2)/2), replace = F),2] <- 0
# for(i in 3:nrow(l1.$sigma)){ if(l1.$sigma[i,2] == 0) l1.$sigma[i,3] <- rbinom(1,1,0.5) }

# Simulation 2
# ~~~~~~~~~~~~
l1.$mu[1,c(3)] <- l1.$mu[2,c(2)] <- 0
l1.$mu[sample(3:p, floor((p-2)/2), replace = F),2] <- 0
for(i in 3:nrow(l1.$mu)){ if(l1.$mu[i,2] != 0) l1.$mu[i,3] <- rbinom(1,1,0.5) }

# Simulation 3
# ~~~~~~~~~~~~
# l1.$mu[1,3] <- l1.$mu[2,2] <- 0
# l1.$mu[sample(3:p, floor((p-3)/3), replace = F),2] <- 0
# for(i in 3:nrow(l1.$mu)){ if(l1.$mu[i,2] != 0) l1.$mu[i,3] <- rbinom(1,1,0.5) } 
# l1.$sigma[1,3] <- l1.$sigma[2,2] <- 0
# l1.$sigma[sample(3:p, floor((p-3)/3), replace = F),2] <- 0
# for(i in 3:nrow(l1.$sigma)){ if(l1.$sigma[i,2] != 0) l1.$sigma[i,3] <- rbinom(1,1,0.5) } 

# Simulation procedure:
# ~~~~~~~~~~~~~~~~~~~~~
lc. <- NULL; for(i in names(l1.)){ lc.[[i]] <- lc[[i]]*l1.[[i]] }
simR <- splvm.sim(n = n, form = s.form, fam = fam, constraints = l1., coefs = lc.)
Y <- simR$Y
Z <- simR$Z
borg <- simR$b
e.form <- s.form

# library(xtable)
# xtable(unname(borg$sigma), digits = 3)

# Restrictions (l1.)
# ~~~~~~~~~~~~~~~~~~ 

# Simulation 1: Quadratic heteroscedastic model
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# restr <- list(c("mu",2,"Z1",0), c("mu",3,"Z1",0), c("mu",5,"Z1",0), c("mu",8,"Z1",0), c("mu",10,"Z1",0),
#               c("mu",1,"Z2",0), c("mu",6,"Z2",0), c("mu",9,"Z2",0),
#               c("sigma",2,"Z1",0), c("sigma",3,"Z1",0), c("sigma",4,"Z1",0), c("sigma",5,"Z1",0), c("sigma",10,"Z1",0),
#               c("sigma",1,"Z2",0), c("sigma",3,"Z2",0), c("sigma",5,"Z2",0), c("sigma",10,"Z2",0))
# 
# irestr <- list(c("mu",1,"Z2",0), c("mu",2,"Z1",0), c("sigma",1,"Z2",0), c("sigma",2,"Z1",0))

# Simulation 2: Interaction binomial
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
restr <- list(c("mu",1,"Z2",0), c("mu",3,"Z2",0), c("mu",8,"Z2",0),
              c("mu",2,"Z1",0), c("mu",4,"Z1",0), c("mu",5,"Z1",0), c("mu",6,"Z1",0), c("mu",9,"Z1",0))

irestr <- list(c("mu",1,"Z2",0), c("mu",2,"Z1",0))
               

# Simulation 3: ZI-Poisson
# restr <- list(c("mu",2,"Z1",0), c("mu",3,"Z1",0), c("mu",6,"Z1",0), c("mu",8,"Z1",0), c("mu",11,"Z1",0),
#               c("mu",1,"Z2",0), c("mu",5,"Z2",0), c("mu",7,"Z2",0), c("mu",12,"Z2",0), c("mu",13,"Z2",0), c("mu",15,"Z2",0),
#               c("sigma",2,"Z1",0), c("sigma",4,"Z1",0), c("sigma",7,"Z1",0), c("sigma",8,"Z1",0), c("sigma",12,"Z1",0),
#               c("sigma",1,"Z2",0), c("sigma",5,"Z2",0), c("sigma",6,"Z2",0), c("sigma",9,"Z2",0))
# 
# irestr <- list(c("mu",1,"Z2",0), c("mu",2,"Z1",0), c("sigma",1,"Z2",0), c("sigma",2,"Z1",0))

testCFa <- splvm.fit(Y,fam,e.form,
           control = list(method = "EM",constraint = restr, silent = F))
testCFb <- splvm.fit(Y,fam,e.form,
           control = list(method = "ML",constraint = restr, silent = F, information = "Hessian"))

testUPa <- splvm.fit(Y,fam,e.form,control = list(method = "EM", silent = F, constraint = irestr))
testUPb <- splvm.fit(Y,fam,e.form,control = list(method = "ML", silent = F, information = "Hessian", constraint = irestr))

testALa <- splvm.fit(Y,fam,e.form, control = list(method = "PEM", constraint = irestr,
           silent = F, pml.control = list(type = "alasso", lambda = "auto", gamma = 2,
           w.alasso = testUPa$b)))
testALb <- splvm.fit(Y,fam,e.form, control = list(method = "PML", information = "Hessian", constraint = irestr,
           silent = F, pml.control = list(type = "alasso", lambda = "auto", gamma = 2,
           w.alasso = testUPb$b)))

# testMCa <- splvm.fit(Y,fam,e.form, control = list(method = "PEM",constraint = irestr,
#            silent = F, pml.control = list(type = "mcp", lambda = 0.05)))
# testMCb <- splvm.fit(Y,fam,e.form, control = list(method = "PML", information = "Hessian", constraint = irestr,
#            silent = F, pml.control = list(type = "mcp", lambda = 0.05)))

round(c(testCFa$log,testUPa$log,testALa$log),3)
round(c(testCFb$log,testUPb$log,testALb$log),3)
for(i in names(s.form)){ 
  print(round(cbind(testCFa$b[[i]],testUPa$b[[i]],testALa$b[[i]],borg[[i]]),4)) }
for(i in names(s.form)){ 
  print(round(cbind(testCFb$b[[i]],testUPb$b[[i]],testALb$b[[i]],borg[[i]]),4)) }

GBIC(testCFa); GBIC(testUPa); GBIC(testALa);# GBIC(testMCa)
GBIC(testCFb); GBIC(testUPb); GBIC(testALb);# GBIC(testMCb)
GIC(testCFa); GIC(testUPa); GIC(testALa);# GIC(testMCa)
GIC(testCFb); GIC(testUPb); GIC(testALb);# GIC(testMCb)

# syntax = 'Z1 =~ Y1 + Y2 + Y3 + Y4 + Y5 + Y6 + Y7 + Y8 + Y9 + Y10 ;
# Z2 =~ Y1 + Y2 + Y3 + Y4 + Y5 + Y6 + Y7 + Y8 + Y9 + Y10 ;
# Y1 + Y2 + Y3 + Y4 + Y5 + Y6 + Y7 + Y8 + Y9 + Y10 ~ 1 ;
# Z1 ~~ 1*Z2'
# # library(penfa)
# alasso_fit <- penfa(model = syntax, data = Y, std.lv = TRUE,
# pen.shrink = "lasso", information = "fisher",
# eta = list(shrink = c("lambda" = 0.1), diff = c("none" = 0)),
# strategy = "auto", meanstructure = T)
# summary(alasso_fit)
# penfaParEstim(alasso_fit)
# penmat(alasso_fit)
# alasso_fit@Options$eta$shrink
# matrm <- array(0,dim = dim(lb2mb(testa2$loadmt)))
# matrm[lb2mb(testa2$loadmt)[,c(2,3,1,4)]] <- coef(alasso_fit)
# ex.lb <- mb2lb(matrm[,c(3,1,2,4)],testa1$b) # [,c(2,1,3)]
# ex.lb <- mb2lb(matrix(coef(alasso_fit),nrow = ncol(Y),byrow = F)[,c(2,3,1,3)],testa1$b) # [,c(2,1,3)]
# ex.lb$sigma <- log(sqrt(ex.lb$sigma))
# round(cbind(testa2$b$mu,ex.lb$mu,borg$mu),3)
 
# round(cbind(testa1$b$mu,testa2$b$mu,borg$mu),4)
# round(cbind(testa1$b$sigma,testa2$b$sigma,borg$sig),3)



####################################
# Test: Comparison vs penfa & lavaan
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 
library(penfa)
data(ccdata)

ccdata.subset <- ccdata[ccdata$country == "LEB",2:8]
### Single-group analysis (no mean-structure, unit factor variances)
syntax = 'Z1 =~ h1 + h2 + h3 + h4 + h5 + h6 + h7 
h1 + h2 + h3 + h4 + h5 + h6 + h7 ~ 1'
alasso_fit <- penfa(model = syntax, data = ccdata.subset, std.lv = TRUE,
pen.shrink = "alasso", information = "fisher",
eta = list(shrink = c("lambda" = 0.1), diff = c("none" = 0)),
strategy = "auto", meanstructure = T)
summary(alasso_fit)
penfaParEstim(alasso_fit)
penmat(alasso_fit)

library(lavaan)
CFA.model <- ' Z1 =~ h1 + h2 + h3 + h4 + h5 + h6 + h7 '
rmtest <- cfa(CFA.model, data = ccdata.subset, orthogonal = T, meanstructure = T, std.lv = TRUE, estimator="ML",
              information = "observed"); # coef(rmtest)

fam.ex <- rep("normal",ncol(ccdata.subset))
# ex.form <- list("mu" = "~ Z1-1", "sigma" = "~ 1")
# test <- t(t(ccdata.subset)-colMeans(ccdata.subset)) 
ex.form <- list("mu" = "~ Z1", "sigma" = "~ 1")
test <- ccdata.subset

ex.test.unp <- splvm.fit(test,fam.ex,ex.form,control = list(method = "EM", silent = F))

# rmtes2 <- ERP::emfa(as.matrix(ccdata.subset), 1, min.err = .Machine$double.eps, verbose = FALSE, svd.method = c("irlba")) #"irlba" "fast.svd"
# ex.EMfa <- NULL
# ex.EMfa$mu <- matrix(rmtes2$B, dimnames = dimnames(ex.test.unp$b$mu))
# ex.EMfa$sigma <- matrix(log(sqrt(rmtes2$Psi)),dimnames = dimnames(ex.test.unp$b$sigma))

# # ~~~~~~~~~~~~~
ex.ulb <- mb2lb(matrix(coef(rmtest),nrow = ncol(test),byrow = F)[,c(3,1,2)],ex.test.unp$b) #[,c(3,1,2)]
ex.ulb$sigma <- log(sqrt(ex.ulb$sigma))
c(logLik(rmtest)); ex.test.unp$loglik
round(cbind(ex.ulb$mu,ex.test.unp$b$mu),3)
round(cbind(ex.ulb$sigma,ex.test.unp$b$sigma),3)
 
ex.lb <- mb2lb(matrix(coef(alasso_fit),nrow = ncol(test),byrow = F)[,c(2,1,3)],ex.test.unp$b) # [,c(2,1,3)]
ex.lb$sigma <- log(sqrt(ex.lb$sigma))

ex.test.pen <- splvm.fit(test,fam.ex,ex.form,control = list(method = "PEM",
               silent = F, pml.control = list(type = "alasso", w.alasso = ex.lb, lambda = "auto", pen.load = T)))

# alasso_fit@Options$eta$shrink
# Y = test; fam = fam.ex; form = ex.form
# ex.test.pen$pml.control$lambda

# Penalised comparison
round(cbind(ex.lb$mu,ex.test.pen$b$mu),3)
round(cbind(ex.lb$sigma,ex.test.pen$b$sigma),3)

# Unpenalised comparison
round(cbind(ex.ulb$mu,ex.test.unp$b$mu),3)
round(cbind(ex.ulb$sigma,ex.test.unp$b$sigma),3)

# Penalised vs. Unpenalised
round(cbind(lb2mb(ex.ulb),lb2mb(ex.test.unp$b)),3)
round(cbind(lb2mb(ex.lb),lb2mb(ex.test.pen$b)),3)

GBIC(ex.test.pen); GIC(ex.test.pen)
GBIC(ex.test.unp); GIC(ex.test.unp)
ex.test.unp$loglik; ex.test.pen$loglik

c(alasso_fit@Optim$logl.unpen);ex.test.unp$loglik
c(alasso_fit@Optim$logl.pen); ex.test.pen$loglik;

#########################
# Test: Comparison vs ltm
# ~~~~~~~~~~~~~~~~~~~~~~~

library(ltm)
test <- WIRS
# test <- Mobility
# data(test)
rmltm <- ltm::ltm(test ~ z1*z2, IRT.param = F,control = list(GHk = 15, iter.em = 350))
# rmltm <- ltm(test ~ z1 + I(z1^2), IRT.param = F,control = list(GHk = 50, iter.em = 350))
# coef(rmltm)
# rmltm$log.Lik
 
fam.ex <- rep("binomial",ncol(test))
ex.form <- list("mu" = "~ Z1*Z2")
# ex.form <- list("mu" = "~ Z1 + I(Z1^2)")
ex.test.unp <- splvm.fit(test,fam.ex,ex.form, control = list(method = "EM", ghQp = 6,
               silent = F,information = "Fisher", iter.lim = 2e3))

ex.lb <- mb2lb(matrix(coef(rmltm),nrow = ncol(test),byrow = F),ex.test.unp$b)

# round(cbind(ex.test.unp$b$mu,ex.lb$mu),3)
# rmltm$log.Lik; ex.test.unp$loglik

ex.test.pen <- splvm.fit(test,fam.ex,ex.form,
           control = list(method = "PEM", silent = F,information = "Fisher", ghQp = 6, iter.lim = 2e3,
           pml.control = list(type = "alasso", w.alasso = ex.test.unp$b, lambda = "auto", pen.load = T)))

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

#######################
# # ~~~~~~~~~~~~~~~~~~~~~~~
# # ~~~~~~~~~~~~~~~~~~~~~~~
# # ~~~~~~~~~~~~~~~~~~~~~~~
# # ~~~~~~~~~~~~~~~~~~~~~~~
# # ~~~~~~~~~~~~~~~~~~~~~~~
# # ~~~~~~~~~~~~~~~~~~~~~~~

zsco <- fscore(testc1) # check function

modt <- testUP
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
legend(0,-40000,legend=c("Analytical", "Numerical (Hessian)"), col=c("gray50", "red",0.6),
       pch = c(16,16), cex=0.8, inset = 0.02, box.col = "white")
rm(modt)

#######################

testm <- testALb

hist(Z$Z1, breaks = 100, freq = F, xlim = c(-5,5), border = "gray", ylim = c(0,0.7))
hist(Z$Z2, breaks = 100, freq = F, xlim = c(-5,5), border = "gray", ylim = c(0,0.7))
# lines(density(rmtestfac$score.dat$z1), col = 5, lwd = 2)
lines(density(Z$Z1), col = 1, lwd = 2)
lines(density(Z$Z2), col = 1, lwd = 2)
lines(density(fscore(testm)$Z1), col = 2, lwd = 2)
lines(density(fscore(testm)$Z2), col = 2, lwd = 2)
plot(Z$Z1, fscore(testm)$Z1); abline(0,1, col = "red", lty = 3)
plot(Z$Z2, fscore(testm)$Z2); abline(0,1, col = "red", lty = 3)
# plot(Z$Z2, fscore(testa0)$Z2)
cor(Z$Z1, fscore(testm)$Z1,method = "kendall")
cor(Z$Z1, fscore(testm)$Z1,method = "pearson")
cor(Z$Z2, fscore(testm)$Z2,method = "pearson")
points(testm$ghQ$points, testm$ghQ$weights, col = "red", pch = 16)
lines(seq(-5,5,length.out = 100), dnorm(seq(-5,5,length.out = 100)), col = "blue", lwd = 2)

##################
#########################
##########################

Ex1 <- splvm.mcsim(3, sil = F, n = n, form = s.form, fam = fam, restr = l1, coefs = lc,
                   saveRes = F, cleanRes = F)
splvm.plot(1,testa1,
           control.plot = list(sim.mod = simR,
           plot.mean = T, plot.sd = F, qtl = c(0.025,0.975), plot.3D = T, bw = F))

splvm.plot(1,testa3,
           control.plot = list(sim.mod = simR,
           plot.mean = T, plot.sd = F, qtl = c(seq(0.1,0.9,by = 0.2)), plot.3D = T, bw = F))
