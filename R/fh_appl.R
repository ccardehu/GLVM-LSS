rm(list = ls())
set.seed(1234)
# rm(.Random.seed, envir=globalenv())

source("fa_main.R")
source("fb_auxf.R")
source("fc_mfit.R")
source("fd_mghq.R")
source("ff_msim.R")
source("fg_grph.R")

library(MASS)
library(abind)
library(splines)
library(mvtnorm)
library(numDeriv)
library(gamlss)
library(haven)

# `%!in%` <- Negate(`%in%`)

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

round(c(GIC(test0),GIC(test1),GIC(test2),GIC(test3)),2)
round(c(GBIC(test0),GBIC(test1),GBIC(test2),GBIC(test3)),2)
round(c(test0$loglik,test1$loglik,test2$loglik,test3$loglik),2)
round(cbind(test0$b$mu,test1$b$mu,test2$b$mu,test3$b$mu),4)
round(cbind(test0$b$sigma,test1$b$sigma,test2$b$sigma,test3$b$sigma),4)


# Grid search on Model 1
# ~~~~~~~~~~~~~~~~~~~~~~

lamt <- seq(0,0.6,length.out = 6)
at <- seq(1.5,3.5,length.out = 5)

testList <- NULL
for(i in 1:length(lamt)){
 if(i == 1){
 testList[[i]] <- splvm.fit(trPOL,famd,formd1,
                  control = list(method = "PEM", constraint = l11., full.hess = F, start.val = test1$b,
                  ghQqp = 30, iter.lim = 1e3, tol = sqrt(.Machine$double.eps), silent = F,information = "Fisher",
                  pml.control = list(type = "mcp", w.alasso = test1$b, lambda = lamt[i], pen.load = T))) } else {
 testList[[i]] <- splvm.fit(trPOL,famd,formd1,
                  control = list(method = "PEM", constraint = l11., full.hess = F, start.val = testList[[i-1]]$b,
                  ghQqp = 30, iter.lim = 1e3, tol = sqrt(.Machine$double.eps), silent = F,information = "Fisher",
                  pml.control = list(type = "mcp", w.alasso = test1$b, lambda = lamt[i], pen.load = T)))                    
  }
}

par <- t(sapply(1:length(lamt), function(r) lb2cb(testList[[r]]$b)))
sel1 <- c(1:7)*4
sel2 <- sel1 - 2
sel <- c(sel1,sel2)
plot(lamt,par[,sel2[1]], xlab = expression(lambda), ylab = expression(alpha), ylim = c(min(par[,sel]),max(par[,sel])),
     type = "l", main = "Trajectories for location & scale penalised loadings (slopes)", col = "royalblue")
for(i in sel2){ lines(lamt,par[,i],col = "royalblue", lwd = 1) };
for(i in sel1){ lines(lamt,par[,i],col = "orange", lwd = 1) };
abline(h = 0, col = "red", lwd = 2, lty = 2)

GBICList <- LLList <- NULL
for(i in 1:length(lamt)){
  GBICList[[i]] <- GBIC(testList[[i]])
  LLList[[i]] <- testList[[i]]$loglik
}

plot(lamt,unlist(GBICList), col = "black",
     main = "GBIC: Result for different Tuning Parameter Values", xlab = expression(lambda), ylab = "GBIC")

plot(lamt,unlist(LLList), col = "black",
     main = "Log-likelihood: Result for different Tuning Parameter Values", xlab = expression(lambda), ylab = "Log-likelihood")


# Grid search on Model 3
# ~~~~~~~~~~~~~~~~~~~~~~

lamt <- seq(0,0.5,by = 0.05)
at <- c(2.5,3.7,4.5)

testList3 <- vector(mode = "list",length = length(at))
for(j in 1:length(at)){
for(i in 1:length(lamt)){
 if(i == 1){
  testList3[[j]][[i]] <- splvm.fit(trPOL,famd,formd3,
                   control = list(method = "PEM", constraint = l31., full.hess = F, start.val = test3$b,
                   ghQqp = 30, iter.lim = 1e3, tol = sqrt(.Machine$double.eps), silent = F,information = "Fisher",
                   pml.control = list(type = "scad", w.alasso = test3$b, lambda = lamt[i], pen.load = T,a = at[j]))) } else {
  testList3[[j]][[i]] <- splvm.fit(trPOL,famd,formd3,
                   control = list(method = "PEM", constraint = l31., full.hess = F, start.val = testList3[[j]][[i-1]]$b,
                   ghQqp = 30, iter.lim = 1e3, tol = sqrt(.Machine$double.eps), silent = F,information = "Fisher",
                   pml.control = list(type = "scad", w.alasso = test3$b, lambda = lamt[i], pen.load = T,a = at[j])))                    
 }
}
}

# save.image("C:/Users/CARDENCA/Dropbox/Camilo and Irini/Research/GitHub/SPLVM/R/fz_exp1.RData")

par3 <- NULL
GBICList3 <- LLList3 <- matrix(NA,ncol=length(at),nrow = length(lamt))
colnames(GBICList3) <- colnames(LLList3) <- at
rownames(GBICList3) <- rownames(LLList3) <- lamt

for(j in 1:length(at)){
par3[[j]] <- t(sapply(1:length(lamt), function(r) lb2cb(testList3[[j]][[r]]$b)))
for(i in 1:length(lamt)){
  GBICList3[i,j] <- GBIC(testList3[[j]][[i]])
  LLList3[i,j] <- testList3[[j]][[i]]$loglik
}}

sel1 <- c(1:7)*5
sel2 <- sel1 - 2
sel3 <- sel2 - 1
sel <- sort(c(sel1,sel2,sel3))

j = 1
plot(lamt,par3[[j]][,sel1[1]], xlab = expression(lambda), ylab = expression(alpha), ylim = c(min(par3[[j]][,sel]),max(par3[[j]][,sel])),
     type = "l",col = "royalblue", lwd = 2) #  main = "Trajectories for location & scale penalised loadings (slopes)", 
for(i in sel1){ lines(lamt,par3[[j]][,i],col = "royalblue", lwd = 2) };
for(i in sel2){ lines(lamt,par3[[j]][,i],col = "orange", lwd = 2) };
for(i in sel3){ lines(lamt,par3[[j]][,i],col = "black", lwd = 2, lty = 1 ) };
abline(h = 0, col = "red", lwd = 2, lty = 2)
abline(v = 0.05, col = "red", lwd = 2, lty = 2)
legend(0.375,0.9,legend=c(expression(alpha[" i1,"~mu]), expression(alpha[" i2,"~mu]), expression(alpha[" i1,"~sigma])),
       col=c("black", "orange","royalblue"),lty = 1, cex = 1, inset = 0.02, box.col = "white")

plot(lamt,unlist(GBICList3[,j]), col = "black",
     main = "GBIC: Result for different Tuning Parameter Values", xlab = expression(lambda), ylab = "GBIC")
abline(h = GBICList3[2,j], col = "red", lwd = 1, lty = 2)
abline(v = 0.05, col = "red", lwd = 1, lty = 2)

plot(lamt,unlist(LLList3[,j]), col = "black",
     main = "Log-likelihood: Result for different Tuning Parameter Values", xlab = expression(lambda), ylab = "Log-likelihood")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# Example No.2 : Rizopoulos & Moustaki (2008) : WIRS #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

rm(list = ls())
set.seed(1234)
# rm(.Random.seed, envir=globalenv())

source("fa_main.R")
source("fb_auxf.R")
source("fc_mfit.R")
source("fd_mghq.R")
source("ff_msim.R")

library(xtable)
library(ltm)
test <- WIRS
rmltm <- ltm::ltm(test ~ z1*z2, IRT.param = F,control = list(GHk = 40, iter.em = 350),constraint = cbind(1:6,4,-1.64))
# constraint = cbind(1:6,4,-1.64)
# rmltm0 <- ltm::ltm(test ~ z1+z2, IRT.param = F,control = list(GHk = 31, iter.em = 350))
# coef(rmltm)
# coef(rmltm0)
# rmltm$log.Lik
# rmltm0$log.Lik
mod <- NULL
mod$b <- ex.lb
mod$Y <- ex.test.pen[[1]][[1]]$Y
mod$ghQ <- ourmod2$ghQ#ex.test.pen[[1]][[1]]$ghQ
mod$fam <- ex.test.pen[[1]][[1]]$fam
mod$pml.control <- ex.test.pen[[1]][[1]]$pml.control
mod$uploglik <- llkf(lb2cb(mod$b),mod$Y,mod$ghQ,mod$b,mod$fam) # Difference with penfa; c(alasso_fit@Optim$logl.unpen)
mod$info = "Hessian"
mod$hessian <- sche.test(lb2cb(mod$b),mod)$hessian
GBIC(mod)

ex.lb <- mb2lb(matrix(coef(rmltm),nrow = ncol(test),byrow = F),ex.test.unp$b)
# ex.lb0 <- NULL; ex.lb0$mu <- ex.lb$mu[,-4]

fam.ex <- rep("binomial",ncol(test))

# ex.form0 <- list("mu" = "~ Z1+Z2")
# ex.test.unp0 <- splvm.fit(test,fam.ex,ex.form0,
#                control = list(method = "EM", full.hess = F, start.val = ex.lb0,
#                ghQqp = 15, iter.lim = 2e3, tol = .Machine$double.eps^0.5, silent = F,information = "Fisher"))

ex.form <- list("mu" = "~ Z1*Z2")
ex.test.unp <- splvm.fit(test,fam.ex,ex.form,
               control = list(method = "EM", full.hess = F, start.val = ex.lb,
               ghQqp = 15, iter.lim = 2e3, tol = .Machine$double.eps^0.5, silent = F,information = "Fisher"))

round(cbind(ex.lb$mu,ex.test.unp$b$mu),3)
c(rmltm$log.Lik,ex.test.unp$loglik)
GBIC(ex.test.unp)

lamt <- seq(0,1,by = 0.1)
at <- seq(0.5,1.5,by = 0.5)

ex.test.pen <- vector(mode = "list",length = length(at))
for(j in 1:length(at)){
 for(i in 1:length(lamt)){
 if(i == 1){
ex.test.pen[[j]][[i]] <- splvm.fit(test,fam.ex,ex.form,
                         control = list(method = "PEM", full.hess = F,start.val = ex.test.unp$b,
                         ghQqp = 15, iter.lim = 2e3, tol = .Machine$double.eps^0.5, silent = F,information = "Fisher",
                         pml.control = list(type = "mcp", w.alasso = ex.test.unp$b, lambda = lamt[i], pen.load = T,
                         a = at[j]))) } else {
ex.test.pen[[j]][[i]] <- splvm.fit(test,fam.ex,ex.form,
                         control = list(method = "PEM", full.hess = F,start.val = ex.test.pen[[j]][[i-1]]$b,
                         ghQqp = 20, iter.lim = 2e3, tol = .Machine$double.eps^0.5, silent = F,information = "Fisher",
                         pml.control = list(type = "mcp", w.alasso = ex.test.unp$b, lambda = lamt[i], pen.load = T,
                         a = at[j]))) }
print(c(j,i))
 }
}
save.image("C:/Users/carde/Dropbox/Camilo and Irini/Research/GitHub/SPLVM/R/fz_exp2.RData")

j = 1

par2 <- NULL
GBICList2 <- LLList2 <- matrix(NA,ncol=length(at),nrow = length(lamt))
colnames(GBICList2) <- colnames(LLList2) <- at
rownames(GBICList2) <- rownames(LLList2) <- lamt

for(j in 1:length(at)){
par2[[j]] <- t(sapply(1:length(lamt), function(r) lb2cb(ex.test.pen[[j]][[r]]$b)))
for(i in 1:length(lamt)){
  GBICList2[i,j] <- GBIC(ex.test.pen[[j]][[i]])
  LLList2[i,j] <- ex.test.pen[[j]][[i]]$loglik
}}

GBICList2[which(GBICList2 < GBIC(ex.test.unp), arr.ind = TRUE)]
LLList2[which(GBICList2 == min(GBICList2), arr.ind = TRUE)]
GBIC(ex.test.unp)
ex.test.unp$loglik

sel1 <- c(1:6)*4
sel2 <- sel1 - 1
sel3 <- sel2 - 1
sel <- sort(c(sel1,sel2,sel3))

j = 2
plot(lamt,par2[[j]][,sel1[1]], xlab = expression(lambda), ylab = expression(alpha), ylim = c(min(par2[[j]][,sel]),max(par2[[j]][,sel])),
     type = "l",col = "royalblue", lwd = 2) #  main = "Trajectories for location & scale penalised loadings (slopes)", 
for(i in sel1){ lines(lamt,par2[[j]][,i],col = "royalblue", lwd = 2) };
for(i in sel2){ lines(lamt,par2[[j]][,i],col = "orange", lwd = 2) };
for(i in sel3){ lines(lamt,par2[[j]][,i],col = "black", lwd = 2, lty = 1 ) };
abline(h = 0, col = "red", lwd = 2, lty = 2)
abline(v = lamt[7], col = "red", lwd = 2, lty = 2)
legend(0.0,-5,legend=c(expression(alpha[" i1,"~pi]), expression(alpha[" i2,"~pi]), expression(alpha[" i3,"~pi])),
       col=c("black", "orange","royalblue"),lty = 1, cex = 1, inset = 0.02, box.col = "white")

# xtable(ex.test.unp$b$mu,digits = 3)
# xtable(ex.test.pen[[2]][[7]]$b$mu,digits = 3)

round(cbind(ex.test.unp$b$mu,ex.lb$mu,ex.test.pen[[j]]$b$mu),3)
ex.test.unp$loglik; names(ex.test.pen[[j]])

GBIC(ex.test.pen[[j]])
GBIC(ex.test.unp)

ourmod <- ex.test.pen[[2]][[7]]
round(cbind(ourmod2$b$mu,mod$b$mu),2)
splvm.plot(2, ourmod2, control.plot = list(plot.3D = T, bw = F))
splvm.plot(2, mod, control.plot = list(plot.3D = T, bw = F))

