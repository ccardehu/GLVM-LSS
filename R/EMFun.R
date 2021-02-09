rm(list = ls())
#rm(list= ls()[!(ls() %in% c("ex1", "simR"))])
set.seed(1234)
# rm(.Random.seed, envir=globalenv())

library(mvtnorm)
library(MASS)
library(foreach)
library(doSNOW)
library(parallel)
library(itertools)
#library(profvis)
#library(ltm)
library(xtable)
# library(gamlss.dist)
library(numDeriv)

`%!in%` <- Negate(`%in%`)

source("ddFun.R")
source("SimFA.R")
source("GHFun.R")
source("scFun.R")
source("GLVMfit.R")
source("graphFun.R")

n = 1000     # Number of individuals
p = 10       # Number of items
nsim = 1000  # Number of simulations
form <- list("mu" = "~ Z1", "sigma" = "~ Z1")
#form1 <- list("mu" = "~ Z1 + I(Z1^2)", "sigma" = "~ 1") #  , "sigma" = "~ Z1"
#form2 <- list("mu" = "~ Z1", "sigma" = "~ 1") #  , "sigma" = "~ Z1"
fam <- rep("normal",p)# c(rep("normal", p/2), rep("ZIpoisson", p/2))

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

# l1 when model misspecification
# ______________________________
#
#l11 <- l1
#l11$sigma <- l11$sigma[,-c(2:ncol(l11$sigma)), drop = F]
#l12 <- l11
#l12$mu <- l12$mu[,-ncol(l12$mu), drop = F]

lc <- NULL
lc$mu <- matrix(runif(length(l1$mu), min = 0.1, max = 0.5),nrow = p)
#lc$mu[,1] <-  runif(p,1,2)
#lc$mu[,2] <-  runif(p,-0.5,0.5)
#lc$mu[,3] <-  runif(p,-0.5,0.5)
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

ex1$b$mu - borg$mu
ex1$b$sigma - borg$sigma


plotGLVM(item = 1, mod = ex1, morg = simR, plot.org = F,
         plot.mean = F, plot.sd = F, quant = c(0.025,0.25,0.75,0.975),
         sep.plots = F, plot.3D = T, plot.dist = T, plot.addpoints = T)

plot.score(ex1)

ex1$b$mu; borg$mu
ex1$b$sigma; borg$sigma
#plot(ex1$cvgRes[,sample(1:ncol(ex1$cvgRes),2,replace = F)])

EMfun <- function(l, silent = T){ # CHANGE HERE!
        
        AA <- simGLVM(n = n, p = p, form = form, dist = fam, coefs = lc, loadmt = l1)
        tmp1 <- GLVM.fit(Y = AA$Y, fam = fam, form = AA$form, silent = T, ghp = 50, iter.lim = 700, tol = 1e-7, loadmt = l1, icoefs = lc, useoptim = F, skipEM = F)
        if(silent == F) cat("Run = ", l, "Iter. to Convg.: ",tmp1$iter,"\n")
        return(c(unlist(tmp1$b),tmp1$iter))
        
}

cores <- detectCores()-2
cl<-makeCluster(cores)
registerDoSNOW(cl)

progress <- function(n) setTxtProgressBar(txtProgressBar(max = nsim, style = 3), n)
opts <- list(progress = progress)

FCOL <- foreach(l = 1:nsim,
                .combine = rbind,
                .packages = c("mvtnorm", "statmod"),
                .options.snow = opts) %dopar% EMfun(l, silent = T)

stopCluster(cl)

save.image(paste0("nsim",nsim,"_n",n, "_p", p,"_ZIP(LinMu_QtSig).RData"))

# FCOLorg <- FCOL # FCOL <- FCOLorg
FCOL1 <- FCOL[FCOL[,ncol(FCOL)] %!in% c(700,-999),-ncol(FCOL)]
#FCOLHT <- FCOL1[,1:40]
#FCOLHM <- FCOL1[,51:90]
lc1 <- NULL
lc1$mu <- lc$mu*l1$mu
lc1$sigma <- lc$sigma*l1$sigma

coefmod(colMeans(FCOL1),borg,gr = ex1$gr,loadmt = l1)$mu; borg$mu
coefmod(colMeans(FCOL1),borg,gr = ex1$gr,loadmt = l1)$sigma; borg$sigma

lcmat <- matrix(unlist(lc), nrow = nrow(FCOL1), ncol = ncol(FCOL1), byrow = T)

bias  <- FCOL1 - lcmat
rbias <- (bias/lcmat)*100
rmse  <- sqrt(colMeans(bias^2))
Bias  <- colMeans(bias)
RBias <- colMeans(rbias)

boxplot(bias[,01:10]); abline(h= 0, col = 2, lwd = 2)
boxplot(bias[,11:20]); abline(h= 0, col = 2, lwd = 2)
boxplot(bias[,21:30]); abline(h= 0, col = 2, lwd = 2)
boxplot(bias[,31:40]); abline(h= 0, col = 2, lwd = 2)
boxplot(bias[,41:50]); abline(h= 0, col = 2, lwd = 2)
boxplot(bias[,51:60]); abline(h= 0, col = 2, lwd = 2)
boxplot(bias[,61:70]); abline(h= 0, col = 2, lwd = 2)
boxplot(bias[,71:80]); abline(h= 0, col = 2, lwd = 2)


length(FCOL[FCOL[,ncol(FCOL)] != 700 & FCOL[,ncol(FCOL)] != -999, ncol(FCOL)])
summary(FCOL[FCOL[,ncol(FCOL)] != 700 & FCOL[,ncol(FCOL)] != -999, ncol(FCOL)])

# For horizontal table in Appendix
# ________________________________
Mrmse <- matrix(rmse, nrow = p, byrow = F)
Mbias <- matrix(Bias, nrow = p, byrow = F)
Mrbias <- matrix(RBias, nrow = p, byrow = F)
MR <- matrix(NA, nrow = p, ncol = ncol(Mbias)*3)
for(i in 1:ncol(Mrmse)){MR[,((i-1)*3+1):((i)*3)] <- round(cbind(Mrmse[,i], Mbias[,i], Mrbias[,i]),3)}
rownames(MR) <- paste("Item",1:p); MR
#xtable(MR,digits = 3)

# For vertical table in Appendix
# ______________________________
Morg <- matrix(unlist(lc)); rownames(Morg) <- names(unlist(lc)); colnames(Morg) <- "True"
Mest <- matrix(colMeans(FCOL1)); colnames(Mest) <- "Avg. Est."
Mbias1 <- matrix(Bias); colnames(Mbias1) <- "Bias"
Mrbias1 <- matrix(RBias); colnames(Mrbias1) <- "RB (%)"
Mrmse1 <- matrix(rmse); colnames(Mrmse1) <- "RMSE"
Msd <- matrix(apply(FCOL1,2,sd)); colnames(Msd) <- "SD"
MR <- round(cbind(Morg, Mest, Msd, Mrmse1),3); head(MR,15)
MR <- round(cbind(Morg, Mbias1, Mrbias1, Mrmse1),3); head(MR,15)
#xtable(unname(MR),digits = 3)

